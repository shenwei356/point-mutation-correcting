#!/usr/bin/env python
# https://github.com/shenwei356/bio_scripts
# Author     : Wei Shen
# Contact    : shenwei356@gmail.com
# LastUpdate : 2015-06-29

import argparse
import sys
from collections import Counter, defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description='Point mutation correcting by k-mer clustering',
                                     epilog='https://github.com/shenwei356/point-mutation-correcting')

    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin,
                        help=' '.join(['String file. One record per line, case ignored. Default: STDIN.',
                                       'Preset counts could be given following string, separated by "\\t".']))
    parser.add_argument('-m', '--max-mutation-sites', type=int, default=1,
                        help='Maximum mutation sites [1]')
    parser.add_argument('-k', '--kmer', type=int, default=9,
                        help='K-mer length for clustering.  Recommending 3 <= k <= 0.5 * length of string [9].')
    parser.add_argument("-s", "--sort", action="store_true",
                        help='Sorting strings in result')
    args = parser.parse_args()

    if args.max_mutation_sites < 0:
        sys.stderr.write("value of option --max-mutation-sites should >= 0")
        sys.exit(1)

    return args


def compute_kmers(s, k):
    kmers = set()
    l = len(s)
    end = l - k if l - k > 0 else 0
    for i in range(0, end + 1):
        kmers.add(s[i:i + k])
    return kmers


def is_similar_key(s1, s2, max_mutation_sites=0):
    if len(s1) != len(s2):
        sys.stderr.write('unequal length of strings: "{}" and "{}"\n'.format(s1, s2))
        sys.exit(1)

    d = 0
    for i in range(0, len(s1)):
        if s1[i:i + 1] != s2[i:i + 1]:
            d += 1
            if d > max_mutation_sites:
                return False, d
    return True, d

def counting(fh, debug=False):
    counts = Counter()
    cnt = 0
    for line in fh:
        if line.isspace() or line[0] == '#':
            continue

        # tick
        cnt += 1
        if debug and cnt % 10000 == 0:
            sys.stderr.write('Counting records: {}\r'.format(cnt))

        # key and preset value
        items = line.rstrip().split('\t')
        if len(items) > 1:
            key, size = items[0].upper(), int(items[1])
        else:
            key, size = items[0].upper(), 1

        counts[key] += size

    return counts

def correct_point_mutation(counts, max_mutation_sites=0, k_len=7, debug=False):
    clusters = defaultdict(Counter)
    similar_keys = dict()  # a cache: key: similar_key.
    kmers_map = defaultdict(set)  # kmer: [key1, key2]

    cnt = 0
    for key in sorted(counts.keys(), key=lambda k: -counts[k]):
        # tick
        cnt += 1
        if debug and cnt % 10000 == 0:
            sys.stderr.write('Checking records: {}\r'.format(cnt))

        if max_mutation_sites == 0:  # do not need the complex computation below
            clusters[key][key] += size
            continue

        size = counts[key]

        # all kmers
        kmers = compute_kmers(key, k_len)

        if key in clusters or len(clusters) == 0:  # existed
            clusters[key][key] += size
        elif key in similar_keys:  # yes, we have known that key is similar to similar_keys[key]
            clusters[similar_keys[key]][key] += size
        else:
            candidates = dict()
            tested = set()
            for kmer in kmers:  # all kmers
                for k in kmers_map[kmer]:  # all keys in a kmer
                    if k == key:  # do not compare with itself
                        continue

                    if k in tested:  # multiple kmers may refer to same keys
                        continue
                    tested.add(k)

                    similar, d = is_similar_key(key, k, max_mutation_sites)
                    if similar:
                        candidates[k] = d
                        break

            if len(candidates) == 0:  # no candidates
                clusters[key][key] += size
            else:  # find the key with minimum d
                k = sorted(candidates.keys(), key=lambda x: candidates[x], reverse=True)[0]
                if k in clusters:
                    clusters[k][key] += size
                    similar_keys[key] = k
                elif k in similar_keys:  # Attention!!
                    clusters[similar_keys[k]][key] += size
                    similar_keys[key] = similar_keys[k]
                else:
                    clusters[k][key] += size
                    similar_keys[key] = k

        for kmer in kmers:
            kmers_map[kmer].add(key)

    return clusters


if __name__ == '__main__':
    args = parse_args()

    counts = counting(args.infile, debug=True)
    clusters = correct_point_mutation(counts, max_mutation_sites=args.max_mutation_sites, k_len=args.kmer, debug=True)

    sys.stderr.write('Original unique strings: {}. After correcting: {}\n'.format(len(counts), len(clusters)))

    keys = clusters.keys()
    if args.sort:
        keys = sorted(keys)
    for key in keys:
        cluster = clusters[key]
        sys.stdout.write('{}\t{}\t{}\n'.format(key, sum(cluster.values()), cluster))
