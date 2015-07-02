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
                        help=' '.join(['K-mer length for clustering.  Recommending 3 <= k <= L / (1+M),',
                                       'L for length of string and M for maximum mutation sites. [9]']))
    parser.add_argument('-s', '--sort', action='store_true',
                        help='Sorting strings in result')
    parser.add_argument("-v", "--verbose", help='Verbosely print information',
                        action="count", default=0)
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


def counting(fh, k_len, verbose=0):
    counts, kmer_counts = Counter(), Counter()
    kmers_map = defaultdict(set)  # kmer: [key1, key2]
    cnt = 0
    for line in fh:
        if line.isspace() or line[0] == '#':
            continue

        # tick
        cnt += 1
        if verbose > 0 and cnt % 10000 == 0:
            sys.stderr.write('Counting records: {}\r'.format(cnt))

        # key and preset value
        items = line.rstrip().split('\t')
        if len(items) > 1:
            key, size = items[0].upper(), int(items[1])
        else:
            key, size = items[0].upper(), 1

        counts[key] += size

        for kmer in compute_kmers(key, k_len):
            # kmer_counts[kmer] += 1
            kmers_map[kmer].add(key)

    return counts, kmers_map, kmer_counts


def correct_point_mutation(counts, kmers_map, kmer_counts, max_mutation_sites=0, k_len=7, verbose=0):
    clusters = defaultdict(Counter)
    similar_keys = dict()  # a cache: key: similar_key.
    stats = Counter()

    cnt = 0
    for key in sorted(counts.keys(), key=lambda k: -counts[k]):
        # tick
        cnt += 1
        if verbose > 0 and cnt % 10000 == 0:
            sys.stderr.write('Checking records: {}\r'.format(cnt))

        size = counts[key]

        if max_mutation_sites == 0:  # do not need the complex computation below
            clusters[key][key] += size
            continue

        if size == 1:
            stats['single'] += 1

        kmers = compute_kmers(key, k_len)
        if key in clusters or len(clusters) == 0:  # existed
            if verbose > 1:
                if len(clusters) == 0:
                    sys.stderr.write('{}: first one\n'.format(key))
                else:
                    sys.stderr.write('{}: already in clusters\n'.format(key))
            clusters[key][key] += size
        elif key in similar_keys:  # yes, we have known that key is similar to similar_keys[key]
            # find the original key
            kk = key
            while True:
                kk = similar_keys[kk]
                if kk not in similar_keys:
                    break
            clusters[kk][key] += size
            if verbose > 1:
                sys.stderr.write('{}: known to be similar to {}\n'.format(key, similar_keys[key]))
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

            if len(candidates) == 0:  # no candidates
                if verbose > 1:
                    sys.stderr.write('{}: no candidates, create new cluster\n'.format(key))
                clusters[key][key] += size
            else:
                cached = False
                ks = sorted(candidates.keys(), key=lambda x: candidates[x])
                for k in ks:
                    if verbose > 1:
                        sys.stderr.write('{}: candi: {}\n'.format(key, k))
                    if k in clusters:
                        clusters[k][key] += size
                        similar_keys[key] = k
                        cached = True
                        if verbose > 1:
                            sys.stderr.write('{}: candi: {} in cluster\n'.format(key, k))
                        break
                    elif k in similar_keys:
                        # find the original key
                        kk = k
                        while True:
                            kk = similar_keys[kk]
                            if kk not in similar_keys:
                                break
                        clusters[kk][key] += size
                        similar_keys[key] = k
                        cached = True
                        if verbose > 1:
                            sys.stderr.write('{}: candi: {} has similar key {}\n'.format(key, k, similar_keys[k]))
                        break
                if not cached:
                    k = ks[0]
                    clusters[k][key] += size
                    similar_keys[key] = k
                    if verbose > 1:
                        sys.stderr.write('{}: candi: {} new cluster\n'.format(key, k))
                if verbose > 1:
                    sys.stderr.write('\n')
    sys.stderr.write('summary: {}\n'.format(stats))
    return clusters


if __name__ == '__main__':
    args = parse_args()

    counts, kmers_map, kmer_counts = counting(args.infile, args.kmer, verbose=args.verbose)
    clusters = correct_point_mutation(counts, kmers_map, kmer_counts, max_mutation_sites=args.max_mutation_sites,
                                      k_len=args.kmer, verbose=args.verbose)

    sys.stderr.write('Original unique strings: {}. After correcting: {}\n'.format(len(counts), len(clusters)))

    keys = clusters.keys()
    if args.sort:
        keys = sorted(keys)
    for key in keys:
        cluster = clusters[key]
        key = sorted(cluster.keys(), key=lambda k: -cluster[k])[0]
        sys.stdout.write('{}\t{}\t{}\n'.format(key, sum(cluster.values()), cluster))
