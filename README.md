Point mutation correcting
========

Point mutation correcting by k-mer clustering.

I used this tool to correct point mutations come from PCR or sequencing error in a de novo sequencing project.

Strings should be in equal length and one record per line, case ignored, can be 
from file or STDIN. Preset count could be given following string, separated by "\t".


Usage
-----

    usage: point_mutation_correcting.py [-h] [-m MAX_MUTATION_SITES] [-k KMER]
                                        [-s] [-v]
                                        [infile]
    
    Point mutation correcting by k-mer clustering
    
    positional arguments:
      infile                String file. One record per line, case ignored.
                            Default: STDIN. Preset counts could be given following
                            string, separated by "\t".
    
    optional arguments:
      -h, --help            show this help message and exit
      -m MAX_MUTATION_SITES, --max-mutation-sites MAX_MUTATION_SITES
                            Maximum mutation sites [1]
      -k KMER, --kmer KMER  K-mer length for clustering. Recommending 3 <= k <= L
                            / (1+M), L for length of string and M for maximum
                            mutation sites. [9]
      -s, --sort            Sorting strings in result
      -v, --verbose         Verbosely print information
    
    https://github.com/shenwei356/point-mutation-correcting


Example
-----



Input
   
    AATGGACGCGATTGCCACGG
    AATGGACGCGATTGCCACGG
    AATGGACGCGATTGtCACGG
    GAATACACACAACACCGTAT
    GAATACACACAACACCGTAT
    GAATACACACAACACCGTAT
    GAATACACACAACgCCGTAT
    GAATACACACAACACtGTAT
    GAATACACACgACACCGTAT
    
Command

    cat example | exampleexamplepoint_mutation_correcting.py -m 1 -k 9
    
Output

    AATGGACGCGATTGCCACGG    3       Counter({'AATGGACGCGATTGCCACGG': 2, 'AATGGACGCGATTGTCACGG': 1})
    GAATACACACAACACCGTAT    6       Counter({'GAATACACACAACACCGTAT': 3, 'GAATACACACAACACTGTAT': 1, 'GAATACACACAACGCCGTAT': 1, 'GAATACACACGACACCGTAT': 1})

Lincence
-------

Copyright (c) 2015, Wei Shen (shenwei356@gmail.com)


[MIT License](https://github.com/shenwei356/point-mutation-correcting/blob/master/LICENSE)