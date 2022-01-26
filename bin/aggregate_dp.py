#!/usr/bin/env python

import argparse
import sys
from collections import Counter

argparser = argparse.ArgumentParser(description = 'Aggregates DP in autosomal chromosomes, chromosome X and Y. Reads input as a stream from `samtools depth` command.')
argparser.add_argument('-t', '--thresholds', metavar = 'number', dest = 'dp_thresholds', type = int, nargs = '+', default = [1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110], help = 'List of DP thresolds for counting coverage. Default: [1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110]')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_filename', required = True, help = 'Output file name prefix.')

if __name__ == '__main__':
    args = argparser.parse_args()
    chromosomes = set()
    chrom_sum_dp = Counter()
    chrom_n_bases = Counter()
    chrom_n_covered = Counter()
    chrom_n_by_dp = dict()
    for dp_threshold in args.dp_thresholds:
        chrom_n_by_dp[dp_threshold] = Counter()
    for line in sys.stdin:
        chrom, pos, dp = line.rstrip().split()
        chromosomes.add(chrom)
        dp = int(dp)
        chrom_sum_dp[chrom] += dp
        chrom_n_bases[chrom] += 1
        if dp > 0:
            chrom_n_covered[chrom] += 1
        for dp_threshold in args.dp_thresholds:
            if dp >= dp_threshold:
                chrom_n_by_dp[dp_threshold][chrom] += 1

    with open(args.out_filename + '.by_chrom.txt', 'wt') as ofile:
        ofile.write('CHROM\tN_BASES\tN_BASES_COVERED')
        for dp_threshold in args.dp_thresholds:
            ofile.write(f'\tN_BASES_DP{dp_threshold}')
        ofile.write('\tSUM_DP\tAVG_DP\n')
        for chrom, sum_dp in chrom_sum_dp.items():
            ofile.write(f'{chrom}\t{chrom_n_bases[chrom]}\t{chrom_n_covered[chrom]}')
            for dp_threshold in args.dp_thresholds:
                ofile.write(f'\t{chrom_n_by_dp[dp_threshold][chrom]}')
            ofile.write(f'\t{sum_dp}\t{sum_dp / chrom_n_bases[chrom]:.6f}\n')


    with open(args.out_filename + '.txt', 'wt') as ofile:
        if any(x.startswith('chr') for x in chromosomes):
            autosomal_chroms = [f'chr{x}' for x in range(1, 23)]
            sex_chroms = ['chrX', 'chrY']
        else:
            autosomal_chroms = [f'{x}' for x in range(1, 23)]
            sex_chroms = ['X', 'Y']
        all_chroms = autosomal_chroms + sex_chroms
        
        header = []
        values = []

        sum_dp = 0
        n_bases = 0
        n_covered = 0
        n_by_dp = dict((dp, 0) for dp in args.dp_thresholds)
        for chrom in chromosomes:
            sum_dp += chrom_sum_dp[chrom]
            n_bases += chrom_n_bases[chrom]
            n_covered += chrom_n_covered[chrom]
            for dp in args.dp_thresholds:
                n_by_dp[dp] += chrom_n_by_dp[dp][chrom]

        header.append('ALL_CHROMS')
        values.append(f'{len(chromosomes)}')
        header.append('ALL_N_BASES')
        values.append(f'{n_bases}')
        header.append('ALL_PRCNT_COVERED')
        values.append(f'{(n_covered / n_bases) * 100:.6f}')
        for dp in args.dp_thresholds:
            header.append(f'ALL_PRCNT_DP{dp}')
            values.append(f'{(n_by_dp[dp] / n_bases) * 100:.6f}')
        header.append('ALL_AVG_DP')
        values.append(f'{sum_dp / n_bases:.6f}')

        n_autosomals = 0
        auto_sum_dp = 0
        auto_n_bases = 0
        auto_n_covered = 0
        auto_n_by_dp = dict((dp, 0) for dp in args.dp_thresholds)
        for chrom in autosomal_chroms:
            if chrom_n_bases[chrom] > 0:
                n_autosomals += 1
            auto_sum_dp += chrom_sum_dp[chrom]
            auto_n_bases += chrom_n_bases[chrom]
            auto_n_covered += chrom_n_covered[chrom]
            for dp in args.dp_thresholds:
                auto_n_by_dp[dp] += chrom_n_by_dp[dp][chrom]
        
        header.append('AUTO_CHROMS')
        values.append(f'{n_autosomals}')
        header.append('AUTO_N_BASES')
        values.append(f'{auto_n_bases}')
        header.append('AUTO_PRCNT_COVERED')
        values.append(f'{(auto_n_covered / auto_n_bases) * 100:.6f}')
        for dp in args.dp_thresholds:
            header.append(f'AUTO_PRCNT_DP{dp}')
            values.append(f'{(auto_n_by_dp[dp] / auto_n_bases) * 100:.6f}')
        header.append('AUTO_AVG_DP')
        values.append(f'{auto_sum_dp / auto_n_bases:.6f}')

        for chrom, column in zip(sex_chroms, ['X', 'Y']):
            header.append(f'{column}')
            header.append(f'{column}_AVG_DP')
            header.append(f'{column}_NORM_AVG_DP')
            if chrom_n_bases[chrom] > 0:
                avg_dp = chrom_sum_dp[chrom] / chrom_n_bases[chrom]
                values.append('1')
                values.append(f'{avg_dp:.6f}')
                values.append(f'{avg_dp / (auto_sum_dp / auto_n_bases):.6f}')
            else:
                values.append('0')
                values.append('NA')
                values.append('NA')

        ofile.write('{}\n'.format('\t'.join(header)))
        ofile.write('{}\n'.format('\t'.join(values)))
