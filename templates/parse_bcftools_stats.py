#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--stats", default="${stats_file}", help="")
parser.add_argument("--dataset", default="${dataset}", help="")
parser.add_argument("--step", default="${step}", help="")
parser.add_argument("--csv_out", default="${csv_out}", help="")

args = parser.parse_args()


def parse_bcftools_stats(stats, dataset, step, csv_out):
    """
    :param :
    :return:
    """

    for sn in open(stats).readlines():
        if sn.startswith('SN'):
            sn = sn.strip().split('\t')
            if 'number of records' in sn[2]:
                records = sn[3]
            elif 'number of no-ALTs' in sn[2]:
                no_alts = sn[3]
            elif 'number of SNPs' in sn[2]:
                snps = sn[3]
            elif 'number of indels' in sn[2]:
                indels = sn[3]
            elif 'number of others' in sn[2]:
                others = sn[3]
            elif 'number of multiallelic sites' in sn[2]:
                multiallelics = sn[3]
            elif 'Sample size' in sn[2]:
                samples = sn[3]
                

    data = [ dataset, step, samples, records, no_alts, snps, indels, others, multiallelics ]
    
    out = open(csv_out, 'w')
    out.writelines('\\t'.join(
        ['Dataset', 'Step', 'Number of samples', 'Number of records', 'Number of no-ALTs', 'Number of SNPs', 'Number of indels', 'Number of others', 'Number of multiallelic sites'])+'\\n')
    out.writelines('\\t'.join(data)+'\\n')
    out.close()

if __name__ == '__main__':
    parse_bcftools_stats(args.stats, args.dataset, args.step, args.csv_out)

