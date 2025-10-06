#!/usr/bin/env python3
import argparse
import sys
import os
import collections


def get_para(options):
    description = '''Read the summary.txt files from MitoZ and output a summary
    across samples.
    '''

    parser = argparse.ArgumentParser(description=description,
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-i", dest='summary_files', metavar="<file>", required=True,
        help="A file containing paths to summary.txt files. Per-line format: sample_code path [required]")

    parser.add_argument("-o", dest='outfile', metavar="<file>", default="out.txt",
        help="outfile name [%(default)s]")

    parser.add_argument("-t", dest='stat_tRNA', action='store_true',
        help="Also do statistics for tRNA genes. By default, also for PCGs and rRNA genes. [%(default)s]")

    if len(options) == 0:
        parser.print_help()
        sys.exit(0)

    return parser.parse_args(options)


def extract_gene_stat(summary_file=None, gene_set=None):

    #Seq_id        Length(bp)     Circularity    Closely_related_species
    seq_title = '#Seq_id        Length(bp)     Circularity'
    seq_title_on = False
    seq_lens = []
    seq_topologies = []

    gene_title = '#Seq_id        Start  End'
    gene_title_on = False
    gene_freq_counter = collections.OrderedDict()
    for k in gene_set:
        gene_freq_counter[k] = 0

    gene_num_stat_title = 'Protein coding genes totally found'
    gene_num_stat_on = False
    gene_num_stat_dict = {}

    with open(summary_file, 'r') as fh:
        for i in fh:
            i = i.strip()
            if not i or i.startswith('-----'):
                continue
            if seq_title in i:
                seq_title_on = True
                continue
            if seq_title_on and gene_title in i:
                gene_title_on = True
                seq_title_on = False
                continue

            if gene_num_stat_title in i:
                gene_title_on = False
                gene_num_stat_on = True

            if seq_title_on:
                line = i.split()            
                seq_lens.append(line[1])
                seq_topologies.append(line[2])
                continue

            if gene_title_on:
                line = i.split()
                gene_type, gene_name = line[5:7]
                freq = line[-1]
                if gene_name in gene_set:
                    gene_freq_counter[gene_name] += 1
                continue

            if gene_num_stat_on:
                gene_freq_counter['PCG_found_num'] = i.split()[-1]
                gene_freq_counter['tRNA_found_num']  = fh.readline().strip().split()[-1]
                gene_freq_counter['rRNA_found_num']  = fh.readline().strip().split()[-1]
                break


    return gene_freq_counter



def main(options):
    args = get_para(options)

    gene_set = ["PCG_found_num", "tRNA_found_num", "rRNA_found_num", "ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "s-rRNA", "l-rRNA"]

    tRNA_set = ["trnA(ugc)", "trnC(gca)", "trnD(guc)", "trnE(uuc)", "trnF(gaa)", "trnG(ucc)", "trnH(gug)", "trnI(gau)", "trnK(uuu)", "trnL(uaa)", "trnL(uag)", "trnM(cau)", "trnN(guu)", "trnP(ugg)", "trnQ(uug)", "trnR(ucg)", "trnS(gcu)", "trnS(uga)", "trnT(ugu)", "trnV(uac)", "trnW(uca)", "trnY(gua)"]
    
    if args.stat_tRNA:
        gene_set.extend(tRNA_set)

    with open(args.summary_files, 'r') as fh, open(args.outfile, 'w') as fhout:
        header = ["Sample", *gene_set]
        print("\t".join(header), file=fhout)
        for i in fh:
            i = i.strip()
            if not i:
                continue
            sample, summary_file = i.split()[0:2]
            gene_freq_counter = extract_gene_stat(
                summary_file=summary_file,
                gene_set=gene_set)

            output = [str(gene_freq_counter[g]) for g in gene_set]
            content = [sample, *output]
            print("\t".join(content), file=fhout)


if __name__ == '__main__':
    main(sys.argv[1:])














