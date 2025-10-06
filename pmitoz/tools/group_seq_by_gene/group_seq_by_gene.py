#!/usr/bin/env python3
import sys
import argparse
from Bio import SeqIO
import os
import re


def add_arguments(parser):
    parser.add_argument('-r', dest='file_info', metavar='<file>', help='the gene file list. Per-line format: Abbreviation geneFilePath. The abbreviation will be added to the seqid to indicate different samples.')

    parser.add_argument('-d', dest='delimiter', metavar='<str>', default=';',
        help='the delimiter between the abbreviation and the seqid [%(default)s]')

    parser.add_argument('-p', dest='prefix', metavar='<str>', default='MitoZ',
        help='the prefix of all result files [%(default)s]')

    parser.add_argument('-clean_header', action='store_true',
        help="Only shows the 'Abbreviation' in the sequence header [%(default)s]")

    return parser


def get_sample_prefix(file_info=None):
    '''
abbreviation filePath
    '''
    file_abbr = {}
    with open(file_info, 'r') as fh:
        for i in fh:
            i = i.strip()
            if not i:
                continue

            abbr, f = i.split()[0:2]
            file_abbr[f] = abbr

    return file_abbr



def patch_gene(file_info=None, prefix='MitoZ', file_abbr=None, delimiter=';', clean_header=False):
    gene_count = {}
    with open(file_info, 'r') as fh:
        for i in fh:
            i = i.strip()
            if not i:
                continue

            abbr, f = i.split()[0:2]
            for rec in SeqIO.parse(f, 'fasta'):
                gene = rec.description.split(';')[1]
                gene_file = prefix + '.gene-' + gene + '.fa'
                if gene not in gene_count:
                    if os.path.exists(gene_file):
                        sys.exit('Please remove the {prefix}.gene-*.fa files first!'.format(prefix=prefix))

                gene_count[gene] = 1
                fhout = open(gene_file, 'a')
                rec.description = re.sub(r'^.*?\;', file_abbr[f]+delimiter, rec.description)
                if clean_header:
                    print('>'+file_abbr[f]+'\n'+str(rec.seq), file=fhout)
                else:
                    print('>'+rec.description+'\n'+str(rec.seq), file=fhout)
                fhout.close()


def main(args):
    file_abbr = get_sample_prefix(file_info=args.file_info)
    patch_gene(
        file_info=args.file_info,
        prefix=args.prefix,
        file_abbr=file_abbr,
        delimiter=args.delimiter,
        clean_header=args.clean_header)


if __name__ == '__main__':
    main()