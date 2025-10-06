#!/usr/bin/env python3
"""
This file is part of MitoZ.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
PERFORMANCE OF THIS SOFTWARE. REVERSE ENGINEERING IS STRICTLY PROHIBITED.

COPYRIGHT © 2019-2022 Guanliang Meng. ALL RIGHTS RESERVED.

"""

import argparse
import sys
import os
import re
import subprocess
import time
from glob import glob



def prepare_common_parser():
    common_parser = argparse.ArgumentParser(add_help=False,
        formatter_class=argparse.RawTextHelpFormatter)

    common_group = common_parser.add_argument_group('Common arguments')

    common_group.add_argument("--outprefix", metavar="<STR>", required=True,
        help="output prefix")

    common_group.add_argument("--thread_number", metavar="<INT>", default="8",
        help="thread number [%(default)s]")

    return common_parser


def prepare_fastq_parser():
    fastq_parser = argparse.ArgumentParser(add_help=False,
        formatter_class=argparse.RawTextHelpFormatter)

    fastq_group = fastq_parser.add_argument_group('Input fastq files')

    fastq_group.add_argument("--fastq1", metavar="<STR>", default='', help="fastq 1 file")

    fastq_group.add_argument("--fastq2", metavar="<STR>", default='', help="fastq 2 file")

    fastq_group.add_argument("--fastq_quality_shift", default=False,
        action="store_true", help="the input is in the Illumina 1.3+" +\
        " FASTQ-like (Q+64) format. (default: Q+33)")

    fastq_group.add_argument("--fastq_read_length", metavar="<INT>", default='150',
        help='''read length of fastq reads, used by SOAPTrans and bwa.
    It must be >= 71 bp [%(default)s]''')

    return fastq_parser


def prepare_fastafile_parser():
    fastafile_parser = argparse.ArgumentParser(add_help=False,
        formatter_class=argparse.RawTextHelpFormatter)

    fastafile_group = fastafile_parser.add_argument_group('Input fasta file')

    fastafile_group.add_argument('--fastafile', metavar="<STR>", help="fasta file")

    return fastafile_parser


def prepare_filter_parser():
    filter_parser = argparse.ArgumentParser(add_help=False,
        formatter_class=argparse.RawTextHelpFormatter)

    filter_group = filter_parser.add_argument_group('Filter arguments')

    return filter_parser


def prepare_assembly_parser():
    assembly_parser = argparse.ArgumentParser(add_help=False,
        formatter_class=argparse.RawTextHelpFormatter)

    assembly_group = assembly_parser.add_argument_group('Assembly arguments')

    assembly_group.add_argument("--assembler", metavar="<STR>", default='mitoassemble',
        choices=['mitoassemble', 'spade', 'megahit'],
        help="Assembler to be used. [%(default)s]")

    assembly_group.add_argument("--insert_size", metavar="<INT>", default='250',
        help="insert size of input fastq files [%(default)s]")

    assembly_group.add_argument("--assembler_thread_number", metavar="<INT>",
        default='8', help=argparse.SUPPRESS)

    assembly_group.add_argument("--kmer", type=int,
        metavar="<INT>", default=71, nargs="+",
        help='kmer size(s) to be used [%(default)s]')

    return assembly_parser


def prepare_search_mito_parser():
    search_mito_parser = argparse.ArgumentParser(add_help=False,
        formatter_class=argparse.RawTextHelpFormatter)

    search_mito_group = search_mito_parser.add_argument_group('Search' +\
        ' mitochondrial sequences arguments')

    search_mito_group.add_argument("--filter_taxa_method", choices=[1,3],
        default=1, type=int,
        help='''1: filter out non-requiring_taxa sequences by mito-PCGs annotation
    to do taxa assignment. 3: do not filter [%(default)s]''')

    search_mito_group.add_argument('--min_abundance', metavar='<float>',
        type=float, default=10,
        help='the minimum abundance of sequence required [%(default)s]')

    search_mito_group.add_argument("--requiring_taxa", metavar="<STR>",
        default="Arthropoda",
        help='''filtering out non-requiring taxa sequences which may be
    contamination [%(default)s]''')

    search_mito_group.add_argument("--requiring_relax", default="0",
        choices=["0", "1", "2", "3", "4", "5", "6"],
        help='''The relaxing threshold for filtering non-target-requiring_taxa.
    The larger digital means more relaxing. [%(default)s]''')

    return search_mito_parser


def prepare_search_and_annot_mito_parser():
    search_and_annot_mito_parser = argparse.ArgumentParser(add_help=False,
        formatter_class=argparse.RawTextHelpFormatter)

    search_and_annot_mito_parser.add_argument("--genetic_code", metavar="<INT>",
        default="auto",
        help='''which genetic code table to use? 'auto' means determined by
    '--clade' option. [%(default)s]''')

    search_and_annot_mito_parser.add_argument("--clade", default="Arthropoda",
        choices=["Chordata", "Arthropoda"],
        # "Echinodermata", "Annelida-segmented-worms", "Bryozoa", "Mollusca", "Nematoda", "Nemertea-ribbon-worms", "Porifera-sponges"],
        help="which clade does your species belong to? [%(default)s]")

    return search_and_annot_mito_parser


def prepare_annotation_parser():
    annotation_parser = argparse.ArgumentParser(add_help=False,
        formatter_class=argparse.RawTextHelpFormatter)

    annotation_group = annotation_parser.add_argument_group('Annotation arguments')

    annotation_group.add_argument("--annotation", default=True,
        action="store_false", help="do annotation or not? [%(default)s]")

    annotation_group.add_argument("--species_name", metavar="<STR>",
        default="Test sp.",
        help='''species name to use in output genbank file ['%(default)s']''')

    return annotation_parser


###############################################################################

def get_parser_all(subparsers):
    ## all subcommand
    parser_all = subparsers.add_parser("all",
        parents=[
            prepare_common_parser(),
            prepare_fastq_parser(),
            prepare_filter_parser(),
            prepare_assembly_parser(),
            prepare_search_mito_parser(),
            prepare_search_and_annot_mito_parser(),
            prepare_annotation_parser()],
        help="run filter, assemble and annotate")

    parser_all.add_argument("--topology", choices=["linear", "circular"],
        default="linear", help=argparse.SUPPRESS)
                        #help="the sequences are circular?" +\
                        #" (seq length >= 12Kbp) (default: %(default)s)")

    parser_all.add_argument("--from_soaptrans", default=True,
                        action="store_true", help=argparse.SUPPRESS)

    return parser_all


def get_parser_all2(subparsers):
    parser_all2 = subparsers.add_parser("all2",
        parents=[
            prepare_common_parser(),
            prepare_fastq_parser(),
            prepare_assembly_parser(),
            prepare_search_mito_parser(),
            prepare_search_and_annot_mito_parser(),
            prepare_annotation_parser()],
        help="run assemble and annotate")

    parser_all2.add_argument("--topology", choices=["linear", "circular"],
        default="linear", help=argparse.SUPPRESS)

    parser_all2.add_argument("--usepre", default=False, action="store_true",
        help='''use previous results (clean data / first assembly by SOAPTrans)
    Otherwise I will delete all previous result files and run from the
    very beginning (default: %(default)s)''')

    parser_all2.add_argument("--from_soaptrans", default=True,
        action="store_true", help=argparse.SUPPRESS)

    return parser_all2


def get_parser_filter(subparsers):
    parser_filter = subparsers.add_parser("filter",
        parents=[
            prepare_common_parser(),
            prepare_fastq_parser(),
            prepare_filter_parser()],
        help="filter raw reads")

    return parser_filter


def get_parser_assemble(subparsers):
    parser_assemble = subparsers.add_parser("assemble",
        parents=[
            prepare_common_parser(),
            prepare_fastq_parser(),
            prepare_assembly_parser(),
            prepare_search_mito_parser(),
            prepare_search_and_annot_mito_parser()],
        help="do assembly from input fastq reads," +\
        " output mitosequences.")

    parser_assemble.add_argument("--usepre", default=False, action="store_true",
        help='''use previous results (clean data / first assembly by SOAPTrans)
    Otherwise I will delete all previous result files and run from the very
    beginning (default: %(default)s)''')

    return parser_assemble


def get_parser_findmitoscaf(subparsers):
    parser_findmitoscaf = subparsers.add_parser("findmitoscaf",
        parents=[
            prepare_common_parser(),
            prepare_fastafile_parser(),
            prepare_fastq_parser(),
            prepare_search_mito_parser(),
            prepare_search_and_annot_mito_parser()],
        help='''Search for mitochondrial sequences from assembly.
    About 2-3 Gbp fastq data is needed to calculate the average
    sequencing depth of each sequences, otherwise,
    '--from_soaptrans' should be used.''')

    parser_findmitoscaf.add_argument("--from_soaptrans",
        default=False,
        action="store_true",
        help='''is the input fasta generated by SOAPTrans?
    if not, '--fastq1' and '--fastq2' must be set! (default: %(default)s)''')

    return parser_findmitoscaf


def get_parser_annotate(subparsers):
    parser_annotate = subparsers.add_parser("annotate",
        parents=[prepare_common_parser(),
            prepare_fastafile_parser(),
            prepare_annotation_parser(),
            prepare_search_and_annot_mito_parser()],
        help="annotate PCGs, tRNA and rRNA genes.")

    parser_annotate.add_argument("--fastq1", metavar="<STR>", help="fastq 1 file, if you want to draw depth.")

    parser_annotate.add_argument("--fastq2", metavar="<STR>", help="fastq 2 file, if you want to draw depth.")

    parser_annotate.add_argument("--depth_file", metavar="<STR>",
        help=argparse.SUPPRESS)
                        #help="file of sequencing depth along the sequences")

    parser_annotate.add_argument("--topology", choices=["linear", "circular"],
        default="linear",
        help='''the sequences are circular? (seq length must be >= 12Kbp)
    (default: %(default)s)''')

    return parser_annotate


def get_parser_visualize(subparsers):
    parser_visualize = subparsers.add_parser("visualize",
        help="visualization of GenBank file")

    return subparsers


################################################################

def get_para(options):
    _version = '3'
    # textwrap.dedent
    description = """
    Description

        MitoZ - A toolkit for mitochondrial genome assembly,
            annotation and visualization

    Version
        {version}

    Citation

    Guanliang Meng, Yiyuan Li, Chentao Yang, Shanlin Liu.
    MitoZ: A toolkit for mitochondrial genome assembly, annotation
    and visualization; doi: https://doi.org/10.1093/nar/gkz173

    THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
    REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
    AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
    INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
    LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
    OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
    PERFORMANCE OF THIS SOFTWARE. REVERSE ENGINEERING IS STRICTLY PROHIBITED.

    COPYRIGHT © 2019-2022 Guanliang Meng. ALL RIGHTS RESERVED.
    """.format(version=_version)


    parser = argparse.ArgumentParser(prog="MitoZ", description=description,
            formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("--version", action="version", version="%(prog)s " + _version)

    subparsers = parser.add_subparsers(dest='command')

    get_parser_all(subparsers)
    get_parser_all2(subparsers)
    get_parser_filter(subparsers)
    get_parser_assemble(subparsers)
    get_parser_findmitoscaf(subparsers)
    get_parser_annotate(subparsers)
    get_parser_visualize(subparsers)

    if len(options) == 0:
        parser.print_help()
        sys.exit()

    args = parser.parse_args(options)


if __name__ == '__main__':
    get_para(sys.argv[1:])