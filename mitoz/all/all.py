#!/usr/bin/env python3

"""
This file is part of MitoZ.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
PERFORMANCE OF THIS SOFTWARE. 

WHEN YOU ADAPT (PART OF) THE SOFTWARE FOR YOUR USE CASES, THE AUTHOR AND
THE SOFTWARE MUST BE EXPLICITLY CREDITED IN YOUR PUBLICATIONS AND SOFTWARE,
AND YOU SHOULD ASK THE USERS OF YOUR SOFTWARE TO CITE THE SOFTWARE IN
THEIR PUBLICATIONS. IN A WORD, 请讲武德.

COPYRIGHT 2019-2022 Guanliang Meng. ALL RIGHTS RESERVED.

"""

import argparse
import sys
import os
import re
import subprocess
import time
from glob import glob
from mitoz.utility import utility
from mitoz import findmitoscaf
from mitoz.utility.utility import gather_result, runcmd, files_exist_0_or_1, file_not_empty
from mitoz.utility.utility import pre_del_cmd, abspath
import copy
from pathlib import Path

import mitoz

python3 = sys.executable

mitoz_pkg_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
findmitoscaf_script_dir = os.path.join(mitoz_pkg_dir, 'findmitoscaf', 'script')
assemble_script_dir = os.path.join(mitoz_pkg_dir, 'assemble', 'script')
profiles_dir = os.path.join(mitoz_pkg_dir, 'profiles')
annotate_script_dir = os.path.join(mitoz_pkg_dir, 'annotate', 'script')
template_sbt = os.path.join(annotate_script_dir, 'template.sbt')

def add_arguments(parser):
    common_group = parser.add_argument_group('Common arguments')

    common_group.add_argument('--outprefix', metavar='<str>', default='out',
        help="output prefix [%(default)s]")

    common_group.add_argument('--thread_number', metavar='<int>', default='8',
        help="thread number [%(default)s]")

    common_group.add_argument('--workdir', metavar='<directory>', default='./',
        help="working directory [%(default)s]")

    '''
    common_group.add_argument('--workdir_done', metavar='<directory>', default='./done',
        help="done directory [%(default)s]")

    common_group.add_argument('--workdir_log', metavar='<directory>', default='./log',
        help="log directory [%(default)s]")
    '''

    common_group.add_argument("--clade", default="Arthropoda",
        choices=["Chordata", "Arthropoda",
            "Echinodermata", "Annelida-segmented-worms",
            "Bryozoa", "Mollusca", "Nematoda",
            "Nemertea-ribbon-worms",
            "Porifera-sponges"],
        help="which clade does your species belong to? [%(default)s]")

    common_group.add_argument("--genetic_code", metavar="<INT>",
        default="auto",
        help='''which genetic code table to use? 'auto' means determined by
            '--clade' option. [%(default)s]''')

    common_group.add_argument("--species_name", metavar="<STR>",
        default="Test sp.",
        help='''species name to use in output genbank file ['%(default)s']''')

    common_group.add_argument("--template_sbt", metavar="<file>",
        default=template_sbt,
        help='''The sqn template to generate the resulting genbank file. 
        Go to https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/#Template
        to generate your own template file if you like. ['%(default)s']''')


    fastq_group = parser.add_argument_group('Input fastq information')
    fastq_group.add_argument('--fq1', metavar='<file>', required=True,
        help="Fastq1 file [required]")

    fastq_group.add_argument('--fq2', metavar='<file>',
        help="Fastq2 file [optional]")

    fastq_group.add_argument('--phred64', action='store_true',
        help="Are the fastq phred64 encoded? [%(default)s]")

    fastq_group.add_argument("--insert_size", metavar="<INT>", default='250',
        help="insert size of input fastq files [%(default)s]")

    fastq_group.add_argument("--fastq_read_length", metavar="<INT>", default='150',
        help='''read length of fastq reads, used by the filter subcommand and mitoAssemble. [%(default)s]''')

    # fastq_group.add_argument("--data_size_for_mt_assembly", metavar="<float>", type=float, default=5,
     #   help='''Data size (GB) used for mitochondrial genome assembly, it is recommmended to be less than 8 GB. Set to 0 if you want to use ALL clean data. [%(default)s]''')
    fastq_group.add_argument("--data_size_for_mt_assembly", metavar="<float1>,<float2>", default="2,0",
        help='''Data size (Gbp) used for mitochondrial genome assembly, usually between 2~8 Gbp is enough.
        The float1 means the size (Gbp) of raw data to be subsampled, while the float2 means the size of 
        clean data must be >= float2 Gbp, otherwise MitoZ will STOP running! When only float1 is set,
        float2 is assumed to be 0.  (1) Set float1 to be 0 if you want to use ALL raw data;
        (2) Set 0,0 if you want to use ALL raw data and do NOT interrupt MitoZ even if you got very little clean data.
        If you got missing mitochondrial genes, try (1) differnt kmers;
        (2)different assembler; (3) increase <float1>,<float2> [%(default)s]''')


    fastq_group.add_argument('--skip_filter', action='store_true',
        help='''Skip the rawdata filtering step, assuming input fastq are clean data.
        To subsample such clean data, set <float2> of the --data_size_for_mt_assembly
        option to be larger than 0 (using all input clean data by default). [%(default)s]''')

    fastq_group.add_argument('--filter_other_para', metavar='<str>', default='',
        help="other parameter for filtering. [%(default)s]")


    #### from assemble

    assembly_group = parser.add_argument_group('Assembly arguments')

    assembly_group.add_argument("--assembler", default='megahit',
        choices=['mitoassemble', 'spades', 'megahit'],
        help="Assembler to be used. [%(default)s]")

    assembly_group.add_argument("--tmp_dir", metavar="<STR>", 
        help="Set temp directory for megahit if necessary (See https://github.com/linzhi2013/MitoZ/issues/176)")

    assembly_group.add_argument("--kmers", type=int,
        metavar="<INT>", default=[71], nargs="+",
        help='kmer size(s) to be used. Multiple kmers can be used, separated by space %(default)s')

    assembly_group.add_argument("--kmers_megahit",
        metavar="<INT>", default=['43', '71', '99'], nargs="+",
        help="kmer size(s) to be used. Multiple kmers can be used, separated by space. Only for megahit [43 71 99]")

    assembly_group.add_argument("--kmers_spades",
        metavar="<INT>", default=['auto'], nargs="+",
        help='kmer size(s) to be used. Multiple kmers can be used, separated by space. Only for spades %(default)s')

    assembly_group.add_argument("--memory", metavar="<INT>", default='50',
        help="memory size limit for spades/megahit, no enough memory will make the two programs halt or exit [%(default)s]")

    assembly_group.add_argument("--resume_assembly", action='store_true',
        help="to resume previous assembly running [%(default)s]")

    search_mito_group = parser.add_argument_group('Search mitochondrial sequences arguments')

    search_mito_group.add_argument("--profiles_dir", metavar="<STR>",
        default=profiles_dir,
        help="Directory cotaining 'CDS_HMM/', 'MT_database/' and 'rRNA_CM/'. [%(default)s]")

    search_mito_group.add_argument("--slow_search", action='store_true',
        help='''By default, we firstly use tiara to perform quick sequence classification (100 times faster than usual!),
however, it is valid only when your mitochondrial sequences are >= 3000 bp.
If you have missing genes, set '--slow_search' to use the tradicitiona search mode. [%(default)s]''')

    search_mito_group.add_argument("--filter_by_taxa", action='store_false',
        help='''filter out non-requiring_taxa sequences by mito-PCGs annotation
    to do taxa assignment.[%(default)s]''')

    search_mito_group.add_argument("--requiring_taxa", metavar="<STR>",
        required=True,
        help='''filtering out non-requiring taxa sequences which may be
    contamination [required]''')


    search_mito_group.add_argument("--requiring_relax", default="0",
        choices=["0", "1", "2", "3", "4", "5", "6"],
        help='''The relaxing threshold for filtering non-target-requiring_taxa.
    The larger digital means more relaxing. [%(default)s]''')

    # search_mito_group.add_argument("--skip_read_mapping", action='store_true',
    #    help='''Skip read-mapping step, assuming we can extract the abundance from seqid line. [%(default)s]''')

    search_mito_group.add_argument('--min_abundance', metavar='<float>',
        default='10',
        help='''the minimum abundance of sequence required. 
        Set this to any value <= 0 if you do NOT want to filter sequences by abundance [%(default)s]''')

    return parser


class Filter_Args():
    def __init__(self, raw_args):
        return self.prepare_new_args(raw_args)

    def prepare_new_args(self, raw_args):
        self.logger = raw_args.logger
        self.fq1 = ''
        self.fq2 = ''
        if raw_args.fq1 and Path(raw_args.fq1).is_file():
            self.fq1 = os.path.abspath(raw_args.fq1)
        if raw_args.fq2 and Path(raw_args.fq2).is_file():
            self.fq2 = os.path.abspath(raw_args.fq2)

        self.phred64 = raw_args.phred64
        self.outprefix = raw_args.outprefix
        self.fastq_read_length = int(raw_args.fastq_read_length)
        self.data_size_for_mt_assembly = raw_args.data_size_for_mt_assembly
        self.filter_other_para = raw_args.filter_other_para
        self.thread_number = raw_args.thread_number
        self.workdir = os.path.join(os.path.abspath(raw_args.workdir), 'clean_data')
        self.workdir_done = self.workdir
        self.workdir_log = self.workdir

    def __str__(self):
        return str(self.__dict__)


class Assemble_Args():
    def __init__(self, raw_args):
        return self.prepare_new_args(raw_args)

    def prepare_new_args(self, raw_args):
        self.logger = raw_args.logger
        self.workdir = os.path.join(os.path.abspath(raw_args.workdir), 'mt_assembly')
        self.outprefix = raw_args.outprefix
        self.thread_number = raw_args.thread_number
        
        self.fq1 = ''
        self.fq2 = ''
        if raw_args.fq1 and Path(raw_args.fq1).is_file():
            self.fq1 = os.path.abspath(raw_args.fq1)
        if raw_args.fq2 and Path(raw_args.fq2).is_file():
            self.fq2 = os.path.abspath(raw_args.fq2)

        self.insert_size = raw_args.insert_size
        self.fastq_read_length = raw_args.fastq_read_length
        self.assembler = raw_args.assembler
        self.tmp_dir =raw_args.tmp_dir
        self.kmers = raw_args.kmers
        self.kmers_megahit = raw_args.kmers_megahit
        self.kmers_spades = raw_args.kmers_spades
        self.memory = raw_args.memory
        self.resume_assembly = raw_args.resume_assembly
        self.profiles_dir = raw_args.profiles_dir
        self.slow_search = raw_args.slow_search
        self.filter_by_taxa = raw_args.filter_by_taxa
        self.requiring_taxa = raw_args.requiring_taxa
        self.requiring_relax = raw_args.requiring_relax
        self.abundance_pattern = r'abun\=([0-9]+\.*[0-9]*)'
        self.min_abundance = raw_args.min_abundance
        # self.skip_read_mapping = raw_args.skip_read_mapping
        # self.skip_read_mapping = True
        self.genetic_code = raw_args.genetic_code
        self.clade = raw_args.clade

    def __str__(self):
        return str(self.__dict__)


class Annotate_Args():
    def __init__(self, raw_args):
        return self.prepare_new_args(raw_args)

    def prepare_new_args(self, raw_args):
        self.logger = raw_args.logger
        self.workdir = os.path.join(os.path.abspath(raw_args.workdir), 'mt_annotation')
        self.outprefix = raw_args.outprefix
        self.thread_number = raw_args.thread_number
        
        self.fq1 = ''
        self.fq2 = ''
        if raw_args.fq1 and Path(raw_args.fq1).is_file():
            self.fq1 = os.path.abspath(raw_args.fq1)
        if raw_args.fq2 and Path(raw_args.fq2).is_file():
            self.fq2 = os.path.abspath(raw_args.fq2)

        self.fastafiles = raw_args.fastafiles
        self.profiles_dir = raw_args.profiles_dir
        self.species_name = raw_args.species_name
        self.template_sbt = raw_args.template_sbt
        self.genetic_code = raw_args.genetic_code
        self.clade = raw_args.clade
        self.mitoz_module_used = 'assemble'

    def __str__(self):
        return str(self.__dict__)


def main(args):
    logger = args.logger
    logger.info("all.main() got args:\n{}".format(args))

    args.workdir = os.path.abspath(args.workdir)

    if args.fq1:
        if Path(args.fq1).is_file():
            args.fq1 = os.path.abspath(args.fq1)
        else:
            logger.error("Can NOT access {} !!".format(args.fq1))
            sys.exit(0)
    else:
        args.fq1 = ''

    if args.fq2:
        if Path(args.fq2).is_file():
            args.fq2 = os.path.abspath(args.fq2)
        else:
            logger.error("Can NOT access {} !!".format(args.fq2))
            sys.exit(0)
    else:
        args.fq2 = ''
    
    if args.tmp_dir:
        args.tmp_dir = os.path.abspath(args.tmp_dir)

    if args.template_sbt:
        if Path(args.template_sbt).is_file():
            args.template_sbt = os.path.abspath(args.template_sbt)
        else:
            logger.error("Can NOT access {} !!".format(args.template_sbt))
            sys.exit(0)

    logger.info("all.main() got updated args:\n{}".format(args))

    if args.assembler == 'spades' and not (args.fq1 and args.fq2):
        logger.error("You have to provide SPAdes with both fq1 and fq2!")
        sys.exit()

    final_result_dir = os.path.join(args.workdir, args.outprefix + '.result')

    if not args.skip_filter:
        filter_args = Filter_Args(args)
        filtered_fq1, filtered_fq2 = mitoz.filter.main(filter_args)
        args.fq1 = filtered_fq1
        args.fq2 = filtered_fq2

    # check if we need to subsample
    line = [float(j) for j in args.data_size_for_mt_assembly.split(',')]
    required_clean_data_size = 0
    if len(line) >= 2:
        required_clean_data_size = float(line[1])
    
    if args.skip_filter and required_clean_data_size:
        sub_fq1, sub_fq2 = mitoz.filter.main_subsampling_fq(
            workdir=args.workdir,
            fastq_read_length=args.fastq_read_length,
            required_data_size=required_clean_data_size,
            fq1=args.fq1,
            fq2=args.fq2,
            logger=args.logger)

        args.fq1 = sub_fq1
        args.fq2 = sub_fq2


    assemble_args = Assemble_Args(args)
    args.fastafiles, assemble_all_result_wdir = mitoz.assemble.main(assemble_args)

    annotate_args = Annotate_Args(args)
    resulting_gb_files, annotate_all_result_wdir = mitoz.annotate.main(annotate_args)

    gather_result(*assemble_all_result_wdir, *annotate_all_result_wdir, logger=args.logger, result_wdir=final_result_dir)

