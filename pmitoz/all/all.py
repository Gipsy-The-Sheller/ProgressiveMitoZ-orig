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
from pmitoz.utility import utility
from pmitoz import findmitoscaf
from pmitoz.utility.utility import gather_result, runcmd, files_exist_0_or_1, file_not_empty
from pmitoz.utility.utility import pre_del_cmd, abspath
import copy
from pathlib import Path

import pmitoz as mitoz  # Keep mitoz alias for compatibility

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

    # Iterative assembly arguments
    iterative_group = parser.add_argument_group('Iterative assembly arguments')
    
    iterative_group.add_argument('--iter', metavar='<int>', type=int, default=1,
        help='''Number of iterative assembly rounds after initial findmitoscaf.
        Each iteration will collect reads mapping to current mitogenome candidates
        and reassemble them to improve assembly quality. [%(default)s]''')

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
        self.iter = getattr(raw_args, 'iter', 1)  # Default to 1 if not specified

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

    # Perform iterative assembly if iter > 1
    if args.iter > 1:
        logger.info(f"Starting iterative assembly with {args.iter} iterations...")
        args.fastafiles = perform_iterative_assembly(args, assemble_all_result_wdir, logger)

    annotate_args = Annotate_Args(args)
    resulting_gb_files, annotate_all_result_wdir = mitoz.annotate.main(annotate_args)

    gather_result(*assemble_all_result_wdir, *annotate_all_result_wdir, logger=args.logger, result_wdir=final_result_dir)


def perform_iterative_assembly(args, assemble_all_result_wdir, logger):
    """
    Perform iterative assembly after initial findmitoscaf.
    
    Args:
        args: Command line arguments
        assemble_all_result_wdir: Result directory from initial assembly
        logger: Logger instance
    
    Returns:
        List of final fasta files after iterative assembly
    """
    from pmitoz.utility.utility import runcmd, abspath
    from Bio import SeqIO
    import shutil
    
    # In iterative mode, we need to run findmitoscaf on the initial assembly
    # Find the assembly file from the initial assembly
    # The assembly file should be in the mt_assembly directory
    initial_assembly_file = None
    
    # First try to find from workdir/mt_assembly
    mt_assembly_dir = os.path.join(args.workdir, 'mt_assembly')
    if os.path.exists(mt_assembly_dir):
        for root, dirs, files in os.walk(mt_assembly_dir):
            for file in files:
                if file.endswith('.reformatted.fa'):
                    initial_assembly_file = os.path.join(root, file)
                    logger.info(f"Found initial assembly file: {initial_assembly_file}")
                    break
            if initial_assembly_file:
                break
    
    # If not found, try to find from result directories
    if not initial_assembly_file:
        for result_dir in assemble_all_result_wdir:
            if os.path.exists(result_dir):
                for root, dirs, files in os.walk(result_dir):
                    for file in files:
                        if file.endswith('.reformatted.fa'):
                            initial_assembly_file = os.path.join(root, file)
                            logger.info(f"Found initial assembly file from result dir: {initial_assembly_file}")
                            break
                    if initial_assembly_file:
                        break
            if initial_assembly_file:
                break
    
    # If still not found, check if args.fastafiles has any valid files
    if not initial_assembly_file and args.fastafiles:
        for fastafile in args.fastafiles:
            if fastafile and os.path.exists(fastafile):
                initial_assembly_file = fastafile
                logger.info(f"Using assembly file from args.fastafiles: {initial_assembly_file}")
                break
    
    if not initial_assembly_file:
        logger.error("No initial assembly file found for iterative assembly")
        logger.error(f"Searched in: {mt_assembly_dir} and result directories")
        return args.fastafiles if args.fastafiles else []
    
    # Run findmitoscaf on the initial assembly
    logger.info(f"Running findmitoscaf on initial assembly: {initial_assembly_file}")
    initial_mt_file = run_initial_findmitoscaf(initial_assembly_file, args, logger)
    
    if not initial_mt_file or not os.path.exists(initial_mt_file):
        logger.error(f"Initial findmitoscaf failed or no mitogenome candidates found")
        return args.fastafiles
    
    logger.info(f"Starting iterative assembly with {args.iter} iterations")
    logger.info(f"Initial mitogenome candidates: {initial_mt_file}")
    
    current_mt_file = initial_mt_file
    
    for iteration in range(1, args.iter):
        logger.info(f"=== Iteration {iteration}/{args.iter-1} ===")
        
        # Create iteration-specific directories
        iter_workdir = os.path.join(args.workdir, f'iterative_assembly_{iteration}')
        iter_reads_dir = os.path.join(iter_workdir, 'collected_reads')
        iter_assembly_dir = os.path.join(iter_workdir, 'assembly')
        
        os.makedirs(iter_reads_dir, exist_ok=True)
        # os.makedirs(iter_assembly_dir, exist_ok=True)
        # MEGAHIT will create the directory itself. If it exists, then MEGAHIT will encounter an error.
        # remove the directory if it exists
        if os.path.exists(iter_assembly_dir):
            logger.info(f"Removing existing reassembly directory: {iter_assembly_dir}")
            shutil.rmtree(iter_assembly_dir)
        
        # Step 1: Collect reads mapping to current mitogenome candidates
        logger.info(f"Collecting reads mapping to mitogenome candidates...")
        collected_reads = collect_mapping_reads(
            current_mt_file, 
            args.fq1, 
            args.fq2, 
            iter_reads_dir, 
            logger,
            args.thread_number
        )
        
        if not collected_reads:
            logger.warning(f"No reads collected in iteration {iteration}, stopping iterative assembly")
            break
            
        # Step 2: Reassemble with collected reads
        logger.info(f"Reassembling with collected reads...")
        new_mt_file = reassemble_with_reads(
            collected_reads,
            iter_assembly_dir,
            args,
            logger
        )
        
        if not new_mt_file or not os.path.exists(new_mt_file):
            logger.warning(f"Reassembly failed in iteration {iteration}, stopping iterative assembly")
            break
            
        # Step 3: Find new mitogenome candidates
        logger.info(f"Finding new mitogenome candidates...")
        current_mt_file = find_new_mitogenome_candidates(
            new_mt_file,
            iter_assembly_dir,
            args,
            logger
        )
        
        if not current_mt_file:
            logger.warning(f"No new mitogenome candidates found in iteration {iteration}, stopping iterative assembly")
            break
            
        logger.info(f"Iteration {iteration} completed. New candidates: {current_mt_file}")
    
    # Return the final mitogenome candidates
    final_fastafiles = [current_mt_file] if current_mt_file else args.fastafiles
    logger.info(f"Iterative assembly completed. Final candidates: {final_fastafiles}")
    
    return final_fastafiles


def run_initial_findmitoscaf(assembly_file, args, logger):
    """
    Run findmitoscaf on the initial assembly to get mitogenome candidates.
    
    Args:
        assembly_file: Path to the initial assembly file
        args: Command line arguments
        logger: Logger instance
    
    Returns:
        Path to mitogenome candidates file
    """
    from pmitoz import findmitoscaf
    
    # Create findmitoscaf arguments
    findmitoscaf_args = type('Args', (), {})()
    findmitoscaf_args.fastafile = assembly_file
    findmitoscaf_args.fq1 = args.fq1
    findmitoscaf_args.fq2 = args.fq2
    findmitoscaf_args.outprefix = args.outprefix + "_initial"
    findmitoscaf_args.workdir = os.path.join(args.workdir, 'initial_findmitoscaf')
    findmitoscaf_args.thread_number = args.thread_number
    findmitoscaf_args.profiles_dir = args.profiles_dir
    findmitoscaf_args.slow_search = args.slow_search
    findmitoscaf_args.filter_by_taxa = args.filter_by_taxa
    findmitoscaf_args.requiring_taxa = args.requiring_taxa
    findmitoscaf_args.requiring_relax = args.requiring_relax
    findmitoscaf_args.min_abundance = args.min_abundance
    findmitoscaf_args.abundance_pattern = r'abun\=([0-9]+\.*[0-9]*)'
    findmitoscaf_args.skip_read_mapping = False  # Need to map reads to calculate abundance
    findmitoscaf_args.genetic_code = getattr(args, 'genetic_code', 'auto')
    findmitoscaf_args.clade = getattr(args, 'clade', 'Arthropoda')
    findmitoscaf_args.logger = logger
    
    logger.info(f"Running initial findmitoscaf on {assembly_file}")
    
    try:
        mt_file = findmitoscaf.main(findmitoscaf_args)
        # Verify the returned file exists
        if mt_file and os.path.exists(mt_file):
            logger.info(f"Initial findmitoscaf completed successfully: {mt_file}")
            return os.path.abspath(mt_file)
        else:
            logger.warning(f"findmitoscaf returned path that doesn't exist: {mt_file}")
            # Try to find the file in the workdir
            workdir = findmitoscaf_args.workdir
            potential_files = [
                os.path.join(workdir, findmitoscaf_args.outprefix + '.mitogenome.fa'),
                os.path.join(workdir, '*.mitogenome.fa'),
            ]
            for pattern in potential_files:
                if '*' in pattern:
                    import glob
                    matches = glob.glob(pattern)
                    if matches:
                        logger.info(f"Found mitogenome file: {matches[0]}")
                        return os.path.abspath(matches[0])
                elif os.path.exists(pattern):
                    logger.info(f"Found mitogenome file: {pattern}")
                    return os.path.abspath(pattern)
            logger.error("Could not find mitogenome file after findmitoscaf")
            return None
    except Exception as e:
        logger.error(f"Initial findmitoscaf failed: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return None


def collect_mapping_reads(mt_file, fq1, fq2, output_dir, logger, thread_number):
    """
    Collect reads mapping to mitogenome candidates using BWA.
    
    Args:
        mt_file: Path to mitogenome candidates file
        fq1: Path to forward reads
        fq2: Path to reverse reads  
        output_dir: Output directory for collected reads
        logger: Logger instance
    
    Returns:
        List of collected read files
    """
    from pmitoz.utility.utility import runcmd
    
    # Build BWA index
    logger.info(f"Building BWA index for {mt_file}")
    index_cmd = f"bwa index {mt_file}"
    runcmd(index_cmd, logger=logger)
    
    # Map reads
    bam_file = os.path.join(output_dir, "mapping.bam")
    sorted_bam = os.path.join(output_dir, "mapping_sorted.bam")
    
    if fq1 and fq2:
        map_cmd = f"bwa mem -t {thread_number} {mt_file} {fq1} {fq2} | samtools view -b -h -F 4 > {bam_file}"
    elif fq1:
        map_cmd = f"bwa mem -t {thread_number} {mt_file} {fq1} | samtools view -b -h -F 4 > {bam_file}"
    else:
        logger.error("No input reads provided")
        return []
    
    logger.info(f"Mapping reads: {map_cmd}")
    runcmd(map_cmd, logger=logger)
    
    # # Convert to BAM and sort
    # logger.info("Converting SAM to BAM and sorting...")
    # runcmd(f"samtools view -bS {sam_file} > {bam_file}", logger=logger)
    # runcmd(f"samtools sort {bam_file} -o {sorted_bam}", logger=logger)
    
    # Extract mapped reads
    collected_fq1 = os.path.join(output_dir, "collected_R1.fq")
    collected_fq2 = os.path.join(output_dir, "collected_R2.fq")
    
    if fq1 and fq2:
        # Extract paired reads
        runcmd(f"samtools fastq -1 {collected_fq1} -2 {collected_fq2} {sorted_bam}", logger=logger)
        return [collected_fq1, collected_fq2]
    else:
        # Extract single reads
        runcmd(f"samtools fastq {sorted_bam} > {collected_fq1}", logger=logger)
        return [collected_fq1]


def reassemble_with_reads(collected_reads, output_dir, args, logger):
    """
    Reassemble with collected reads using the same assembler as initial assembly.
    
    Args:
        collected_reads: List of collected read files
        output_dir: Output directory for reassembly
        args: Command line arguments
        logger: Logger instance
    
    Returns:
        Path to reassembled contigs file
    """
    from pmitoz.utility.utility import runcmd
    
    # Prepare reads
    fq1 = collected_reads[0] if len(collected_reads) > 0 else None
    fq2 = collected_reads[1] if len(collected_reads) > 1 else None
    
    # Use the same assembler as initial assembly
    if args.assembler == 'megahit':
        return reassemble_with_megahit(fq1, fq2, output_dir, args, logger)
    elif args.assembler == 'spades':
        return reassemble_with_spades(fq1, fq2, output_dir, args, logger)
    elif args.assembler == 'mitoassemble':
        return reassemble_with_mitoassemble(fq1, fq2, output_dir, args, logger)
    else:
        logger.error(f"Unsupported assembler for iterative assembly: {args.assembler}")
        return None


def reassemble_with_megahit(fq1, fq2, output_dir, args, logger):
    """Reassemble with MEGAHIT"""
    from pmitoz.utility.utility import runcmd
    import sys
    
    megahit_script_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'assemble', 'script', 'megahit')
    python3 = sys.executable
    
    cmd = ["megahit", "--out-dir", output_dir, "--num-cpu-threads", str(args.thread_number)]
    
    if fq1 and fq2:
        cmd.extend(["-1", fq1, "-2", fq2])
    elif fq1:
        cmd.extend(["-r", fq1])
    
    # Use k-mer list for iterative assembly (typically smaller k-mers for better assembly)
    if hasattr(args, 'kmers_megahit') and args.kmers_megahit:
        k_list = ",".join(map(str, args.kmers_megahit))
        cmd.extend(["--k-list", k_list])
    
    cmd_str = " ".join(cmd)
    logger.info(f"Reassembling with MEGAHIT: {cmd_str}")
    runcmd(cmd_str, logger=logger)
    
    final_contigs = os.path.join(output_dir, "final.contigs.fa")
    if not os.path.exists(final_contigs):
        return None
    
    # Reformat the contigs file for findmitoscaf
    assembly_file_reformated = final_contigs + '.reformatted.fa'
    soft = os.path.join(megahit_script_dir, 'reformat_megahit_scafSeq.py')
    if os.path.exists(soft):
        command = python3 + " " + soft +\
            ' {0} '.format(final_contigs) +\
            ' {0} '.format(assembly_file_reformated)
        runcmd(command, logger=logger)
        if os.path.exists(assembly_file_reformated):
            return assembly_file_reformated
    
    # If reformat script not found, return original file
    return final_contigs


def reassemble_with_spades(fq1, fq2, output_dir, args, logger):
    """Reassemble with SPAdes"""
    from pmitoz.utility.utility import runcmd
    
    cmd = ["spades.py", "--threads", str(args.thread_number), "-o", output_dir]
    
    if fq1 and fq2:
        cmd.extend(["-1", fq1, "-2", fq2])
    elif fq1:
        cmd.extend(["-s", fq1])
    
    # Use k-mer list for iterative assembly
    if hasattr(args, 'kmers_spades') and args.kmers_spades and args.kmers_spades != ['auto']:
        k_list = ",".join(map(str, args.kmers_spades))
        cmd.extend(["-k", k_list])
    
    cmd_str = " ".join(cmd)
    logger.info(f"Reassembling with SPAdes: {cmd_str}")
    runcmd(cmd_str, logger=logger)
    
    final_contigs = os.path.join(output_dir, "contigs.fasta")
    if not os.path.exists(final_contigs):
        return None
    
    # Reformat the contigs file for findmitoscaf
    assembly_file_reformated = final_contigs + '.reformatted.fa'
    soft = os.path.join(spades_script_dir, 'reformat_spades_scafSeq.py')
    if os.path.exists(soft):
        command = python3 + " " + soft +\
            ' {0} '.format(final_contigs) +\
            ' {0} '.format(assembly_file_reformated)
        runcmd(command, logger=logger)
        if os.path.exists(assembly_file_reformated):
            return assembly_file_reformated
    
    # If reformat script not found, return original file
    return final_contigs


def reassemble_with_mitoassemble(fq1, fq2, output_dir, args, logger):
    """Reassemble with mitoAssemble"""
    from pmitoz.utility.utility import runcmd
    
    # mitoAssemble requires specific input format
    # This is a simplified implementation
    cmd = ["mitoAssemble", "--workdir", output_dir]
    
    if fq1 and fq2:
        cmd.extend(["--fq1", fq1, "--fq2", fq2])
    elif fq1:
        cmd.extend(["--fq1", fq1])
    
    cmd_str = " ".join(cmd)
    logger.info(f"Reassembling with mitoAssemble: {cmd_str}")
    runcmd(cmd_str, logger=logger)
    
    # mitoAssemble output location may vary
    final_contigs = os.path.join(output_dir, "final.contigs.fa")
    if not os.path.exists(final_contigs):
        final_contigs = os.path.join(output_dir, "contigs.fa")
    
    if not os.path.exists(final_contigs):
        return None
    
    # Reformat the contigs file for findmitoscaf
    assembly_file_reformated = final_contigs + '.reformatted.fa'
    soft = os.path.join(mitoassemble_script_dir, 'reformat_mitoassemble_scafSeq.py')
    if os.path.exists(soft):
        command = python3 + " " + soft +\
            ' {0} '.format(final_contigs) +\
            ' {0} '.format(assembly_file_reformated)
        runcmd(command, logger=logger)
        if os.path.exists(assembly_file_reformated):
            return assembly_file_reformated
    
    # If reformat script not found, return original file
    return final_contigs


def find_new_mitogenome_candidates(contigs_file, output_dir, args, logger):
    """
    Find new mitogenome candidates from reassembled contigs using findmitoscaf.
    
    Args:
        contigs_file: Path to reassembled contigs
        output_dir: Output directory
        args: Command line arguments
        logger: Logger instance
    
    Returns:
        Path to new mitogenome candidates file
    """
    from pmitoz import findmitoscaf
    from pmitoz.utility.utility import runcmd
    
    # Create findmitoscaf arguments
    findmitoscaf_args = type('Args', (), {})()
    findmitoscaf_args.fastafile = contigs_file
    findmitoscaf_args.fq1 = args.fq1
    findmitoscaf_args.fq2 = args.fq2
    findmitoscaf_args.outprefix = args.outprefix + f"_iter"
    findmitoscaf_args.workdir = output_dir
    findmitoscaf_args.thread_number = args.thread_number
    findmitoscaf_args.profiles_dir = args.profiles_dir
    findmitoscaf_args.slow_search = args.slow_search
    findmitoscaf_args.filter_by_taxa = args.filter_by_taxa
    findmitoscaf_args.requiring_taxa = args.requiring_taxa
    findmitoscaf_args.requiring_relax = args.requiring_relax
    findmitoscaf_args.min_abundance = args.min_abundance
    findmitoscaf_args.abundance_pattern = r'abun\=([0-9]+\.*[0-9]*)'
    findmitoscaf_args.skip_read_mapping = False  # Need to map reads to calculate abundance
    findmitoscaf_args.genetic_code = getattr(args, 'genetic_code', 'auto')
    findmitoscaf_args.clade = getattr(args, 'clade', 'Arthropoda')
    findmitoscaf_args.logger = logger
    
    logger.info(f"Finding mitogenome candidates from {contigs_file}")
    
    try:
        mt_file = findmitoscaf.main(findmitoscaf_args)
        # Verify the returned file exists
        if mt_file and os.path.exists(mt_file):
            logger.info(f"findmitoscaf iteration completed successfully: {mt_file}")
            return os.path.abspath(mt_file)
        else:
            logger.warning(f"findmitoscaf returned path that doesn't exist: {mt_file}")
            # Try to find the file in the workdir
            workdir = findmitoscaf_args.workdir
            potential_files = [
                os.path.join(workdir, findmitoscaf_args.outprefix + '.mitogenome.fa'),
            ]
            for pattern in potential_files:
                import glob
                if '*' in pattern:
                    matches = glob.glob(pattern)
                    if matches:
                        logger.info(f"Found mitogenome file: {matches[0]}")
                        return os.path.abspath(matches[0])
                elif os.path.exists(pattern):
                    logger.info(f"Found mitogenome file: {pattern}")
                    return os.path.abspath(pattern)
            logger.error("Could not find mitogenome file after findmitoscaf")
            return None
    except Exception as e:
        logger.error(f"findmitoscaf failed: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return None

