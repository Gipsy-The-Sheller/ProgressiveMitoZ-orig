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
from mitoz.utility import utility
from mitoz import findmitoscaf
from mitoz.utility.utility import gather_result, runcmd, files_exist_0_or_1, file_not_empty
from mitoz.utility.utility import pre_del_cmd, abspath, find_subdirs_with_suffix
import copy
from pathlib import Path

python3 = sys.executable

mitoz_pkg_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
findmitoscaf_script_dir = os.path.join(mitoz_pkg_dir, 'findmitoscaf', 'script')
assemble_script_dir = os.path.join(mitoz_pkg_dir, 'assemble', 'script')
profiles_dir = os.path.join(mitoz_pkg_dir, 'profiles')

def add_arguments(parser):
    common_group = parser.add_argument_group('Common arguments')

    common_group.add_argument("--workdir", metavar="<STR>", default='./',
        help="workdir [%(default)s]")

    common_group.add_argument("--outprefix", metavar="<STR>", required=True,
        help="output prefix")

    common_group.add_argument("--thread_number", metavar="<INT>", default="8",
        help="thread number. Caution: For spades, --thread_number 32 can take 150 GB RAM! Setting this to 8 to 16 is typically good. [%(default)s]")

    fastq_group = parser.add_argument_group('Input fastq files')

    fastq_group.add_argument("--fq1", metavar="<file>", required=True, help="fastq 1 file. Set only this option but not --fastq2 means SE data. [required]")

    fastq_group.add_argument("--fq2", metavar="<file>", default='', help="fastq 2 file (optional for mitoassemble and megahit, required for spades)")

    fastq_group.add_argument("--insert_size", metavar="<INT>", default='250',
        help="insert size of input fastq files [%(default)s]")

    fastq_group.add_argument("--fastq_read_length", metavar="<INT>", type=int, default=150,
        help='''read length of fastq reads, used by mitoAssemble. [%(default)s]''')


    assembly_group = parser.add_argument_group('Assembly arguments')

    assembly_group.add_argument("--assembler", default='megahit',
        choices=['mitoassemble', 'spades', 'megahit'],
        help="Assembler to be used. [%(default)s]")

    assembly_group.add_argument("--tmp_dir", metavar="<STR>", 
        help="Set temp directory for megahit if necessary (See https://github.com/linzhi2013/MitoZ/issues/176)")

    assembly_group.add_argument("--kmers", type=int,
        metavar="<INT>", default=[71], nargs="+",
        help='kmer size(s) to be used. Multiple kmers can be used, separated by space. Only for mitoassemble %(default)s')

    assembly_group.add_argument("--kmers_megahit",
        metavar="<INT>", default=['21', '29', '39', '59', '79', '99', '119', '141'], nargs="+",
        help="kmer size(s) to be used. Multiple kmers can be used, separated by space. Only for megahit [21 29 39 59 79 99 119 141]")

    assembly_group.add_argument("--kmers_spades",
        metavar="<INT>", default=['auto'], nargs="+",
        help='kmer size(s) to be used. Multiple kmers can be used, separated by space. Only for spades %(default)s')

    assembly_group.add_argument("--memory", metavar="<INT>", default='50',
        help="memory size limit for spades/megahit, no enough memory will make the two programs halt or exit [%(default)s]")

    assembly_group.add_argument("--resume_assembly", action='store_true',
        help="to resume previous assembly running [%(default)s]")

    search_mito_group = parser.add_argument_group('Search' +\
        ' mitochondrial sequences arguments')

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

    search_mito_group.add_argument('--min_abundance', metavar='<float>',
        type=float, default=10,
        help='''the minimum abundance of sequence required. 
        Set this to any value <= 0 if you do NOT want to filter sequences by abundance [%(default)s]''')

    search_mito_group.add_argument("--abundance_pattern", metavar="<STR>",
        default=r'abun\=([0-9]+\.*[0-9]*)',
        help='''the regular expression pattern to capture the abundance information
            in the header of sequence ['%(default)s']''')

    # search_mito_group.add_argument("--skip_read_mapping", action='store_true',
    #    help='''Skip read-mapping step, assuming we can extract the abundance from seqid line. [%(default)s]''')

    search_mito_group.add_argument("--genetic_code", metavar="<INT>",
        default="auto",
        help='''which genetic code table to use? 'auto' means determined by
            '--clade' option. [%(default)s]''')

    search_mito_group.add_argument("--clade", default="Arthropoda",
        choices=["Chordata", "Arthropoda",
            "Echinodermata", "Annelida-segmented-worms",
            "Bryozoa", "Mollusca", "Nematoda",
            "Nemertea-ribbon-worms",
            "Porifera-sponges"],
        help="which clade does your species belong to? [%(default)s]")

    return parser



def use_mitoAssemble(args=None, workdir=None, outprefix=None, thread_number=None, fq1=None, fq2=None, insert_size=None, fastq_read_length=150, kmer=71, resume=True, logger=None):
    """
    de novo with soaptrans, then find mitosequences

    """
    workdir = os.path.abspath(workdir)
    logger.info('use_mitoAssemble() workdir: ' + workdir)

    mitoAssemble_dir = os.path.join(assemble_script_dir, 'mitoAssemble')
    mitoAssemble = os.path.join(mitoAssemble_dir, 'mitoAssemble')
    assembly_dir = os.path.join(workdir, 'mitoAssemble', 'K{kmer}'.format(kmer=kmer))

    logger.info('use_mitoAssemble() mitoAssemble_dir: ' + mitoAssemble_dir)
    logger.info('use_mitoAssemble() mitoAssemble: ' + mitoAssemble)
    logger.info('use_mitoAssemble() assembly_dir: ' + assembly_dir)

    prefix= outprefix + ".mitoAssemble.K{kmer}".format(kmer=kmer)
    assembly_file = os.path.join(assembly_dir, prefix + ".scafSeq")
    logger.info('use_mitoAssemble() prefix: ' + prefix)
    logger.info('use_mitoAssemble() expected assembly_file: ' + assembly_file)

    if not os.path.exists(assembly_dir):
        logger.info('use_mitoAssemble() creating: ' + assembly_dir)
        os.makedirs(assembly_dir, exist_ok=True)

    skip_this_kmer_assembly = False
    if resume and os.path.exists(assembly_file):
        logger.info('use_mitoAssemble() skip_this_kmer_assembly: ' + assembly_dir)
        skip_this_kmer_assembly = True

    logger.info('use_mitoAssemble() chdir: ' + assembly_dir)
    os.chdir(assembly_dir)

    if not skip_this_kmer_assembly:
        # print("\n\n", '-'*30, file=sys.stdout)
        logger.info("use_mitoAssemble() starting assembly with mitoAssemble and kmer {kmer}".format(kmer=kmer))
        trans_lib_f =prefix+".mitoAssemble.lib"
        logger.info("mitoAssemble configure file: " + trans_lib_f)
        fh_tmp = open(trans_lib_f, 'w')
        tmp_out = ""
        tmp_out = "max_rd_len=%s\n" % fastq_read_length + \
            "[LIB]\n" + \
            "avg_ins=%s\n" % insert_size + \
            "reverse_seq=0\n" + \
            "asm_flags=3\n" + \
            "map_len=32\n"

        # pair end data
        if fq1 and fq2 and os.path.isfile(fq1) and os.path.isfile(fq2):
            tmp_out += "q1=" + fq1 + "\n" + "q2=" + fq2
        elif fq1 and os.path.isfile(fq1):
            # single end data
            tmp_out += "q=" + fq1
        elif fq2 and os.path.isfile(fq2):
            # single end data
            tmp_out += "q=" + fq2
        else:
            logger.error('must set at least one of the fq1 and fq2!')

        print(tmp_out, file=fh_tmp)
        fh_tmp.close()

        command = mitoAssemble + " all " +\
                    " -K " + str(kmer) +\
                    " -o " + prefix +\
                    " -s " + trans_lib_f +\
                    " -p " + str(thread_number)

        runcmd(command, logger=logger)

        files_to_del = ".Arc .ContigIndex .PEreadOnContig.gz" +\
                    " .agp .contig .ctg2Read .edge.gz .gapSeq .kmerFreq" +\
                    " .newContigIndex .peGrads .preArc .preGraphBasic" +\
                    " .readInGap .readInformation .readOnContig .readOnScaf" +\
                    " .scaf_gap .shortreadInGap.gz .updated.edge" +\
                    " .vertex"

        command = pre_del_cmd(prefix=prefix, filestr=files_to_del)
        runcmd(command, logger=logger)
    else:
        logger.warning("Using existing assembly file: " + assembly_file)

    if os.path.exists(assembly_file):
        assembly_file_reformated = assembly_file + '.reformatted.fa'
        soft = os.path.join(mitoAssemble_dir, 'reformat_mitoAssemble_scafSeq.py')
        command = python3 + " " + soft +\
            ' {0} '.format(assembly_file) +\
            ' {0} '.format(assembly_file_reformated)
        runcmd(command, logger=logger)

        args.outprefix = prefix
        args.fastafile = assembly_file_reformated
        findmitoscaf_args=copy.copy(args)

        # m = """use_mitoAssemble():\nSince we are using mitoAssemble, there is no need to calculate
        # sequencing abundance. So we set args.fq1 = '' and args.fq2 = ''"""
        # logger.warning(m)
        # args.fq1 = ''
        # args.fq2 = ''
        
        # Check if we're in iterative mode
        if hasattr(args, 'iter') and args.iter > 1:
            # In iterative mode, don't call findmitoscaf here
            # Return the assembly file for iterative processing
            logger.info('use_mitoAssemble() in iterative mode, skipping findmitoscaf')
            logger.info('use_mitoAssemble() found assembly_file: {assembly_file}'.format(assembly_file=assembly_file))
            os.chdir('../../')
            return assembly_file, None
        else:
            # Normal mode: call findmitoscaf
            mt_file = findmitoscaf.main(args)
            os.chdir('../../')
            logger.info('use_mitoAssemble() found:\nassembly_file: {assembly_file}\nmt_file: {mt_file}'.format(assembly_file=assembly_file, mt_file=mt_file))
            return assembly_file, mt_file

    else:
        os.chdir('../../')
        logger.warning("Can not find assembly_file: " + assembly_file)
        
        return None, None


def use_spades(args=None, workdir=None, outprefix=None, thread_number=None, fq1=None, fq2=None, kmers='auto', resume=False, memory=50, logger=None):
    '''
    Currently metaSPAdes supports only a single short-read
    library which has to be paired-end
    '''
    if not (fq1 and fq2):
        logger.error("You have to provide SPAdes with both fq1 and fq2!")
        sys.exit()
    spades_dir = os.path.join(assemble_script_dir, 'spades')
    prefix= outprefix + ".metaspades"

    workdir = os.path.abspath(workdir)
    logger.info('use_spades() workdir: ' + workdir)
    assembly_dir = os.path.join(workdir, 'metaspades')
    if not os.path.exists(assembly_dir):
        logger.info('use_spades() creating: ' + assembly_dir)
        os.makedirs(assembly_dir, exist_ok=True)

    os.chdir(assembly_dir)

    assembly_file = os.path.join(assembly_dir, "scaffolds.fasta")
    logger.info('use_spades() expected assembly_file: ' + assembly_file)

    if resume:
        command = "spades.py --continue -o " + assembly_dir
        runcmd(command, logger=logger)

    else:
        command = "spades.py " +\
            " --checkpoints all " +\
            " -o " + assembly_dir +\
            " --meta " +\
            " --pe1-1 " + fq1 +\
            " --pe1-2 " + fq2 +\
            " --threads " + thread_number +\
            " -k " + kmers +\
            " --memory " + memory
        # Problem: package spades-3.15.4-h95f258a_0 requires python <=3.9, but none of the providers can be installed
        # == Error ==  you cannot specify --careful, --mismatch-correction or --cov-cutoff in metagenomic mode!
        runcmd(command, logger=logger)

    if os.path.exists(assembly_file):
        if resume:
            logger.info("The 'spades.py --continue -o {assembly_dir}' looks good because we found the assembly_file {assembly_file}".format(assembly_dir=assembly_dir, assembly_file=assembly_file))

        assembly_file_reformated = assembly_file + '.reformatted.fa'
        soft = os.path.join(spades_dir, 'reformat_spades_scafSeq.py')
        command = python3 + " " + soft +\
            ' {0} '.format(assembly_file) +\
            ' {0} '.format(assembly_file_reformated)
        runcmd(command, logger=logger)

        args.outprefix = prefix
        args.fastafile = assembly_file_reformated
        findmitoscaf_args=copy.copy(args)

        # Check if we're in iterative mode
        if hasattr(args, 'iter') and args.iter > 1:
            # In iterative mode, don't call findmitoscaf here
            # Return the assembly file for iterative processing
            logger.info('use_spades() in iterative mode, skipping findmitoscaf')
            logger.info('use_spades() found assembly_file: {assembly_file}'.format(assembly_file=assembly_file))
            os.chdir('../')
            return assembly_file, None
        else:
            # Normal mode: call findmitoscaf
            mt_file = findmitoscaf.main(args)
            os.chdir('../')
            logger.info('use_spades() found:\nassembly_file: {assembly_file}\nmt_file: {mt_file}'.format(assembly_file=assembly_file, mt_file=mt_file))
            return assembly_file, mt_file

    else:
        os.chdir('../')
        logger.warning("Can not find assembly_file: " + assembly_file)
        
        return None, None



def use_megahit(args=None, workdir=None, outprefix=None, thread_number=None, fq1=None, fq2=None, kmers='21,29,39,59,79,99,119,141', resume=False, memory=50, tmp_dir=None, logger=None):
    '''
    '''
    megahit_dir = os.path.join(assemble_script_dir, 'megahit')
    prefix= outprefix + ".megahit"

    workdir = os.path.abspath(workdir)
    logger.info('use_megahit() workdir: ' + workdir)
    assembly_dir = os.path.join(workdir, 'megahit')

    if not os.path.exists(assembly_dir):
        logger.info('use_megahit() creating: ' + assembly_dir)
        os.makedirs(assembly_dir, exist_ok=True)

    assembly_file = os.path.join(assembly_dir, "megahit_out", "final.contigs.fa")
    logger.info('use_megahit() expected assembly_file: ' + assembly_file)

    # --continue: continue a MEGAHIT run from its last available check point.
    # please set the output directory correctly when using this option.
    
    os.chdir(assembly_dir)

    if resume:
        command = "megahit --continue --out-dir ./megahit_out "
        if tmp_dir:
            command += " --tmp-dir " + tmp_dir
        runcmd(command, logger=logger)
    else:
        if fq1 and fq2:
            command = "megahit " +\
                " --out-dir ./megahit_out "+\
                " --num-cpu-threads " + thread_number +\
                " --k-list " + kmers +\
                " --memory " + str(memory) +\
                " -1 " + fq1 +\
                " -2 " + fq2
        else:
            command = "megahit " +\
                " --out-dir ./megahit_out " +\
                " --num-cpu-threads " + thread_number +\
                " --k-list " + kmers +\
                " --memory " + str(memory) +\
                " --read " + fq1
        
        if tmp_dir:
            command += " --tmp-dir " + tmp_dir

        runcmd(command, logger=logger)

    if os.path.exists(assembly_file):
        if resume:
            logger.info("The 'megahit --continue --out-dir ./megahit_out' looks good because we found the assembly_file {assembly_file}".format(assembly_file=assembly_file))
        # header meaning: https://github.com/voutcn/megahit/issues/54
        # 'flag' is a tag to represent the connectivity of a contig in the assembly graph. 
        # 'flag=1' means the contig is standalone, 
        # 'flag=2' a looped path and 'flag=0' for other contigs.
        # 'multi' is roughly the average kmer coverage. 
        # But the figure is not precise and if you want to quantify
        # the coverage of a contig, I would suggest you align the
        # reads back to the contigs by short read aligners.
        # 
        # Basically, you could just ignore the header and put the contigs to the subsequent analysis.

        assembly_file_reformated = assembly_file + '.reformatted.fa'
        soft = os.path.join(megahit_dir, 'reformat_megahit_scafSeq.py')
        command = python3 + " " + soft +\
            ' {0} '.format(assembly_file) +\
            ' {0} '.format(assembly_file_reformated)
        runcmd(command, logger=logger)

        args.outprefix = prefix
        args.fastafile = assembly_file_reformated
        findmitoscaf_args=copy.copy(args)

        # Check if we're in iterative mode
        if hasattr(args, 'iter') and args.iter > 1:
            # In iterative mode, don't call findmitoscaf here
            # Return the assembly file for iterative processing
            logger.info('use_megahit() in iterative mode, skipping findmitoscaf')
            logger.info('use_megahit() found assembly_file: {assembly_file}'.format(assembly_file=assembly_file))
            return assembly_file, None
        else:
            # Normal mode: call findmitoscaf
            mt_file = findmitoscaf.main(args)
            logger.info('use_megahit() found:\nassembly_file: {assembly_file}\nmt_file: {mt_file}'.format(assembly_file=assembly_file, mt_file=mt_file))
            return assembly_file, mt_file

    else:
        os.chdir('../')
        logger.warning("Can not find assembly_file: " + assembly_file)
        
        return None, None

def main(args):
    args.skip_read_mapping = True
    logger = args.logger
    logger.info("assemble.main() got args:\n{}".format(args))
    ## genetic code selection
    genetic_code_dict = {"Chordata":"2", "Arthropoda":"5", "Echinodermata":"9",
        "Annelida-segmented-worms":"5", "Bryozoa":"5", "Mollusca":"5",
        "Nematoda":"5", "Nemertea-ribbon-worms":"5", "Porifera-sponges":"4"}

    if hasattr(args, "genetic_code") and args.genetic_code == "auto":
        args.genetic_code = genetic_code_dict[args.clade]

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

    logger.info("assemble.main() got updated args:\n{}".format(args))

    fq1 = args.fq1
    fq2 = args.fq2

    logger.info("assemble.main() got\nfq1: {0}\nfq2: {1}".format(fq1, fq2))

    insert_size = args.insert_size
    fastq_read_length = args.fastq_read_length

    assembler = args.assembler
    thread_number = args.thread_number
    kmers = args.kmers
    kmers_spades = ' '.join(args.kmers_spades)
    kmers_megahit = ','.join(args.kmers_megahit)
    memory = args.memory
    memory_in_bytes = int(memory) * 10e9 
    resume_assembly = args.resume_assembly
    outprefix = args.outprefix
    workdir = os.path.abspath(args.workdir)
    tmp_dir = args.tmp_dir
    logger.info("assemble.main() working in {}".format(workdir))

    resulting_mt_files = []

    if assembler == 'mitoassemble':
        for kmer in kmers:
            kmer = int(kmer)
            if kmer > int(fastq_read_length):
                logger.warning("Kmer {kmer} is larger than the fastq_read_length {fastq_read_length}!Skip over this kmer assembly!".format(kmer=kmer, fastq_read_length=fastq_read_length))
            assembly_file, mt_file = use_mitoAssemble(
                args=copy.copy(args),
                workdir=workdir,
                outprefix=outprefix,
                thread_number=thread_number,
                fq1=fq1,
                fq2=fq2,
                insert_size=insert_size,
                fastq_read_length=fastq_read_length,
                kmer=kmer,
                resume=True,
                logger=logger)

            resulting_mt_files.append(mt_file)

    elif assembler == 'spades':
        assembly_file, mt_file = use_spades(
            args=copy.copy(args),
            workdir=workdir,
            outprefix=outprefix,
            thread_number=thread_number,
            fq1=fq1,
            fq2=fq2,
            kmers=kmers_spades,
            memory=memory,
            resume=resume_assembly,
            logger=logger)
        resulting_mt_files.append(mt_file)

    elif assembler == 'megahit':
        assembly_file, mt_file = use_megahit(
            args=copy.copy(args),
            workdir=workdir,
            outprefix=outprefix,
            thread_number=thread_number,
            fq1=fq1,
            fq2=fq2,
            kmers=kmers_megahit,
            resume=resume_assembly,
            memory=memory_in_bytes,
            tmp_dir=tmp_dir,
            logger=logger)
        resulting_mt_files.append(mt_file)

    # Filter out None values (which occur in iterative mode)
    resulting_mt_files = [f for f in resulting_mt_files if f is not None]
    
    if resulting_mt_files:
        logger.info('assemble.main() found resulting_mt_files:\n' + '\n'.join(resulting_mt_files))
    else:
        logger.info('assemble.main() found no mitogenome files (iterative mode)')

    all_result_wdir = []
    for f in resulting_mt_files:
        if f:  # Check if f is not None
            d = find_subdirs_with_suffix(indir=os.path.dirname(f), suffix='.result')
            if len(d) > 0:
                all_result_wdir.extend(d)

    d = find_subdirs_with_suffix(indir=workdir, suffix='.result')
    for f in d:
        if f not in all_result_wdir:
            all_result_wdir.append(f)

    if all_result_wdir:
        logger.info('And related resulting files can be found from:\n' + '\n'.join(all_result_wdir))

    return resulting_mt_files, all_result_wdir

