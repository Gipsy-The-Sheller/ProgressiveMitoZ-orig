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


# https://www.gnu.org/licenses/gpl-faq.zh-cn.html#GPLPlugins


import sys
import os
import re
import subprocess
import time
from glob import glob

import Bio
from Bio import SeqIO
from ete3 import NCBITaxa
from mitoz.utility.utility import gather_result, runcmd, files_exist_0_or_1, file_not_empty
from pathlib import Path


# from mitoz import findmitoscaf
# from mitoz import annotate
# from mitoz import profiles

python3 = sys.executable

mitoz_pkg_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
findmitoscaf_script_dir = os.path.join(mitoz_pkg_dir, 'findmitoscaf', 'script')
annotate_script_dir = os.path.join(mitoz_pkg_dir, 'annotate', 'script')
profiles_dir = os.path.join(mitoz_pkg_dir, 'profiles')

def add_arguments(parser):
    parser.add_argument('--fastafile', metavar="<file>", required=True, help="Input fasta file. Gzip supported. [required]")

    parser.add_argument("--fq1", metavar="<file>", default='',
        help='''Input fastq 1 file. use this option if the headers of your '--fastafile' 
        does NOT have abundance information BUT you WANT to filter
         sequence by their sequencing abundances [optional]''')

    parser.add_argument("--fq2", metavar="<file>", default='',
        help='''Input fastq 2 file. use this option if the headers of your '--fastafile' 
        does NOT have abundance information BUT you WANT to filter
         sequence by their sequencing abundances [optional]''')

    parser.add_argument("--outprefix", metavar="<STR>", required=True,
        help="output prefix")

    parser.add_argument("--workdir", metavar="<STR>", default='./',
        help="workdir [%(default)s]")

    parser.add_argument("--thread_number", metavar="<INT>", default="8",
        help="thread number [%(default)s]")

    parser.add_argument("--profiles_dir", metavar="<STR>",
        default=profiles_dir,
        help="Directory cotaining 'CDS_HMM/', 'MT_database/' and 'rRNA_CM/'. [%(default)s]")

    parser.add_argument("--slow_search", action='store_true',
        help='''By default, we firstly use tiara to perform quick sequence classification (100 times faster than usual!),
however, it is valid only when your mitochondrial sequences are >= 3000 bp.
If you have missing genes, set '--slow_search' to use the tradicitiona search mode. [%(default)s]''')

    parser.add_argument("--filter_by_taxa", action='store_false',
        help='''filter out non-requiring_taxa sequences by mito-PCGs annotation
    to do taxa assignment.[%(default)s]''')

    parser.add_argument("--requiring_taxa", metavar="<STR>",
        required=True,
        help='''filtering out non-requiring taxa sequences which may be
    contamination [required]''')

    parser.add_argument("--requiring_relax", default="0",
        choices=["0", "1", "2", "3", "4", "5", "6"],
        help='''The relaxing threshold for filtering non-target-requiring_taxa.
    The larger digital means more relaxing. [%(default)s]''')

    parser.add_argument('--min_abundance', metavar='<float>',
        type=float, default=10,
        help='''the minimum abundance of sequence required. 
        Set this to any value <= 0 if you do NOT want to filter sequences by abundance [%(default)s]''')

    parser.add_argument("--abundance_pattern", metavar="<STR>",
        default=r'abun\=([0-9]+\.*[0-9]*)',
        help='''the regular expression pattern to capture the abundance information
            in the header of sequence ['%(default)s']''')

    parser.add_argument("--skip_read_mapping", action='store_true',
        help='''Skip read-mapping step, assuming we can extract the abundance from seqid line. [%(default)s]''')

    parser.add_argument("--genetic_code", metavar="<INT>",
        default="auto",
        help='''which genetic code table to use? 'auto' means determined by
            '--clade' option. [%(default)s]''')

    parser.add_argument("--clade", default="Arthropoda",
        choices=["Chordata", "Arthropoda",
            "Echinodermata", "Annelida-segmented-worms",
            "Bryozoa", "Mollusca", "Nematoda",
            "Nemertea-ribbon-worms",
            "Porifera-sponges"],
        help="which clade does your species belong to? [%(default)s]")


    return parser


def filter_by_tiara(assembly_file=None, thread_number=None, prefix='tiara', tiara='tiara', prob_cutoff="0.65 0.60", min_len=1000, logger=None):
    mt_expected_result_file = "mitochondrion_" + os.path.basename(assembly_file)
    logger.info("tiara mt_expected_result_file: " + mt_expected_result_file)
    command = "{tiara} --input {assembly_file} --output {prefix}.tiara.out.txt ".format(tiara=tiara, assembly_file=assembly_file, prefix=prefix)
    command += " --to_fasta mit pla pro --threads {thread_number} ".format(thread_number=thread_number)
    command += " --prob_cutoff {prob_cutoff} --probabilities --min_len {min_len}".format(prob_cutoff=prob_cutoff, min_len=min_len)
    runcmd(command, logger=logger)

    return mt_expected_result_file


def filter_by_nhmmer(assembly_file=None, thread_number=None, prefix=None, nhmmer='nhmmer', nhmmer_profile=None, logger=None):
    # nhmmer needs uncompress input fasta file when multiple hmm profiles
    # in one file
    if assembly_file.endswith(".gz"):
        hmm_infile = prefix + ".tmp.nhmmer.infile"
        command = "gzip -dc " + assembly_file + " > " + hmm_infile
        runcmd(command, logger=logger)
    else:
        hmm_infile = assembly_file

    ## do nhmmer
    hmm_outf = prefix + ".hmmout"
    hmm_outf_tbl = prefix + ".hmmtblout"
    command = nhmmer + " -o " + hmm_outf + " --tblout " + hmm_outf_tbl +\
            " --cpu " + str(thread_number) + " " + nhmmer_profile +\
            " " + hmm_infile
    runcmd(command, logger=logger)


    ## get besthit of HMM genes for each sequence
    hmm_outf_tbl_besthit = prefix + ".hmmtblout.besthit"
    soft = os.path.join(findmitoscaf_script_dir,
        "get_besthit_of_each_CDS_from_nhmmer.py")
    command = python3 + " " + soft +\
        " " + hmm_outf_tbl +\
        " " + hmm_outf_tbl_besthit
    runcmd(command, logger=logger)

    ## discard some fields of nhmmer result
    hmm_outf_tbl_besthit_sim = prefix + ".hmmtblout.besthit.sim"
    soft = os.path.join(findmitoscaf_script_dir, "simlify_nhmmer_tbl_besthit.py")
    command = python3 + " " + soft +\
        " " + hmm_outf_tbl_besthit +\
        " " + hmm_outf_tbl_besthit_sim
    runcmd(command, logger=logger)

    # extract all hmm besthit sequences
    hmm_besthit_fa = hmm_outf_tbl_besthit_sim + ".fa"
    soft = os.path.join(findmitoscaf_script_dir, "extract_fasta.py")
    command = python3 + " " + soft +\
        " -i " + hmm_infile +\
        " -q " + hmm_outf_tbl_besthit_sim +\
        " -o " + hmm_besthit_fa
    runcmd(command, logger=logger)

    logger.info("filter_by_nhmmer() return:\n{hmm_infile}, {hmm_outf_tbl_besthit_sim}, {hmm_besthit_fa}".format(hmm_infile=hmm_infile, hmm_outf_tbl_besthit_sim=hmm_outf_tbl_besthit_sim, hmm_besthit_fa=hmm_besthit_fa))
    return hmm_infile, hmm_outf_tbl_besthit_sim, hmm_besthit_fa


def do_filter_by_taxa(hmm_besthit_fa=None, hmm_outf_tbl_besthit_sim=None, MT_annotation_BGI_pro_db=None, thread_number=None, require_genetic_code=None, requiring_taxa=None, requiring_relax=None, logger=None):
    MT_annotation_BGI = os.path.join(annotate_script_dir, 'MT_annotation_BGI_V1.32', 'MT_annotation_BGI.pl')
    # WISECONFIGDIR = os.path.join(annotate_script_dir, "wisecfg")
    command = 'which genewise'
    genewise_binary = subprocess.check_output(command, shell=True, encoding='utf8')
    genewise_binary = genewise_binary.strip()
    env_dir = os.path.dirname(os.path.dirname(genewise_binary))
    WISECONFIGDIR = os.path.join(env_dir, 'share/wise2/wisecfg')
    
    # 2. filter taxonomy by MT PCGs annotation
    hmm_besthit_filtered_by_taxa_fa = hmm_outf_tbl_besthit_sim + ".filtered-by-taxa.fa"
    soft = os.path.join(findmitoscaf_script_dir,
        "filter_taxonomy_by_CDS_annotation.py")
    command = python3 + " " + soft +\
        " -fa " + hmm_besthit_fa +\
        " -MTsoft " + MT_annotation_BGI +\
        " -db " + MT_annotation_BGI_pro_db +\
        " -thread " + str(thread_number) +\
        " -genetic_code " + require_genetic_code +\
        " -requiring_taxa " + "'" + requiring_taxa + "'" +\
        " -relax " + requiring_relax +\
        " -WISECONFIGDIR " + WISECONFIGDIR +\
        " -outf " + hmm_besthit_filtered_by_taxa_fa
    runcmd(command, logger=logger)

    # 3. filter the 'hmm_outf_tbl_besthit_sim' file
    hmm_outf_tbl_besthit_sim_filtered_by_taxa = \
        hmm_outf_tbl_besthit_sim + ".filtered-by-taxa"
    soft = os.path.join(findmitoscaf_script_dir, "filter_hmm-besthit-sim.py")
    command = python3 + " " + soft +\
        " " + hmm_besthit_filtered_by_taxa_fa +\
        " " + hmm_outf_tbl_besthit_sim +\
        " " + hmm_outf_tbl_besthit_sim_filtered_by_taxa
    runcmd(command, logger=logger)

    logger.info("do_filter_by_taxa() return:\nhmm_outf_tbl_besthit_sim_filtered_by_taxa:{hmm_outf_tbl_besthit_sim_filtered_by_taxa}\nhmm_besthit_filtered_by_taxa_fa: {hmm_besthit_filtered_by_taxa_fa}".format(hmm_outf_tbl_besthit_sim_filtered_by_taxa=hmm_outf_tbl_besthit_sim_filtered_by_taxa, hmm_besthit_filtered_by_taxa_fa=hmm_besthit_filtered_by_taxa_fa))

    return hmm_outf_tbl_besthit_sim_filtered_by_taxa, hmm_besthit_filtered_by_taxa_fa


def filter_by_abundance(assembly_file=None, fq1=None, fq2=None, skip_read_mapping=True, abundance_pattern=None, min_abundance=None, hmm_outf_tbl_besthit_sim=None, thread_number=8, result_wdir=None, logger=None):
    if os.path.exists(fq1) and not skip_read_mapping:
        logger.info("We are mapping reads to the assembly file to calculate sequence abundance...")
        soft = os.path.join(findmitoscaf_script_dir, "cal_bwa_abundance.py")
        abun_out_file = os.path.basename(assembly_file) + ".abun.fa"
        command = python3 + " " + soft +\
                " -fa " + assembly_file +\
                " -fq1 " + fq1 +\
                " -fq2 '{}' ".format(fq2) +\
                " -out " + abun_out_file +\
                " -bwa bwa " +\
                " -samtools samtools" +\
                " -thread " + thread_number

        runcmd(command, logger=logger)
        assembly_file = abun_out_file

        # append the abundance to the file 'hmm_outf_tbl_besthit_sim'
        hmm_outf_tbl_besthit_sim_new_abun = hmm_outf_tbl_besthit_sim + '.new_abun'
        # get abundance dict
        abun_dict = {}
        with open(abun_out_file, 'r') as fh:
            for i in fh:
                i = i.strip()
                if not i.startswith('>'):
                    continue
                seqid, seq_info = i.split(maxsplit=1)
                seqid = re.sub('^>', '', seqid)
                seq_info = re.sub(r'\t+', ' ', seq_info)
                abun_dict[seqid] = seq_info
        with open(hmm_outf_tbl_besthit_sim, 'r') as fh_in, open(hmm_outf_tbl_besthit_sim_new_abun, 'w') as fh_out:
            for i in fh_in:
                i = i.strip()
                if not i:
                    continue
                seqid = i.split()[0]
                seq_info = abun_dict[seqid]
                print(i, seq_info, sep='\t', file=fh_out)

        hmm_outf_tbl_besthit_sim = hmm_outf_tbl_besthit_sim_new_abun


    ## filter sequences with abundance < 10X
    # high_abun_sim = '{0}.high_abundance_{1}X'.format(hmm_outf_tbl_besthit_sim, min_abundance)
    soft = os.path.join(findmitoscaf_script_dir, 'filter_by_abundance.py')
    command = python3 + " " + soft +\
        " '{0}' ".format(abundance_pattern) +\
        ' {0} '.format(min_abundance) +\
        ' {0} '.format(hmm_outf_tbl_besthit_sim) +\
        ' {0} '.format(assembly_file)
    runcmd(command, logger=logger)
    
    # and gather the low abundance sequences
    low_abundance_files = glob('*low_abundance*')
    gather_result(*low_abundance_files, result_wdir=result_wdir, logger=logger)

    # must after running the above command, check the result file
    high_abun_sim = glob('{0}.high_abundance_*'.format(hmm_outf_tbl_besthit_sim))[0]
    if file_not_empty(high_abun_sim):
        ## reformat to fasta-like format
        hmm_outf_tbl_besthit_sim_fasta_like = hmm_outf_tbl_besthit_sim + ".fasta-like"
        soft = os.path.join(findmitoscaf_script_dir, "reformat_nhmmer_besthit-sim.py")
        command = python3 + " " + soft +\
                " " + high_abun_sim +\
                " " + hmm_outf_tbl_besthit_sim_fasta_like
        runcmd(command, logger=logger)

        logger.info("filter_by_abundance() returns hmm_outf_tbl_besthit_sim_fasta_like: {}".format(hmm_outf_tbl_besthit_sim_fasta_like))
        return hmm_outf_tbl_besthit_sim_fasta_like

    else:
        logger.warning('filter_by_abundance():\nAll sequences are low abundance (<{0}X)'.format(min_abundance))
        return None


def select_sequences_by_gene_completeness(defined_gene_length_f=None, hmm_outf_tbl_besthit_sim_fasta_like=None, hmm_infile=None, result_wdir=None, prefix=None, abundance_pattern=None, logger=None):
    ## scoring and sorting
    hmm_outf_tbl_besthit_sim_fasta_like_sorted = \
        hmm_outf_tbl_besthit_sim_fasta_like + ".sorted"
    soft = os.path.join(findmitoscaf_script_dir,
                        "scoring_nhmmer_besthit_sim_reformat.py")
    command = python3 + " " + soft +\
            " " + defined_gene_length_f +\
            " " + hmm_outf_tbl_besthit_sim_fasta_like +\
            " '{0}' ".format(abundance_pattern) +\
            " " + hmm_outf_tbl_besthit_sim_fasta_like_sorted
    runcmd(command, logger=logger)

    ## picking by dotting method
    hmm_outf_tbl_besthit_sim_fasta_like_sorted_dotted = \
        hmm_outf_tbl_besthit_sim_fasta_like + ".sorted.dotted"
    hmm_outf_tbl_besthit_sim_fasta_like_sorted_picked = \
        hmm_outf_tbl_besthit_sim_fasta_like + ".sorted.picked"
    hmm_outf_tbl_besthit_sim_fasta_like_sorted_picked_stat = \
        hmm_outf_tbl_besthit_sim_fasta_like + ".sorted.picked.stat"
    soft = os.path.join(findmitoscaf_script_dir,
        "pick_seq_from_scoing_result_base_arr_v4.py")
    command = python3 + " " + soft +\
            " " + hmm_outf_tbl_besthit_sim_fasta_like_sorted +\
            " " + hmm_outf_tbl_besthit_sim_fasta_like_sorted_dotted +\
            " " + hmm_outf_tbl_besthit_sim_fasta_like_sorted_picked +\
            " " + hmm_outf_tbl_besthit_sim_fasta_like_sorted_picked_stat
    runcmd(command, logger=logger)

    ## extract picked fasta sequences
    soft = os.path.join(findmitoscaf_script_dir, "extract_nhmmer_seq_v2.py")
    selected_mt_fa = hmm_outf_tbl_besthit_sim_fasta_like_sorted_picked + ".fa"
    command = python3 + " " + soft +\
            " -f " + hmm_outf_tbl_besthit_sim_fasta_like_sorted_picked +\
            " -d " + hmm_infile +\
            " -o " + selected_mt_fa
    runcmd(command, logger=logger)

    # check result file
    if not file_not_empty(selected_mt_fa):
        logger.warning("select_sequences_by_gene_completeness():\ncan not find any mitosequences by nhmmer!")
        return None

    ## extract the not-picked fasta sequences
    soft = os.path.join(findmitoscaf_script_dir, "extract_Not-picked-seq.py")
    hmm_outf_tbl_besthit_sim_fasta_like_sorted_Not_picked = \
        hmm_outf_tbl_besthit_sim_fasta_like + ".sorted.Not-picked"
    hmm_outf_tbl_besthit_sim_fasta_like_sorted_Not_picked_fa = hmm_outf_tbl_besthit_sim_fasta_like_sorted_Not_picked + '.fa'
    command = python3 + " " + soft +\
        " " + hmm_outf_tbl_besthit_sim_fasta_like_sorted +\
        " " + hmm_outf_tbl_besthit_sim_fasta_like_sorted_picked +\
        " " + hmm_infile +\
        " " + hmm_outf_tbl_besthit_sim_fasta_like_sorted_Not_picked +\
        " " + hmm_outf_tbl_besthit_sim_fasta_like_sorted_Not_picked_fa
    runcmd(command, logger=logger)

    gather_result(hmm_outf_tbl_besthit_sim_fasta_like_sorted_Not_picked, hmm_outf_tbl_besthit_sim_fasta_like_sorted_Not_picked_fa, result_wdir=result_wdir, logger=logger)

    logger.info("select_sequences_by_gene_completeness() return:\n" + os.path.abspath(selected_mt_fa))

    return os.path.abspath(selected_mt_fa)


def quick_circle_check(mitoscaf_file=None, prefix='ZZZ', logger=None):
    """
    check if sequences are circular (seq length >= 12Kbp), do not map reads
    to mitoscaf with bwa

    """
    mt_file = "{0}.mitogenome.fa".format(prefix)
    overlap_information_file = "{0}.overlap_information".format(prefix)

    soft = os.path.join(findmitoscaf_script_dir, "circle_check.py")
    command = python3 + " " + soft +\
            " " + mitoscaf_file +\
            " " + prefix +\
            " 3"
    runcmd(command, logger=logger)

    logger.info("quick_circle_check() return:\n" + os.path.abspath(mt_file))
    return os.path.abspath(mt_file), os.path.abspath(overlap_information_file)


def main(args):
    """
    find candidate mitochondrial sequences, about 2-3 Gbp fastq data is
    needed to calculate the average sequencing depth of each sequences. Unless
    the input assembly is from SOAPTrans (specify with --from_soaptrans option)

    """

    # step 1: HMMER
    # step 2 (optionally): calculate the abundance
    # step 3: choose the candidate mt sequences

    logger = args.logger
    logger.info("findmitoscaf.main() got args:\n{}".format(args))

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

    logger.info("findmitoscaf.main() got updated args:\n{}".format(args))
    

    skip_read_mapping = args.skip_read_mapping

    # findmitoscaf_script_dir
    # annotate_script_dir
    # profiles_dir

    nhmmer_profile = os.path.join(args.profiles_dir, 'CDS_HMM', args.clade+"_CDS.hmm")
    defined_gene_length_f = os.path.join(args.profiles_dir, 'CDS_HMM', args.clade+"_CDS_length_list")
    # MT_annotation_BGI_pro_db = os.path.join(args.profiles_dir, "MT_database", args.clade +"_CDS_protein.fa") ## for test only (to finish quick)
    MT_annotation_BGI_pro_db = os.path.join(args.profiles_dir, "MT_database", "Animal_CDS_protein.fa") # change to this after test!!

    logger.info("findmitoscaf.main() nhmmer_profile:\n{}".format(nhmmer_profile))
    logger.info("findmitoscaf.main() defined_gene_length_f:\n{}".format(defined_gene_length_f))
    logger.info("findmitoscaf.main() MT_annotation_BGI_pro_db:\n{}".format(MT_annotation_BGI_pro_db))

    workdir = os.path.abspath(args.workdir)
    result_wdir = os.path.join(workdir, args.outprefix + ".result")
    logger.info("findmitoscaf.main() workdir:\n{}".format(workdir))
    logger.info("findmitoscaf.main() result_wdir:\n{}".format(result_wdir))

    if Path(args.fastafile).is_file():
        assembly_file = os.path.abspath(args.fastafile)
    else:
        logger.error("Can NOT access {} !!".format(args.fastafile))
        sys.exit(0)
    # assembly_file = os.path.abspath(args.fastafile)
    logger.info("findmitoscaf.main() got assembly_file:\n{}".format(assembly_file))

    fq1 = args.fq1
    fq2 = args.fq2

    logger.info("findmitoscaf.main() got\nfq1: {0}\nfq2: {1}".format(fq1, fq2))

    slow_search = args.slow_search
    filter_by_taxa = args.filter_by_taxa
    requiring_taxa = args.requiring_taxa
    require_genetic_code = args.genetic_code
    requiring_relax = args.requiring_relax

    thread_number = args.thread_number
    prefix = args.outprefix

    min_abundance = args.min_abundance
    abundance_pattern = args.abundance_pattern
    # test if it works
    if float(min_abundance) > 0 and abundance_pattern and not skip_read_mapping:
        has_fq = fq1 or fq2
        if not has_fq:
            with open(assembly_file, 'r') as fh:
                for i in fh:
                    i = i.strip()
                    if not i:
                        continue
                    if i.startswith('>'):
                        m = re.search(abundance_pattern, i)
                        try:
                            abun = m.group(1)
                        except AttributeError:
                            m = """Failed to extract the abundance information from the seqid!

    abundance_pattern: {abundance_pattern}
    The seqid of your input assembly file is like: {seqid}

    You may need to adapt the regular expression using the '--abundance_pattern' option,

    Or set '--min_abundance 0' if you do NOT want to filter sequences by abundance!
                            """.format(abundance_pattern=abundance_pattern, seqid=i)
                            logger.error(m)
                            sys.exit(0)
                        break
        else:
            logger.info("We will calculate sequence abundances based on the input fastq file!")

    else:
        logger.info("We will skip read-mapping during findmitoscaf().")


    if not slow_search:
        assembly_file = filter_by_tiara(
            assembly_file=assembly_file,
            thread_number=thread_number,
            prefix=prefix,
            tiara='tiara',
            prob_cutoff="0.65 0.60",
            min_len=1000,
            logger=logger)


    hmm_infile, hmm_outf_tbl_besthit_sim, hmm_besthit_fa = filter_by_nhmmer(
        assembly_file=assembly_file,
        thread_number=thread_number,
        prefix=prefix,
        nhmmer_profile=nhmmer_profile,
        logger=logger)


    if filter_by_taxa:
        hmm_outf_tbl_besthit_sim_filtered_by_taxa, hmm_besthit_filtered_by_taxa_fa = do_filter_by_taxa(
            hmm_besthit_fa=hmm_besthit_fa,
            hmm_outf_tbl_besthit_sim = hmm_outf_tbl_besthit_sim,
            MT_annotation_BGI_pro_db=MT_annotation_BGI_pro_db,
            thread_number=thread_number,
            require_genetic_code=require_genetic_code,
            requiring_taxa=requiring_taxa,
            requiring_relax=requiring_relax,
            logger=logger)
   
    if float(min_abundance) > 0:
        if filter_by_taxa:
            hmm_outf_tbl_besthit_sim_fasta_like = filter_by_abundance(
                assembly_file=hmm_besthit_filtered_by_taxa_fa,
                fq1=fq1,
                fq2=fq2,
                skip_read_mapping=skip_read_mapping,
                abundance_pattern=abundance_pattern, 
                min_abundance=min_abundance,
                hmm_outf_tbl_besthit_sim=hmm_outf_tbl_besthit_sim_filtered_by_taxa,
                thread_number=thread_number,
                result_wdir=result_wdir,
                logger=logger)
        else:
            hmm_outf_tbl_besthit_sim_fasta_like = filter_by_abundance(
                assembly_file=hmm_infile,
                fq1=fq1,
                fq2=fq2,
                skip_read_mapping=skip_read_mapping,
                abundance_pattern=abundance_pattern, 
                min_abundance=min_abundance,
                hmm_outf_tbl_besthit_sim=hmm_outf_tbl_besthit_sim,
                thread_number=thread_number,
                result_wdir=result_wdir,
                logger=logger)

        selected_mt_fa = select_sequences_by_gene_completeness(
            defined_gene_length_f=defined_gene_length_f,
            hmm_outf_tbl_besthit_sim_fasta_like=hmm_outf_tbl_besthit_sim_fasta_like,
            hmm_infile=hmm_infile,
            result_wdir=result_wdir,
            abundance_pattern=abundance_pattern,
            logger=logger)

        mt_file, overlap_information_file = quick_circle_check(
            mitoscaf_file=selected_mt_fa,
            prefix=prefix,
            logger=logger)

        gather_result(selected_mt_fa, mt_file, overlap_information_file, result_wdir=result_wdir, logger=logger)
        
        return mt_file

    else:
        if filter_by_taxa:
            mt_file,overlap_information_file = quick_circle_check(
                mitoscaf_file=hmm_besthit_filtered_by_taxa_fa,
                prefix=prefix,
                logger=logger)
            gather_result(hmm_besthit_filtered_by_taxa_fa, mt_file, overlap_information_file, result_wdir=result_wdir, logger=logger)

            return mt_file

        else:
            mt_file, overlap_information_file = quick_circle_check(
                mitoscaf_file=hmm_besthit_fa,
                prefix=prefix,
                logger=logger)

            gather_result(hmm_besthit_fa, mt_file, overlap_information_file, result_wdir=result_wdir, logger=logger)

            return mt_file



