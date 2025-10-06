
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
from mitoz import visualize
from mitoz.utility.utility import gather_result, runcmd, runcmd2, files_exist_0_or_1, file_not_empty
from mitoz.utility.utility import pre_del_cmd, abspath, find_subdirs_with_suffix, check_fafile_seqid_len
import copy
from Bio import SeqIO
from pathlib import Path


python3 = sys.executable

mitoz_pkg_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
findmitoscaf_script_dir = os.path.join(mitoz_pkg_dir, 'findmitoscaf', 'script')
annotate_script_dir = os.path.join(mitoz_pkg_dir, 'annotate', 'script')
profiles_dir = os.path.join(mitoz_pkg_dir, 'profiles')
template_sbt = os.path.join(annotate_script_dir, 'template.sbt')


def add_arguments(parser):
    parser.add_argument("--workdir", metavar="<STR>", default='./',
        help="workdir [%(default)s]")

    parser.add_argument("--outprefix", metavar="<STR>", required=True,
        help="output prefix [required]")

    parser.add_argument("--thread_number", metavar="<INT>", default="8",
        help="thread number [%(default)s]")

    parser.add_argument('--fastafiles', metavar="<STR>", nargs='+', required=True,
        help='''fasta file(s). The length of sequence id should be <= 13 characters, and each sequence should have 'topology=linear' or
         'topology=circular' at the header line, otherwise it is assumbed to be 'topology=linear'. For example, '>Contig1 topology=linear' [required]''')

    parser.add_argument('--fq1', metavar='<file>', default='',
        help="Fastq1 file if you want to visualize the depth distribution")

    parser.add_argument('--fq2', metavar='<file>', default='',
        help="Fastq2 file if you want to visualize the depth distribution")

    parser.add_argument("--profiles_dir", metavar="<STR>",
        default=profiles_dir,
        help="Directory cotaining 'CDS_HMM/', 'MT_database/' and 'rRNA_CM/'. [%(default)s]")

    parser.add_argument("--species_name", metavar="<STR>",
        default="Test sp.",
        help='''species name to use in output genbank file ['%(default)s']''')

    parser.add_argument("--template_sbt", metavar="<file>",
        default=template_sbt,
        help='''The sqn template to generate the resulting genbank file. 
        Go to https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/#Template
        to generate your own template file if you like. ['%(default)s']''')

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

    parser.add_argument("--mitoz_module_used", metavar="<STR>",
        default="annotate",
        help=argparse.SUPPRESS)

    return parser


class Visualize_Args():
    def __init__(self, raw_args):
        return self.prepare_new_args(raw_args)

    def prepare_new_args(self, raw_args):
        self.logger = raw_args.logger
        self.circos = 'circos'
        self.gb = ''
        self.gc = 'no'
        self.win = 50
        self.gc_fill = '128,177,211'
        self.depth_file = ''
        self.run_map = 'yes'

        self.fq1 = ''
        self.fq2 = ''
        if raw_args.fq1 and Path(raw_args.fq1).is_file():
            self.fq1 = os.path.abspath(raw_args.fq1)
        if raw_args.fq2 and Path(raw_args.fq2).is_file():
            self.fq2 = os.path.abspath(raw_args.fq2)

        self.bwa = 'bwa'
        self.thread = raw_args.thread_number
        self.samtools = 'samtools'
        self.opts_samtools = "-a -a"
        self.depth_fill = '190,186,218'
        self.cds_color = '141,211,199'
        self.trna_color = '251,128,114'
        self.rrna_color = '253,192,134'
        self.label_color = 'black'
        self.locus_color = 'black'
        self.base = 'no'
        self.bgc = 'white'
        self.outdir = './'
        self.png = 'yes'
        self.svg = 'yes'

    def __str__(self):
        return str(self.__dict__)


def cds_annotation(workdir=None, mt_file=None, MT_annotation_BGI_pro_db=None, genetic_code=None, thread_number=None, logger=None):
    MT_annotation_BGI = os.path.join(annotate_script_dir, 'MT_annotation_BGI_V1.32', 'MT_annotation_BGI.pl')
    # WISECONFIGDIR = os.path.join(annotate_script_dir, "wisecfg")
    
    mt_file_base = os.path.basename(mt_file)

    current_dir = os.getcwd()

    if not os.path.exists(workdir):
        os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)    
    logger.info('cds_annotation() chdir to ' + workdir)

    command = 'which genewise'
    genewise_binary = subprocess.check_output(command, shell=True, encoding='utf8')
    genewise_binary = genewise_binary.strip()
    env_dir = os.path.dirname(os.path.dirname(genewise_binary))
    WISECONFIGDIR = os.path.join(env_dir, 'share/wise2/wisecfg')

    command = "export WISECONFIGDIR=%s\n" % WISECONFIGDIR
    # must cat export and perl command together here
    soft = MT_annotation_BGI
    command += "perl " + soft +\
            " -i " + mt_file +\
            " -d " + MT_annotation_BGI_pro_db +\
            " -o ./" +\
            " -g " + genetic_code +\
            " -cpu " + thread_number
    runcmd(command, logger=logger)

    command = "rm -rf *.nhr *.nin *.nsq *.ndb *.nto *.not *.ntf formatdb.log *.tblastn.shell *.length"
    command += " *.blast *.blast.filter *.shell *.blast.solar"
    command += " *.blast.solar.filter *.blast.solar.filter.table"
    command += " *.blast.solar.filter.table.nonredundance"
    command += " *.blast.solar.filter.nr *.genewise *.solar.genewise.gff"
    command += " *.solar.genewise.gff.cds"
    runcmd(command, logger=logger)

    soft = os.path.join(annotate_script_dir, 'prepare_CDSposition.py')
    MT_annotation_BGI_outfile = mt_file_base +\
        ".solar.genewise.gff.cds.position.cds"
    cds_position = mt_file_base + ".cds.position"
    command = python3 + " " + soft +\
            " " + MT_annotation_BGI_outfile +\
            " " + cds_position
    runcmd(command, logger=logger)

    # get the most related species information
    soft = os.path.join(annotate_script_dir, 'pick_most_related_sp_from_position_cds.py')
    annt_most_related_sp_file = mt_file_base + '.most_related_species.txt'
    annt_most_related_sp_file = os.path.join(workdir, annt_most_related_sp_file)
    command = python3 + " " + soft +\
            " " + MT_annotation_BGI_outfile +\
            " " + mt_file +\
            " {0} ".format(annt_most_related_sp_file)
    runcmd(command, logger=logger)
    # gather_result(annt_most_related_sp_file)

    cds_position_sorted = cds_position + ".sorted"
    command = "sort -k1,1 -k6,6n " + cds_position + " >" + cds_position_sorted
    runcmd(command, logger=logger)

    ## find positions of start codon and stop codon
    soft = os.path.join(annotate_script_dir, 'revise_CDS_pos_v6.py')
    cds_position_revised = cds_position_sorted + ".revised"
    command = python3 + " " + soft +\
            " " + mt_file +\
            " " + cds_position_sorted +\
            " " + genetic_code +\
            " " + cds_position_revised
    runcmd(command, logger=logger)


    # filter out some cds with quite short annotated regions
    cds_position_revised_filtered = cds_position_revised + ".filtered"
    soft = os.path.join(annotate_script_dir, "filter_cds_by_annotated-len.py")
    command = python3 + " " + soft +\
        " " + cds_position_revised +\
        " " + cds_position_revised_filtered +\
        " 0.8 "
    runcmd(command, logger=logger)

    cds_position_revised = cds_position_revised_filtered

    ## prepare CDS feature table file
    soft = os.path.join(annotate_script_dir, 'cds_ft_v2.py')
    mt_file_cdsft = mt_file_base + ".cds.ft"
    command = python3 + " " + soft +\
            " " + cds_position_revised +\
            " " + genetic_code +\
            " " + mt_file_cdsft
    runcmd(command, logger=logger)

    logger.info('cds_annotation() chdir back to ' + current_dir)
    os.chdir(current_dir)

    m = 'cds_annotation() return:\nmt_file_cdsft: {mt_file_cdsft}\nannt_most_related_sp_file: {annt_most_related_sp_file}'.format(mt_file_cdsft=mt_file_cdsft, annt_most_related_sp_file=annt_most_related_sp_file)
    logger.info(m)
    return mt_file_cdsft, annt_most_related_sp_file


def tRNA_annotation(workdir=None, mt_file=None, genetic_code=None, logger=None):
    mitfipath = os.path.join(annotate_script_dir, "mitfi")
    mt_file_base = os.path.basename(mt_file)

    current_dir = os.getcwd()
    if not os.path.exists(workdir):
        os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)   
    logger.info('tRNA_annotation() chdir to ' + workdir)

    trna_out = mt_file_base + ".trna"
    trna_out = abspath(trna_out)
    command = "## begin tRNA annotation...\n" +\
            "## remove previous tRNA result file (if any)\n" +\
            "rm -rf " + trna_out
    runcmd(command, logger=logger)

    # one sequence by one sequence
    for rec in SeqIO.parse(mt_file, 'fasta'):
        tmp_file = mt_file_base + "." + rec.id
        tmp_file = os.path.abspath(tmp_file)
        fh_tmp = open(tmp_file, 'w')
        SeqIO.write(rec, fh_tmp, 'fasta')
        fh_tmp.close()

        ## MiTFi will get no results when use multiple cores!!
        ## must cat cd command to java command together here
        command = "cd " + mitfipath + "\n" +\
                  "# and find the cmsearch progrm \n" +\
                  "./detect_platform.sh \n"
        # we have to use bash instead of dash!
        # to fix the problem https://github.com/linzhi2013/MitoZ/issues/187
        # https://stackoverflow.com/questions/12230690/string-comparison-in-bash-not-found
        # and https://app.deepsource.com/directory/analyzers/docker/issues/DOK-DL4005
        
        command += "java -Xmx2048m -jar mitfi.jar -cores 1 " +\
            " -code " + genetic_code +\
            " -evalue 0.001" +\
            " -onlycutoff" +\
            " " + tmp_file + " >>" + trna_out
        runcmd(command, logger=logger)

        command = "rm -rf " + tmp_file
        runcmd(command, logger=logger)

    ## prepare tRNA feature table file
    soft = os.path.join(annotate_script_dir, "tRNA_ft_mitfi.py")
    mt_file_trnaft = trna_out + ".ft"

    command = python3 + " " + soft +\
            " " + trna_out +\
            " " + mt_file_trnaft
    runcmd(command, logger=logger)

    os.chdir(current_dir)
    logger.info('tRNA_annotation() chdir back to ' + current_dir)
    return mt_file_trnaft



def check_cmsearch_global(glresult=None):
    with open(glresult, 'r') as fh:
        fh.readline()
        fh.readline()
        i = fh.readline()
        if i.startswith("#"):
            return False
        else:
            return True


def rRNA_annotation(workdir=None, mt_file=None, thread_number=None, s_rRNA_CM=None, l_rRNA_CM=None, logger=None):
    mt_file_base = os.path.basename(mt_file)
    cmsearch_s_rRNA_result = mt_file_base + ".s-rRNA.out"
    cmsearch_s_rRNA_tbl_result = mt_file_base + ".s-rRNA.tbl"

    current_dir = os.getcwd()
    if not os.path.exists(workdir):
        os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)   
    logger.info('rRNA_annotation() chdir to ' + workdir)

    # firstly, use global searching
    command = 'cmsearch' +\
            " -g --tblout " + cmsearch_s_rRNA_tbl_result +\
            " --cpu " + thread_number +\
            " " + s_rRNA_CM +\
            " " + mt_file +\
            " >" + cmsearch_s_rRNA_result
    runcmd(command, logger=logger)

    # then, check if the result is empty.
    # If empty, use local searching,
    # otherwise, process the result
    if not check_cmsearch_global(cmsearch_s_rRNA_tbl_result):
        print("global cmsearch could not find results! Now try local cmsearch")
        command = 'cmsearch' +\
            " --tblout " + cmsearch_s_rRNA_tbl_result +\
            " --cpu " + thread_number +\
            " " + s_rRNA_CM +\
            " " + mt_file +\
            " >" + cmsearch_s_rRNA_result
        runcmd(command, logger=logger)


    s_rRNA_ft = mt_file + ".s-rRNA.ft"
    soft = os.path.join(annotate_script_dir, "rRNA_ft.py")
    command = python3 + " " + soft +\
            " " + cmsearch_s_rRNA_tbl_result +\
            " " + s_rRNA_ft
    runcmd(command, logger=logger)


    # l-rRNA
    cmsearch_l_rRNA_result = mt_file + ".l-rRNA.out"
    cmsearch_l_rRNA_tbl_result = mt_file + ".l-rRNA.tbl"

    # firstly, use global searching
    command = 'cmsearch' +\
            " -g --tblout " + cmsearch_l_rRNA_tbl_result +\
            " --cpu " + thread_number +\
            " " + l_rRNA_CM +\
            " " + mt_file +\
            " >" + cmsearch_l_rRNA_result
    runcmd(command, logger=logger)

    # then, check if the result is empty.
    # If empty, use local searching,
    # otherwise, process the result
    if not check_cmsearch_global(cmsearch_l_rRNA_tbl_result):
        print("global cmsearch could not find results! Now try local cmsearch")
        command = 'cmsearch' +\
            " --tblout " + cmsearch_l_rRNA_tbl_result +\
            " --cpu " + thread_number +\
            " " + l_rRNA_CM +\
            " " + mt_file +\
            " >" + cmsearch_l_rRNA_result
        runcmd(command, logger=logger)


    l_rRNA_ft = mt_file + ".l-rRNA.ft"
    soft = os.path.join(annotate_script_dir, "rRNA_ft.py")
    command = python3 + " " + soft +\
            " " + cmsearch_l_rRNA_tbl_result +\
            " " + l_rRNA_ft
    runcmd(command, logger=logger)

    os.chdir(current_dir)
    logger.info('rRNA_annotation() chdir back to ' + current_dir)

    m = 'rRNA_annotation() return:\ns_rRNA_ft: {s_rRNA_ft}\nl_rRNA_ft: {l_rRNA_ft}'.format(s_rRNA_ft=s_rRNA_ft, l_rRNA_ft=l_rRNA_ft)
    logger.info(m)
    return s_rRNA_ft, l_rRNA_ft


def combine_annotations_and_find_control_region(workdir=None, mt_file=None, mt_file_cdsft=None, mt_file_trnaft=None, mt_file_s_rRNA_ft=None, mt_file_l_rRNA_ft=None, defined_gene_length_f=None, s_rRNA_CM=None, l_rRNA_CM=None, genetic_code=5, species_name='test sp.', template_sbt=None, annt_most_related_sp_file=None, mitoz_module_used='annotate', logger=None):
    current_dir = os.getcwd()
    if not os.path.exists(workdir):
        os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir) 
    logger.info('combine_annotations_and_find_control_region() chdir to ' + workdir)

    mt_file_base = os.path.basename(mt_file)

    ## combine CDS and RNA annotation into a feature table file
    ## check file existing before join!!
    annot_fail_count = 0
    soft = os.path.join(annotate_script_dir, 'join_ft.py')
    tbl_file = mt_file + ".tbl"
    tbl_file = os.path.abspath(tbl_file)
    command = python3 + " " + soft +\
            " " + tbl_file

    if file_not_empty(mt_file_cdsft):
        command += " " + mt_file_cdsft
    else:
        print("can not find CDS!", file=sys.stderr)
        annot_fail_count += 1

    if file_not_empty(mt_file_trnaft):
        command += " " + mt_file_trnaft
    else:
        print("can not find tRNAs!", file=sys.stderr)
        annot_fail_count += 1

    if file_not_empty(mt_file_s_rRNA_ft):
        command += " " + mt_file_s_rRNA_ft
    else:
        print("can not find s-rRNA!", file=sys.stderr)
        annot_fail_count += 1

    if file_not_empty(mt_file_l_rRNA_ft):
        command += " " + mt_file_l_rRNA_ft
    else:
        print("can not find l-rRNA!", file=sys.stderr)
        annot_fail_count += 1

    if annot_fail_count == 4:
        sys.exit("Can not find any genes!")

    runcmd(command, logger=logger)


    # find control region
    control_region_ft = mt_file_base + ".control_region.ft"
    if file_not_empty(mt_file_cdsft) and file_not_empty(mt_file_trnaft) and file_not_empty(mt_file_s_rRNA_ft) and file_not_empty(mt_file_l_rRNA_ft):
        soft = os.path.join(annotate_script_dir, 'find_contro_region.py')
        command = python3 + " " + soft +\
                " -fa_file " + mt_file +\
                " -PCG_cutoff_file " + defined_gene_length_f +\
                " -PCG_len_ratio 0.9 " +\
                " -s_rRNA_CM_file " + s_rRNA_CM +\
                " -l_rRNA_CM_file " + l_rRNA_CM +\
                " -rRNA_len_ratio 0.9 "  +\
                " -tRNA_num_min 22 "  +\
                " -fea_files " + " ".join([mt_file_cdsft, mt_file_trnaft, mt_file_s_rRNA_ft, mt_file_l_rRNA_ft]) +\
                " -CR_len_min 600 " +\
                " -outfile " + control_region_ft
        runcmd(command, logger=logger)

        command = 'cat {0} >> {1}'.format(control_region_ft, tbl_file)
        runcmd(command, logger=logger)

    # prepare a *.fsa file
    soft = os.path.join(annotate_script_dir, 'prepare_fsa.py')
    fsa_file = mt_file + ".fsa"
    command = python3 + " " + soft +\
                " " + mt_file +\
                " '%s'" % species_name +\
                " " + genetic_code +\
                " " + fsa_file
    runcmd(command, logger=logger)

    # use the NCBI sequin soft
    #if not template_sbt:
    #    template_sbt = os.path.join(annotate_script_dir, 'template.sbt')
    command = 'tbl2asn' +\
            " -t " + template_sbt +\
            " -a s -p . -V vb"
    # runcmd2(command)
    runcmd2(command, logger=logger)

    errorsummary_val_file = os.path.abspath('errorsummary.val')

    # for the GenBank file,
    # 1. get a gene stats
    tbl2asn_gbf = mt_file_base + ".gbf"
    tbl2asn_gbf = os.path.abspath(tbl2asn_gbf)
    soft = os.path.join(annotate_script_dir, 'genbank_gene_stat_v2.py')
    command = python3 + " " + soft +\
                " " + tbl2asn_gbf +\
                " " + annt_most_related_sp_file +\
                " " + mitoz_module_used +\
                " > summary.txt"
    runcmd(command, logger=logger)

    summary_file = os.path.abspath('summary.txt')
    

    os.chdir(current_dir)
    logger.info('combine_annotations_and_find_control_region() chdir back to ' + current_dir)

    m = '''combine_annotations_and_find_control_region() return:
    tbl_file:{tbl_file}
    errorsummary_val_file: {errorsummary_val_file}
    tbl2asn_gbf: {tbl2asn_gbf}
    summary_file: {summary_file}'''.format(
        tbl_file=tbl_file,
        errorsummary_val_file=errorsummary_val_file,
        tbl2asn_gbf=tbl2asn_gbf,
        summary_file=summary_file)
    logger.info(m)

    return tbl_file, errorsummary_val_file, tbl2asn_gbf, summary_file


def extract_gene_sequences(workdir='./', gbf=None, logger=None):
    current_dir = os.getcwd()
    if not os.path.exists(workdir):
        os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir) 
    logger.info('extract_gene_sequences() chdir to ' + workdir)

    gbf_base = os.path.basename(gbf)
    command = 'mitoz-tools gbseqextractor -f {gbf} -prefix {prefix} '.format(gbf=gbf, prefix=gbf_base)
    command += '-types CDS rRNA tRNA wholeseq gene '
    command += '-cds_translation -p -l'
    runcmd(command, logger=logger)

    os.chdir(current_dir)
    logger.info('extract_gene_sequences() chdir back to ' + current_dir)


def main(args):
    logger = args.logger
    logger.info("annotate.main() got args:\n{}".format(args))

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

    if args.template_sbt:
        if Path(args.template_sbt).is_file():
            args.template_sbt = os.path.abspath(args.template_sbt)
        else:
            logger.error("Can NOT access {} !!".format(args.template_sbt))
            sys.exit(0)

    ## genetic code selection
    genetic_code_dict = {"Chordata":"2", "Arthropoda":"5", "Echinodermata":"9",
        "Annelida-segmented-worms":"5", "Bryozoa":"5", "Mollusca":"5",
        "Nematoda":"5", "Nemertea-ribbon-worms":"5", "Porifera-sponges":"4"}

    if hasattr(args, "genetic_code") and args.genetic_code == "auto":
        args.genetic_code = genetic_code_dict[args.clade]

    logger.info("annotate.main() got updated args:\n{}".format(args))

    workdir = os.path.abspath(args.workdir)
    logger.info('annotate.main() will be working in:' + workdir)

    if not os.path.exists(workdir):
        logger.info('Creating workdir:' + workdir)
        os.makedirs(workdir, exist_ok=True)

    os.chdir(workdir)
    logger.info('annotate.main() chdir:' + workdir)

    mt_files = args.fastafiles
    genetic_code = args.genetic_code
    thread_number = args.thread_number
    prefix = args.outprefix
    species_name = args.species_name
    template_sbt = os.path.abspath(args.template_sbt) # overwrite

    mitoz_module_used = args.mitoz_module_used

    MT_annotation_BGI_pro_db = os.path.join(args.profiles_dir, "MT_database", args.clade +"_CDS_protein.fa")
    defined_gene_length_f = os.path.join(args.profiles_dir, 'CDS_HMM', args.clade+"_CDS_length_list")
    logger.info('annotate.main() MT_annotation_BGI_pro_db:' + MT_annotation_BGI_pro_db)
    logger.info('annotate.main() defined_gene_length_f:' + defined_gene_length_f)

    s_rRNA_CM = os.path.join(args.profiles_dir, 'rRNA_CM', "v1.1_12snew.cm")
    l_rRNA_CM = os.path.join(args.profiles_dir, 'rRNA_CM',  "v1.1_16snew.cm")
    logger.info('annotate.main() s_rRNA_CM:' + s_rRNA_CM)
    logger.info('annotate.main() l_rRNA_CM:' + l_rRNA_CM)


    total_check_status = 0
    logger.info("Start to check if the sequence ids are too long (which will not work with BioPython) later.")
    for mt_file in mt_files:
        mt_file = os.path.abspath(mt_file)
        if not os.path.exists(mt_file):
            logger.error("Can NOT access {} !!!".format(mt_file))
            sys.exit(0)

        check_status = check_fafile_seqid_len(fafile=mt_file, logger=logger)
        total_check_status += check_status
    
    if total_check_status > 0:
        logger.critical("Some sequence ids are too long! Please rename them first!\nBut REMEMBER to keep 'topology=linear' or 'topology=circular' at the header line!!")


    resulting_gb_files = []
    for mt_file in mt_files:
        mt_file = os.path.abspath(mt_file)
        mt_file_base = os.path.basename(mt_file)
        if not os.path.exists(mt_file):
            logger.error("Error: can NOT find {} !!!".format(mt_file))

        logger.info('Starting to annotate '+ mt_file)

        result_wdir = os.path.join(workdir, prefix + '.' + mt_file_base +".result")
        m = 'result_wdir for {0}: {1}'.format(mt_file, result_wdir)
        logger.info(m)

        annotation_infile_base = prefix + "_" + mt_file_base + "_mitoscaf.fa"

        tmp_workdir = os.path.join(workdir, "tmp_"+annotation_infile_base)
        logger.info('working in '+ tmp_workdir)
        if not os.path.exists(tmp_workdir):
            os.makedirs(tmp_workdir, exist_ok=True)

        annotation_infile = os.path.join(tmp_workdir, annotation_infile_base)
        logger.info('annotation_infile: '+ annotation_infile)

        fh_out = open(annotation_infile, 'w')
        for rec in SeqIO.parse(mt_file, 'fasta'):
            topology = 'topology=linear'
            m = re.search(r'(topology=\w+)', rec.description)
            if m:
                topology = m.group(1)
                logger.info('Got topology information for {seqid}: {topology}'.format(seqid=rec.id, topology=topology))
            else:
                logger.warning('Can NOT find topology information for {seqid}. Now I assume it is {topology}!!!'.format(seqid=rec.id, topology=topology))

            if len(rec.id) > 13:
                logger.critical("Error: the sequence id {} is too long! Please rename it!".format(rec.id))
            else:
                # do not miss the circularity information!
                # print(">"+rec.id+"\n"+str(rec.seq), file=fh_out)
                print(">{seqid} {topology}\n{seq}".format(seqid=rec.id, topology=topology, seq=str(rec.seq)), file=fh_out)

        fh_out.close()

        mt_file_cdsft, annt_most_related_sp_file = cds_annotation(
            workdir=tmp_workdir,
            mt_file=annotation_infile,
            MT_annotation_BGI_pro_db=MT_annotation_BGI_pro_db,
            genetic_code=genetic_code,
            thread_number=thread_number,
            logger=logger)

        mt_file_trnaft = tRNA_annotation(
            workdir=tmp_workdir,
            mt_file=annotation_infile,
            genetic_code=genetic_code,
            logger=logger)

        mt_file_s_rRNA_ft, mt_file_l_rRNA_ft = rRNA_annotation(
            workdir=tmp_workdir,
            mt_file=annotation_infile,
            thread_number=thread_number,
            s_rRNA_CM=s_rRNA_CM,
            l_rRNA_CM=l_rRNA_CM,
            logger=logger)

        tbl_file, errorsummary_val_file, tbl2asn_gbf, summary_file = combine_annotations_and_find_control_region(
            workdir=tmp_workdir,
            mt_file=annotation_infile,
            mt_file_cdsft=mt_file_cdsft,
            mt_file_trnaft=mt_file_trnaft,
            mt_file_s_rRNA_ft=mt_file_s_rRNA_ft,
            mt_file_l_rRNA_ft=mt_file_l_rRNA_ft,
            defined_gene_length_f=defined_gene_length_f,
            s_rRNA_CM=s_rRNA_CM,
            l_rRNA_CM=l_rRNA_CM,
            genetic_code=genetic_code,
            species_name=species_name,
            template_sbt=template_sbt,
            annt_most_related_sp_file=annt_most_related_sp_file,
            mitoz_module_used=mitoz_module_used,
            logger=logger)

        extract_gene_sequences(workdir=tmp_workdir, gbf=tbl2asn_gbf, logger=logger)

        gather_result(annt_most_related_sp_file, tbl_file, errorsummary_val_file, tbl2asn_gbf+'*', summary_file, result_wdir=result_wdir, logger=logger)
        resulting_gb_files.append(os.path.abspath(tbl2asn_gbf))

        visualize_args = Visualize_Args(args)
        visualize_args.gb = tbl2asn_gbf
        visualize_args.outdir = os.path.join(tmp_workdir, 'mt_visualization')
        visualize_result_dir = visualize.main(visualize_args)
        circos_png = os.path.join(visualize_result_dir, 'circos.png')
        circos_svg = os.path.join(visualize_result_dir, 'circos.svg')
        circos_depth = os.path.join(visualize_result_dir, 'circos.depth.txt')
        if not os.path.exists(circos_depth):
            circos_depth = ''
        gather_result(circos_png, circos_svg, circos_depth, result_wdir=result_wdir, logger=logger)

    logger.info('\nFound resulting_gb_files:\n' + '\n'.join(resulting_gb_files))
    
    all_result_wdir = find_subdirs_with_suffix(indir=workdir, suffix='.result')
    logger.info('And related resulting files can be found from:\n' + '\n'.join(all_result_wdir))

    return resulting_gb_files, all_result_wdir
















