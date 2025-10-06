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

COPYRIGHT © 2019-2022 Guanliang Meng. ALL RIGHTS RESERVED.

"""

import sys
import os
import argparse
import time
import subprocess
from string import Template # Python >=3.6
from pathlib import Path
import glob
import json


def add_arguments(parser):
    parser.add_argument('--fq1', metavar='<file>', required=True,
        help="Fastq1 file")

    parser.add_argument('--fq2', metavar='<file>',
        help="Fastq2 file")

    parser.add_argument('--phred64', action='store_true',
        help="Are the fastq phred64 encoded? [%(default)s]")

    parser.add_argument('--outprefix', metavar='<str>', default='out',
        help="output prefix [%(default)s]")

    parser.add_argument("--fastq_read_length", metavar="<INT>",  type=int, default=150,
        help='''read length of fastq reads, used to split clean fastq files. [%(default)s]''')

    parser.add_argument("--data_size_for_mt_assembly", metavar="<float1>,<float2>", default="5,0",
        help='''Data size (Gbp) used for mitochondrial genome assembly, usually between 3~8 Gbp is enough.
        The float1 means the size (Gbp) of raw data to be subsampled, while the float2 means the size of 
        clean data should be >= float2 Gbp, otherwise MitoZ will stop to run. When only float1 is set,
        float2 is assumed to be 0.  Set float1 to be 0 if you want to use ALL raw data. [%(default)s]''')

    parser.add_argument('--filter_other_para', metavar='<str>', default='',
        help="other parar. [%(default)s]")

    parser.add_argument('--thread_number', metavar='<int>', type=int, default=4,
        help="thread number [%(default)s]")

    parser.add_argument('--workdir', metavar='<directory>', default='./',
        help="working directory [%(default)s]")

    parser.add_argument('--workdir_done', metavar='<directory>', default='./done',
        help="done directory [%(default)s]")

    parser.add_argument('--workdir_log', metavar='<directory>', default='./log',
        help="log directory [%(default)s]")

    return parser


def filter_rawdata(fq1=None, fq2=None, phred64=False, outprefix='out', reads_to_process=0, other_para='', thread=4, workdir='./', workdir_done='./done/', workdir_log='./log/', logger=None):
    out1 = outprefix + ".clean_R1.fq.gz"
    out2 = outprefix + ".clean_R2.fq.gz"
    report_json = outprefix + ".QC.json"
    report_html = outprefix + ".QC.html"

    filtered_fq1 = os.path.join(workdir, out1)
    filtered_fq2 = os.path.join(workdir, out2)

    if phred64:
        phred64_value = '--phred64'
    else:
        phred64_value = ''

    command_template = Template('''
    # set -vex
    # set -o pipefail
    set -e

    mkdir -p ${workdir_done}
    mkdir -p ${workdir_log}
    mkdir -p ${workdir}
    cd ${workdir}

    task_id="FastpTask"
    done_file="${workdir_done}/$task_id.done"
    log_o_file="${workdir_log}/$task_id.o.log"
    log_e_file="${workdir_log}/$task_id.e.log"

    if [ -f $done_file ]; then
        exit 0
    fi

    if [ -f "${fq1}" ] && [ -f "${fq2}" ];
    then
        fastp ${other_para} \
        -i ${fq1} \
        -I ${fq2} ${phred64} \
        --length_required 50 \
        --reads_to_process ${reads_to_process} \
        -o ${out1} \
        -O ${out2} \
        -R "${outprefix}" \
        -j ${report_json} \
        -h ${report_html} \
        -w ${thread} && 1>$log_o_file 2>$log_e_file && touch $done_file   
 
    else
        fastp ${other_para} \
        -i ${fq1} ${phred64} \
        --length_required 50 \
        --reads_to_process ${reads_to_process} \
        -o ${out1} \
        -R "${outprefix}" \
        -j ${report_json} \
        -h ${report_html} \
        -w ${thread} && 1>$log_o_file 2>$log_e_file && touch $done_file
    
    fi
    ''')

    command = command_template.safe_substitute(
        workdir=workdir,
        workdir_done=workdir_done,
        workdir_log=workdir_log,
        fq1=fq1,
        fq2=fq2,
        out1=out1,
        out2=out2,
        other_para=other_para,
        outprefix=outprefix,
        reads_to_process=reads_to_process,
        report_json=report_json,
        report_html=report_html,
        thread=thread,
        phred64=phred64_value,
        )

    # logger.info("Starting to filter raw data {}".format(command))
    runcmd(command, logger=logger)

    if fq1 and not fq2:
        return os.path.abspath(filtered_fq1), '', os.path.join(workdir, report_json)
    else:
        return os.path.abspath(filtered_fq1), os.path.abspath(filtered_fq2), os.path.join(workdir, report_json)


def runcmd(command, logger=None):
    logger.info(command)
    subprocess.run(command, shell=True)

    '''
    try:
        current_time = time.strftime("%Y-%m-%d %H:%M:%S",
                        time.localtime(time.time()))
        print(current_time, "\n", command, "\n", sep="", flush=True)
        subprocess.check_call(command, shell=True)
    except:
        sys.exit("Error occured when running command:\n%s" % command)
    '''

def get_single_fq_file_lines(fastq_read_length=150, data_size=5, fq1=None, fq2=None, logger=None):
    # do not split
    if not data_size:
        logger.info("Not to subsample raw data!")
        return 0

    # pair-end data
    if fq1 and fq2:
        fq_records_per_file = int(data_size * 1e9 / (fastq_read_length * 2))
        
        single_fq_file_lines = fq_records_per_file * 4
        m = '''
        Paired-end data:
        fastq_read_length: {fastq_read_length}
        required raw data (Gbp): {data_size}
        fq_records_per_file: {fq_records_per_file:,}
        lines_per_file * 4: {line:,}
        '''.format(fastq_read_length=fastq_read_length,
            data_size=data_size, 
            fq_records_per_file=fq_records_per_file,
            line=single_fq_file_lines)

        logger.info(m)

        return fq_records_per_file

    # SE data
    elif (fq1 or fq2):
        fq_records_per_file = int(data_size * 1e9 / fastq_read_length)
        single_fq_file_lines = fq_records_per_file * 4
        m = '''
        Single-end data:
        fastq_read_length: {fastq_read_length}
        required raw data (Gbp): {data_size}
        fq_records_per_file: {fq_records_per_file:,}
        fq_records_per_file * 4: {line:,}
        '''.format(fastq_read_length=fastq_read_length,
            data_size=data_size, 
            fq_records_per_file=fq_records_per_file,
            line=single_fq_file_lines)

        logger.info(m)
        return fq_records_per_file


def get_clean_data_size(fastp_json_file=None, logger=None):
    with open(fastp_json_file, 'r') as fh:
        data = json.load(fh)

    message_template = Template('''
    sequencing: ${sequencing}

    before_filtering:
        total_reads: ${raw_total_reads}
        total_bases: ${raw_total_bases}

    after_filtering:
        total_reads: ${clean_total_reads}
        total_bases: ${clean_total_bases}

    ''')

    message = message_template.safe_substitute(
        sequencing=data['summary']['sequencing'],
        raw_total_reads="{num:,}".format(num=data['summary']['before_filtering']['total_reads']),
        raw_total_bases="{num:,}".format(num=data['summary']['before_filtering']['total_bases']),
        clean_total_reads="{num:,}".format(num=data['summary']['after_filtering']['total_reads']),
        clean_total_bases="{num:,}".format(num=data['summary']['after_filtering']['total_bases']),
        )

    logger.info(message)

    return data['summary']['after_filtering']['total_bases']



def main(args):
    logger = args.logger
    logger.info("filter.main() get args: {}".format(args))

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

    args.workdir = os.path.abspath(args.workdir)

    logger.info("filter.main() get updated args: {}".format(args))

    line = [float(j) for j in args.data_size_for_mt_assembly.split(',')]
    raw_data_size = required_clean_data_size = 0
    if len(line) >= 2:
        raw_data_size = float(line[0])
        required_clean_data_size = float(line[1])
    elif len(line) == 1:
        raw_data_size = float(line[0])

    reads_to_process = get_single_fq_file_lines(
        fastq_read_length=args.fastq_read_length,
        data_size=raw_data_size,
        fq1=args.fq1,
        fq2=args.fq2,
        logger=logger)

    filtered_fq1, filtered_fq2, report_json = filter_rawdata(
        fq1=args.fq1,
        fq2=args.fq2,
        phred64=args.phred64,
        outprefix=args.outprefix,
        reads_to_process=reads_to_process,
        other_para=args.filter_other_para,
        thread=args.thread_number,
        workdir=args.workdir,
        workdir_done=args.workdir_done,
        workdir_log=args.workdir_log,
        logger=logger)

    actual_clean_data_size = get_clean_data_size(fastp_json_file=report_json, logger=logger)
    if required_clean_data_size > 0:
        if actual_clean_data_size < (required_clean_data_size * 1e9):
            m = """
The clean data has size of {actual_clean_data_size:,} bp  < required {required_clean_data_size} Gb.

Exit now!

See the --data_size_for_mt_assembly <float1>,<float2> option!

If you want to continue MitoZ with such volume of clean data, set 

--data_size_for_mt_assembly {raw_data_size},0 
OR,
increase the size of <float1>.
            """.format(raw_data_size=raw_data_size,
                actual_clean_data_size=actual_clean_data_size,
                required_clean_data_size=required_clean_data_size)
            logger.error(m)
            sys.exit(0)
    '''
    sub_fq1, sub_fq2 = extract_certain_amount_fq_records(
        workdir=args.workdir,
        fastq_read_length=args.fastq_read_length,
        data_size_for_mt_assembly=args.data_size_for_mt_assembly,
        fq1=filtered_fq1,
        fq2=filtered_fq2,
        logger=logger)

    '''
    # logger.info('filter.main() returns:\n{0}\n{1}'.format(sub_fq1, sub_fq2))
    # return sub_fq1, sub_fq2

    logger.info('filter.main() returns:\n{0}\n{1}'.format(filtered_fq1, filtered_fq2))
    
    return filtered_fq1, filtered_fq2

