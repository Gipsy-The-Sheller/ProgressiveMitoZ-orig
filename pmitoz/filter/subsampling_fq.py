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
from pathlib import Path
from pmitoz.utility import utility


def runcmd(command, logger=None):
    logger.info(command)
    subprocess.run(command, shell=True)


def get_single_fq_file_lines(fastq_read_length=150, required_data_size=5, fq1=None, fq2=None, logger=None):
    required_data_size = float(required_data_size)
    fastq_read_length = int(fastq_read_length)

    # do not subsample
    if not required_data_size:
        logger.info("Not to subsample raw data!")
        return ''

    # pair-end data
    if fq1 and fq2:
        fq_records_per_file = int(required_data_size * 1e9 / (fastq_read_length * 2))
        
        single_fq_file_lines = fq_records_per_file * 4
        m = '''
        Paired-end data:
        fastq_read_length: {fastq_read_length}
        required raw data (Gbp): {required_data_size}
        fq_records_per_file: {fq_records_per_file:,}
        lines_per_file * 4: {line:,}
        '''.format(fastq_read_length=fastq_read_length,
            required_data_size=required_data_size, 
            fq_records_per_file=fq_records_per_file,
            line=single_fq_file_lines)

        logger.info(m)

        return fq_records_per_file

    # SE data
    elif (fq1 or fq2):
        fq_records_per_file = int(required_data_size * 1e9 / fastq_read_length)
        single_fq_file_lines = fq_records_per_file * 4
        m = '''
        Single-end data:
        fastq_read_length: {fastq_read_length}
        required raw data (Gbp): {required_data_size}
        fq_records_per_file: {fq_records_per_file:,}
        fq_records_per_file * 4: {line:,}
        '''.format(fastq_read_length=fastq_read_length,
            required_data_size=required_data_size, 
            fq_records_per_file=fq_records_per_file,
            line=single_fq_file_lines)

        logger.info(m)
        return fq_records_per_file


def get_fq_total_bases(workdir='./', fq1=None, fq2=None, logger=None):
    input_fq_stat_file = 'input.fq.stat.txt'
    
    command = "cd {workdir}\n".format(workdir=workdir)
    command += "seqkit stats --tabular -o {input_fq_stat_file} {fq1} {fq2} ".format(
        input_fq_stat_file=input_fq_stat_file,
        fq1=fq1,
        fq2=fq2)

    runcmd(command, logger=logger)

    sum_len = 0
    with open(os.path.join(workdir, input_fq_stat_file)) as fh:
        for i in fh:
            if 'sum_len' in i:
                continue
            i = i.strip()
            line = i.split()
            sum_len += int(line[4])

    return sum_len


def extract_certain_amount_fq_records(workdir='./', fastq_read_length=150, required_data_size=5, fq1=None, fq2=None, logger=None):
    fq_records_per_file = get_single_fq_file_lines(
        fastq_read_length=fastq_read_length,
        required_data_size=required_data_size,
        fq1=fq1,
        fq2=fq2,
        logger=logger)
   
    # no need to subsampling
    if not fq_records_per_file:
        logger.info("extract_certain_amount_fq_records(): since there is no need to subsampling, returns input {fq1} and {fq2} directly.".format(fq1=fq1, fq2=fq2))
        return os.path.abspath(fq1), os.path.abspath(fq2)

    sub_fq1 = "sampling_{}GB.R1.fq.gz".format(required_data_size)
    sub_fq2 = "sampling_{}GB.R2.fq.gz".format(required_data_size)

    sub_fq1 = os.path.join(workdir, sub_fq1)
    sub_fq2 = os.path.join(workdir, sub_fq2)

    # already sampled, return existing resulting files directly
    if Path(sub_fq1).is_file() and Path(sub_fq2).is_file():
        logger.info("extract_certain_amount_fq_records(): Found existing {sub_fq1} and {sub_fq2}, there is no need to further subsampling, returns them directly.".format(sub_fq1=sub_fq1, sub_fq2=sub_fq2))
        return sub_fq1, sub_fq2

    # perform sampling
    # pair-end data
    if fq1 and fq2:
        command = "cd {workdir}\n".format(workdir=workdir)
        command += "seqkit head -n {fq_records_per_file} -o {sub_fq1} {fq1} \n".format(
            fq_records_per_file=fq_records_per_file,
            fq1=fq1,
            sub_fq1=sub_fq1)
        command += "seqkit head -n {fq_records_per_file} -o {sub_fq2} {fq2} \n".format(
            fq_records_per_file=fq_records_per_file,
            fq2=fq2,
            sub_fq2=sub_fq2)
        runcmd(command, logger=logger)

        if Path(sub_fq1).is_file() and Path(sub_fq2).is_file():
            logger.info("extract_certain_amount_fq_records(): return subsampling results:\n{sub_fq1}\n{sub_fq2}".format(sub_fq1=sub_fq1, sub_fq2=sub_fq2))
            return sub_fq1, sub_fq2
        else:
            logger.error("extract_certain_amount_fq_records(): Can not find sub_fq1: {sub_fq1} and sub_fq2 {sub_fq2}".format(sub_fq1=sub_fq1, sub_fq2=sub_fq2))

    elif fq1:
        command = "cd {workdir}\n".format(workdir=workdir)
        command += "seqkit head -n {fq_records_per_file} -o {sub_fq1} {fq1} \n".format(
            fq_records_per_file=fq_records_per_file,
            fq1=fq1,
            sub_fq1=sub_fq1)
        runcmd(command, logger=logger)

        if Path(sub_fq1).is_file():
            logger.info("extract_certain_amount_fq_records(): return subsampling results:\n{sub_fq1}".format(sub_fq1=sub_fq1))
            return sub_fq1, ''
        else:
            logger.error("extract_certain_amount_fq_records(): Can not find sub_fq1: {sub_fq1}".format(sub_fq1=sub_fq1))

    else:
        logger.error("extract_certain_amount_fq_records(): Can not find input fq1 and/or fq2")



def main_subsampling_fq(workdir='./', fq1=None, fq2=None, outfq1=None, outfq2=None, fastq_read_length=150, thread_number=4, required_data_size=5, logger=None):
    required_data_size = float(required_data_size)

    input_fq_sum_len = get_fq_total_bases(workdir=workdir, fq1=fq1, fq2=fq2, logger=logger)

    if input_fq_sum_len < (required_data_size * 1e9):
        m = """
        Input fastq size: {input_fq_sum_len:,} bp
        required_data_size: {required_data_size} Gb
        """.format(input_fq_sum_len="{num:,}".format(num=input_fq_sum_len),
            required_data_size=required_data_size)
        logger.warn(m)

        return fq1, fq2

    sub_fq1, sub_fq2 = extract_certain_amount_fq_records(
        workdir=workdir,
        fastq_read_length=fastq_read_length,
        required_data_size=required_data_size,
        fq1=fq1,
        fq2=fq2,
        logger=logger)

    return sub_fq1, sub_fq2 


if __name__ == '__main__':
    main_subsampling_fq(sys.argv)