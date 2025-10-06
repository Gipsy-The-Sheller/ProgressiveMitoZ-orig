#!/usr/bin/env python3

import os
import re
import sys
import time
import subprocess
from pathlib import Path

import logging

def runcmd(command, logger=None):
    try:
        if logger:
            logger.info(command)
        subprocess.check_call(command, shell=True)
    except Exception as e:
        m = "Error occured when running command:\n{0}\n".format(command, e)
        if logger:
            logger.error(m)
        else:
            sys.exit(m)

def runcmd2(command, logger=None):
    try:
        if logger:
            logger.info(command)
        subprocess.call(command, shell=True)
    except Exception as e:
        m = "Error occured when running command:\n{0}\n".format(command, e)
        if logger:
            logger.error(m)
        else:
            sys.exit(m)

def pre_del_cmd(prefix=None, filestr=None):
    file_list = [prefix + '*' + j for j in filestr.split()]
    file_list = " ".join(file_list)
    command = "rm -rf " + file_list
    return command

def gather_result(*files, logger=None, result_wdir=None):
    if not os.path.isdir(Path(result_wdir)):
        logger.info("gather_result() creating: " + result_wdir)
        os.mkdir(result_wdir)
    flist = " ".join(files)
    command = "cp -r " + flist + " " + result_wdir
    runcmd(command, logger=logger)

def file_not_empty(file=None):
    """
    check if file is empty

    """
    if os.stat(file).st_size > 0:
        return True
    else:
        return False

def check_program_invoked(cmd):
    result = subprocess.call('type %s' % cmd, shell = True,
            stdout = subprocess.PIPE, stderr = subprocess.PIPE) == 0
    if result:
        return 0
    else:
        print(cmd, " not found!", file=sys.stderr)
        return 1


def files_exist_0_or_1(filelist):
    num = 0
    for file in filelist:
        if os.path.exists(file):
            num += 1
        else:
            print("%s doesn't exist!" % file, file=sys.stderr)
    if len(filelist) == num:
        return 0
    else:
        return 1


def directory_exist_check(*directories):
    err = 0
    for directory in directories:
        if os.path.exists(directory):
            print("Directory %s exists, please delete it!" % directory,
                file=sys.stderr)
            err += 1
    if err > 0:
        sys.exit('Exit!')


def abspath(*files):
    flist = []
    for file in files:
        file = os.path.abspath(file)
        flist.append(file)
    if len(flist) > 1:
        return flist[:]
    else:
        return flist[0]


def find_subdirs_with_suffix(indir='./', suffix='.result'):
    all_subdirs = []
    for p in Path(indir).iterdir():
        if str(p).endswith('.result') and p.is_dir():
            all_subdirs.append(str(p))

    return all_subdirs



def get_logger(debug=False):
    # 级别排序:CRITICAL > ERROR > WARNING > INFO > DEBUG
    formatter = logging.Formatter(
        "\n%(asctime)s - %(name)s - %(levelname)s - \n%(message)s")

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)  # must be DEBUG, then 'ch' below works.
    # logger.setFormatter(formatter)

    fh = logging.FileHandler(os.path.basename(sys.argv[0]) + '.log')
    if debug:
        fh.setLevel(logging.DEBUG)
    else:
        fh.setLevel(logging.INFO)  # INFO level goes to the log file
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    ch = logging.StreamHandler()
    if debug:
        ch.setLevel(logging.DEBUG)
    else:
        ch.setLevel(logging.INFO)  # only INFO level will output on screen
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    return logger


def check_fafile_seqid_len(fafile=None, logger=None):
    error = 0
    with open(fafile, 'r') as fh:
        for i in fh:
            i = i.strip()
            if not i or not i.startswith('>'):
                continue
            seqid = i.split()[0]

            m = re.search(r'(topology=\w+)', i)

            if len(seqid) > 13:
                error += 1
                logger.warning("Error: file {} the sequence id '{}' is too long! Please rename it! But REMEMBER to keep 'topology=linear' or 'topology=circular' at the header line!!".format(fafile, seqid))
            else:
                if m:
                    logger.warning("File {}, the sequence id '{}': length is okay! And there is topology information '{}' at the header line!!".format(fafile, seqid, m.group(1)))
                else:
                    logger.warning("File {}, the sequence id {}: length is okay! But NO topology information ('topology=linear' or 'topology=circular') found at the header line!!".format(fafile, seqid))

    return error


