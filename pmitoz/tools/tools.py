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

COPYRIGHT 2019-2022 Guanliang Meng. ALL RIGHTS RESERVED.

"""

import argparse
import sys
import os
import re
import subprocess
import time
from glob import glob
import logging

import Bio
from Bio import SeqIO
from ete3 import NCBITaxa

# we must make 'mitoz' directory a python package, otherwise, the subpackages won't work!
# and during develpment, we also need to 'pip3 install -e .' the whole 'mitoz' package.

from p  mitoz.utility import utility
import importlib

__version__ = '1.0'

COMMANDS = ('bold_identification', 'taxonomy_ranks', 'msaconverter', 'gbseqextractor', 'gbfiletool', 'circle_check', 'group_seq_by_gene')

def main():
    '''
    The code below is from the secapr program and modified a bit.
    '''
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    parser = argparse.ArgumentParser(description=__doc__, prog='mitoz-tools', formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

    subparsers = parser.add_subparsers()
    for command_name in COMMANDS:
        module = importlib.import_module('mitoz.tools.' + command_name)
        subparser = subparsers.add_parser(command_name,
            help=module.__doc__.split('\n')[1], description=module.__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
        subparser.set_defaults(func=module.main)
        module.add_arguments(subparser)
    '''
    ArgumentParser.parse_args(args=None, namespace=None)
    args - List of strings to parse. The default is taken from sys.argv.
    '''
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    args = parser.parse_args()

    if not hasattr(args, 'func'):
        parser.error('Please provide the name of a subcommand to run')
    else:        
        args.func(args)

if __name__ == '__main__':
    main()
