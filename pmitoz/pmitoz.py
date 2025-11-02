#!/usr/bin/env python3
"""
MitoZ 3.6

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

When you use MitoZ, please cite:
Meng, G., Li, Y., Yang, C., & Liu, S. (2019). MitoZ: a toolkit for
animal mitochondrial genome assembly, annotation and visualization.
Nucleic acids research, 47(11), e63-e63. https://doi.org/10.1093/nar/gkz173

Additionaly, you SHOULD also CITE related software invoked by MitoZ,
without them, MitoZ will not make its way. For more details, please
refer to the MitoZ paper and https://github.com/linzhi2013/MitoZ.

"""

__version__ = '3.6'

import argparse
import sys
import os
import re
import subprocess
import time
from glob import glob
import logging
import shutil

import Bio
from Bio import SeqIO
from ete3 import NCBITaxa

# we must make 'mitoz' directory a python package, otherwise, the subpackages won't work!
# and during develpment, we also need to 'pip3 install -e .' the whole 'mitoz' package.

from pmitoz.utility import utility
import importlib


COMMANDS = ('filter', 'assemble', 'findmitoscaf', 'annotate', 'visualize', 'all')

def main():
    '''
    The code below is from the secapr program and modified a bit.
    '''
    # logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    parser = argparse.ArgumentParser(description=__doc__, prog='mitoz', formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument('--debug', dest='debug', action='store_true', default=False,
        help='debug mode output [%(default)s]')
    
    subparsers = parser.add_subparsers()
    for command_name in COMMANDS:
        module = importlib.import_module('pmitoz.' + command_name)
        subparser = subparsers.add_parser(command_name,
            help=module.__doc__.split('\n')[1], description=module.__doc__)
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
    args.logger = utility.get_logger(debug=args.debug)
    # logger.info()
    # logger.debug(e)

    if not hasattr(args, 'func'):
        parser.error('Please provide the name of a subcommand to run')
    else:
        # check if main executables all exist!
        required_executables = ['spades.py', 'megahit', ]
        for exe in required_executables:
            exe_full_path = shutil.which(exe)
            if exe_full_path:
                args.logger.info("Found: " + exe_full_path)
            else:
                args.logger.warn(exe + " NOT found ") 

        # all executable found, we are ready to go
        args.func(args)

if __name__ == '__main__':
    main()
