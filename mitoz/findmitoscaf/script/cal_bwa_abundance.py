#!/usr/bin/python3
"""
cal_bwa_abundance.py

Copyright (c) 2017-2018 Guanliang Meng <mengguanliang@foxmail.com>.

This file is part of MitoZ.

MitoZ is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MitoZ is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MitoZ.  If not, see <http://www.gnu.org/licenses/>.

"""

import argparse
import sys
import os
import subprocess
from Bio import SeqIO
import time

description = "do BWA alignment and calculate average sequencing depth."


parser = argparse.ArgumentParser(description=description)

parser.add_argument("-fa", metavar="<STR>", required=True,
				help="input fasta file")

parser.add_argument("-fq1", metavar="<STR>", required=True,
				help="input fastq1 file")

parser.add_argument("-fq2", metavar="<STR>",
				help="input fastq2 file")

parser.add_argument("-out", metavar="<STR>", required=True,
				help="output file")

parser.add_argument("-bwa", metavar="<STR>", default="bwa",
				help="bwa command [%(default)s]")

parser.add_argument("-thread", metavar="<INT>", default="1",
				help="bwa thread number [%(default)s]")

parser.add_argument("-samtools", metavar="<STR>", default="samtools",
				help="samtools command [%(default)s]")

if len(sys.argv) == 1:
	parser.print_help()
	parser.exit()
else:
	args = parser.parse_args()

python3 = sys.executable

def reads_mapping(mitoscaf_file=None, fq1=None, fq2=None, bwa=None, samtools=None):
	# if file are not in current directory
	dirname = os.path.dirname(os.path.abspath(mitoscaf_file))
	basename = os.path.basename(mitoscaf_file)
	current_dir = os.getcwd()
	if dirname != current_dir:
		command = "ln -s " + mitoscaf_file + " "  + basename
		runcmd(command)
		mitoscaf_file = basename

	## bwa indexing
	command = bwa + " index " + mitoscaf_file
	runcmd(command)

	bwa_bam = basename + ".bwa.bam"
	command = bwa + " mem" +\
			" -a -t " + args.thread +\
			" " +  mitoscaf_file +\
			" " + fq1 +\
			" " + fq2 +\
			" | " + samtools + " view -b -h -F 4 - " +\
			" | samtools sort --threads " + args.thread +\
			" --output-fmt BAM -o " + bwa_bam
	runcmd(command)

	# must be sorted first, otherwise got error:
	# samtools depth: Data is not position sorted
	depth_file = basename + ".depth"
	command = "samtools depth -aa {bam} > {depth_file}".format(bam=bwa_bam, depth_file=depth_file)
	# The output is a simple tab-separated table with three columns: reference name, position, and coverage depth.
	runcmd(command)

	files_to_delete = [basename+j for j in (".amb", ".ann", ".bwt", ".pac", ".rbwt", ".rsa", ".sa", ".rpac")]
	files_to_delete = " ".join(files_to_delete)
	command = "rm -rf " + files_to_delete
	runcmd(command)

	return depth_file

def cal_avg_depth(fasta=None, depth_file=None, outfile=None):
	total_depth_dict = {}
	seq_len_dict = {}
	with open(depth_file, 'r') as fh_in:
		for i in fh_in:
			i = i.strip()
			seqid, pos, depth = i.split("\t")
			if seqid not in total_depth_dict:
				total_depth_dict[seqid] = 0
			total_depth_dict[seqid] += float(depth)

			seq_len_dict[seqid] = int(pos) # the last pos is the seq length

	avg_depth_dict = {}
	for seqid in total_depth_dict:
		tot_len = total_depth_dict[seqid]
		seq_len = seq_len_dict[seqid]

		avg_depth = tot_len / seq_len
		avg_depth_dict[seqid] = str(round(avg_depth, 2))

	fh_out = open(outfile, 'w')
	for rec in SeqIO.parse(fasta, 'fasta'):
		if rec.id in avg_depth_dict:
			print(">"+rec.id, "abun="+str(avg_depth_dict[rec.id]), "length="+str(len(rec)), file=fh_out)
			print(rec.seq, file=fh_out)
		else:
			print("can not find average depth of " + rec.id, file=sys.stderr)
			print("So we set it to be 0X!!", file=sys.stderr)
			print(">"+rec.id, "abun=0", "length="+str(len(rec)), file=fh_out)
			print(rec.seq, file=fh_out)
	fh_out.close()


def runcmd(command):
	try:
		current_time = time.strftime("%Y-%m-%d %H:%M:%S",
						time.localtime(time.time()))
		print(current_time, "\n", command, "\n", sep="", flush=True)
		subprocess.check_call(command, shell=True)
	except:
		sys.exit("Error occured when running command:\n%s" % command)


depth_file = reads_mapping(mitoscaf_file=args.fa, fq1=args.fq1,
				fq2=args.fq2, bwa=args.bwa, samtools=args.samtools)

cal_avg_depth(fasta=args.fa, depth_file=depth_file, outfile=args.out)
