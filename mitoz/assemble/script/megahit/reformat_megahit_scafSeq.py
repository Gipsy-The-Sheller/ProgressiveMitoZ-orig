#!/usr/bin/env python3
import os
import sys

usage = '''
python3 {} <raw_scafSeq_file> <outfile>
'''.format(sys.argv[0])

if len(sys.argv) != 3:
    sys.exit(usage)


raw_scafSeq_file, reformatted_scafSeq_file = sys.argv[1:]

#>k141_1709277 flag=0 multi=1.0000 len=507

with open(raw_scafSeq_file, 'r') as fh, open(reformatted_scafSeq_file, 'w') as fhout:
    for i in fh:
        i = i.rstrip()
        if i.startswith('>'):
            line = i.split()
            seqid = line[0]
            abun = 'NA'

            if i.startswith('>k'):
                abun = line[2].replace('multi=', '')

            desc = "{seqid} abun={abun}".format(seqid=seqid, abun=abun)
            print(desc, file=fhout)
        else:
            print(i, file=fhout)