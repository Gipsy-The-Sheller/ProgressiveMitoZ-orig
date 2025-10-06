#!/usr/bin/env python3
import os
import sys

usage = '''
python3 {} <raw_scafSeq_file> <outfile>
'''.format(sys.argv[0])

if len(sys.argv) != 3:
    sys.exit(usage)


raw_scafSeq_file, reformatted_scafSeq_file = sys.argv[1:]

# >NODE_43646_length_2959_cov_2.862603 abun=NA

with open(raw_scafSeq_file, 'r') as fh, open(reformatted_scafSeq_file, 'w') as fhout:
    for i in fh:
        i = i.rstrip()
        if i.startswith('>'):
            line = i.split('_')
            seqid = line[0] + '_' + line[1]
            abun = 'NA'

            if i.startswith('>NODE_'):
                abun = line[5]

            desc = "{seqid} abun={abun}".format(seqid=seqid, abun=abun)
            print(desc, file=fhout)
        else:
            print(i, file=fhout)