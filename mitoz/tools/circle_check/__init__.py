"""
Check whether the sequences are circular when the sequences have length >= 12Kbp

output files:

1. <outPrefix>.mitogenome.fa
All the sequences from <in.fasta>.

The sequence id line will be like:
>C1 topology=circular
>C2 topology=linear

For the circular mt sequence, the overlapping region (the second `ATGCNN`
below) has been removed (below is an example)

ATGCNNNNN[ATGCNN]

Assuming `ATGCNNNNN` is a circular mt sequence, `ATGCNN` are the overlapping
regions.

2. <outPrefix>.start2end_for-circular-mt-only
This file contains the circular sequences only, and the first 300 bp of each
has been moved to the end of the sequence, just for better reads mapping. You
can check the sequencing depth around the 'joining site' (-300 bp) using the
`annotate` module of MitoZ, to confirm if the sequence is really circular.

3. <outPrefix>.overlap_information
The overlapping sequence detected for the circular sequences.


Author
    mengguanliang@genomics.cn, BGI.

Please cite:
Guanliang Meng, Yiyuan Li, Chentao Yang, Shanlin Liu,
MitoZ: a toolkit for animal mitochondrial genome assembly, annotation
and visualization, Nucleic Acids Research, https://doi.org/10.1093/nar/gkz173

"""
name='circle_check'
from .circle_check import main, add_arguments