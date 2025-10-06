"""
Extract any CDS or rNRA or tRNA DNA sequences of genes from Genbank file.

Seqid will be the value of '/gene=' or '/product=', if they both were not
present, the gene will not be output!

version 20201128:
    Now we can handle compounlocation (feature location with "join")!
    We can also output the translation for each CDS (retrived from '/translation=')

Please cite:
Guanliang Meng, Yiyuan Li, Chentao Yang, Shanlin Liu,
MitoZ: a toolkit for animal mitochondrial genome assembly, annotation
and visualization, Nucleic Acids Research, https://doi.org/10.1093/nar/gkz173

"""
name='gbseqextractor'
from .gbseqextractor import main, add_arguments