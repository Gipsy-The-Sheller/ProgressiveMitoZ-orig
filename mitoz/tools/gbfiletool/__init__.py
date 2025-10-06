"""
A tool to deal with genbank records: cut, comrev, sort and select

Version
    0.0.1

Author
    mengguanliang@genomics.cn, BGI.

Please cite:
Guanliang Meng, Yiyuan Li, Chentao Yang, Shanlin Liu,
MitoZ: a toolkit for animal mitochondrial genome assembly, annotation
and visualization, Nucleic Acids Research, https://doi.org/10.1093/nar/gkz173

"""
name='gbfiletool'
from .genbank_file_tool import main, add_arguments