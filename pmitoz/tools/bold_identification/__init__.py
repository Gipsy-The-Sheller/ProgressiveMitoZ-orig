'''
To identify taxa of given sequences from BOLD (http://www.boldsystems.org/).

Some sequences can fail to get taxon information, which can be caused by
TimeoutException if your network to the BOLD server is bad.
Those sequences will be output in the file '*.TimeoutException.fasta'.

You can:
1) run another searching with the same command directly (but add -c option);
2) lengthen the time to wait for each query (-t option);
3) increase submission times (-r option) for a sequence.

Also, the sequences without BOLD matches will be output in the
file '*.NoBoldMatchError.fasta'

By Guanliang Meng.
See https://github.com/linzhi2013/bold_identification.

Citation:
Yang C, Zheng Y, Tan S, Meng G, et al.
Efficient COI barcoding using high throughput single-end 400 bp sequencing.
https://doi.org/10.1186/s12864-020-07255-w
'''

name='bold_identification'
from .BOLD_identification import main, add_arguments