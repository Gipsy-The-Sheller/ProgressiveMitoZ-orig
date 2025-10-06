'''
To get the lineage information of input taxid, species name, or higher ranks
(e.g., Family, Order) with ETE3 package.

The ete3 package will automatically download the NCBI Taxonomy database during
 the first time using of this program.

Please be informed:

(1) A rank name may have more than one taxids, e.g., Pieris can means:
Pieris <butterfly> and Pieris <angiosperm>. I will search the lineages for
both of them.

(2) When you give a species name, if I can not find the taxid for the species
name, I will try to search the first word (Genus).

(3) Any input without result found will be output in outfile.err ('-o' option).
'''
name='taxonomy_ranks'
from .get_taxonomy_rank_with_ete3_with_super_and_sub import TaxonomyRanks, main, add_arguments
