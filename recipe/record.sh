# conda install mamba -n base -c conda-forge
# see https://github.com/mamba-org/mamba
# and then do:
# ~/.conda/envs/mybase/bin/mamba  env create -n mitozEnv2  -f mitozEnv.yaml
#
git clone git@github.com:linzhi2013/bioconda-recipes.git

cd bioconda-recipes
git remote add upstream git@github.com:bioconda/bioconda-recipes.git

# now prepare the meta.yaml in the 'recipes/mitoz/' directory
# But WARNING: in the 'build' and 'run' section, we must use the same version, e.g. '- python 3.9'
# Otherwise, the bioconda-utils command below always fail!!!!
# 
bioconda-utils build recipes config.yml --packages  mitoz

# see https://anaconda.org/guanliangmeng/dashboard
# if upload:

conda deactivate
# out of any 'base' env
/home/gmeng/.conda/envs/mybase/bin/mamba env create -n mitozEnv.test3.6 -f ../recipe/mitozEnv.yml
