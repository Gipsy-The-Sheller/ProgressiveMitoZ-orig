The easiest way to install MitoZ is to use Singularity and Docker MitoZ images, which can help you avoid a lot of potential installation troubles.

# 1. Docker

## 1.1 Install Docker

Any platform (e.g. Linux, Mac or Windows) on which Docker is able to run should be able to run MitoZ via the MitoZ Docker image. This also applies to [Singularity](#2-singularity).

Please refer to https://docs.docker.com/.

## 1.2 Download the MitoZ container

    $ docker pull guanliangmeng/mitoz:3.4

- with the Docker image, you don't need to install the **etetoolkit (NCBI Taxonomy) database** by yourself, everything has been packaged into the Docker image.

## 1.3 Run the container

In your working directory (i.e. `$PWD`) (the fastq files should be in there),
shell into the container:

    $ docker run -v $PWD:$PWD -w $PWD --rm guanliangmeng/mitoz:3.4 mitoz -h
    $ docker run -v $PWD:$PWD -w $PWD --rm guanliangmeng/mitoz:3.4 mitoz-tools -h

The `-v $PWD:$PWD` here means mounting your current host directory into `$PWD` of the Docker container. Multiple `-v` options can be used at the same time. **Only in this way, you can access the files under the `$PWD` within the docker container.** But within the docker container, we won't be able to access any other files (or soft-links or maybe hard-links) outside the `$PWD` directory.

## 1.4 You can also `shell` into  the container

In your **host** working directory (i.e. `$PWD`) (**the fastq files should be in there**),
shell into the container:

    $ docker run -v $PWD:$PWD -w $PWD --rm -it guanliangmeng/mitoz:3.4 

In the container,

    $ mitoz

    $ mitoz-tools


To learn more about Docker usage, please go to https://docs.docker.com/.

## 1.5 Installation location within the Docker image

With the Docker image, MitoZ is installed `/app/anaconda/bin/mitoz` and `/app/anaconda/lib/python3.9/site-packages/mitoz`:
```bash
$ docker run -it -v  $PWD:$PWD -w $PWD --rm guanliangmeng/mitoz:3.4
root@cb99de738f74:/Users/gmeng# ls -lhrt /app/anaconda/lib/python3.9/site-packages/mitoz
total 48K
-rw-rw-r-- 2 root root   12 Jun 10 08:54 __init__.py
-rw-rw-r-- 2 root root 3.2K Jun 10 08:54 MitoZ.py
drwxr-xr-x 4 root root 4.0K Jul  1 13:47 annotate
drwxr-xr-x 3 root root 4.0K Jul  1 13:47 utility
drwxr-xr-x 4 root root 4.0K Jul  1 13:47 assemble
drwxr-xr-x 3 root root 4.0K Jul  1 13:47 all
drwxr-xr-x 3 root root 4.0K Jul  1 13:47 visualize
drwxr-xr-x 7 root root 4.0K Jul  1 13:47 tools
drwxr-xr-x 6 root root 4.0K Jul  1 13:47 profiles
drwxr-xr-x 4 root root 4.0K Jul  1 13:47 findmitoscaf
drwxr-xr-x 3 root root 4.0K Jul  1 13:47 filter
drwxr-xr-x 2 root root 4.0K Jul  1 13:47 __pycache__
```

And MitoZ's database is at:
```bash
root@cb99de738f74:/Users/gmeng# ls -lhrt /app/anaconda/lib/python3.9/site-packages/mitoz/profiles/
total 16K
-rw-rw-r-- 2 root root    0 Jun 10 08:54 __init__.py
drwxr-xr-x 2 root root 4.0K Jul  1 13:47 rRNA_CM
drwxr-xr-x 2 root root 4.0K Jul  1 13:47 MT_database
drwxr-xr-x 2 root root 4.0K Jul  1 13:47 CDS_HMM
drwxr-xr-x 2 root root 4.0K Jul  1 13:47 __pycache__
```

If you want to copy this database out of the Docker image, do:
```bash
$ cd ~
$ mkdir mitoz_custom_db
$ docker run -v $PWD:$PWD -w $PWD --rm -it guanliangmeng/mitoz:3.4 
root@cb99de738f74:/Users/gmeng# cp -a /app/anaconda/lib/python3.9/site-packages/mitoz/profiles mitoz_custom_db
$ exit

# this way, the 'profiles' directory is copied to the  'mitoz_custom_db' of your host machine.

# Later, if you want to use the '--profiles_dir' option, you need to use the Docker's '-v' option
# to map this host's 'mitoz_custom_db' directory into the Docker container via

$ docker run -v $PWD:$PWD -v ~/mitoz_custom_db:/mitoz_custom_db/ -w $PWD --rm -it guanliangmeng/mitoz:3.4
# then within the Docker container:
root@cb99de738f74:/Users/gmeng# mitoz --profiles_dir /mitoz_custom_db/profiles <other options>

```

See also https://github.com/linzhi2013/MitoZ/wiki/Extending-MitoZ%27s-database.

# 2. Singularity

## 2.1 Install Singularity
See [https://www.sylabs.io/docs/](https://www.sylabs.io/docs/) for instructions to install Singularity.

Any platform (e.g. Linux, Mac or Windows) on which Singularity is able to run should be able to run MitoZ via the MitoZ Singularity image. This also applies to [Docker](#3-docker).

For the installation of Singularity on Mac or Windows, please refer to https://docs.sylabs.io/guides/3.2/user-guide/installation.html#install-on-windows-or-mac.

Note: according to the official documentation (Oct. 2019), the Singularity must be installed with root privilege.

And the Singularity installed via conda (e.g. `conda install -c bioconda singularity`) may not work (at least when installing as normal users)!

_How about Singularity on Mac and Windows?_

MitoZ only runs on Linux systems, although some of its functions can now run on Mac or Windows.

*Why do we want to run MitoZ on Mac and Windows?* There are two main reasons: 
(1) With the two new [*de novo* assemblers](https://github.com/linzhi2013/MitoZ/wiki/New-Features#2-two-new-assemblers) and small datasets, it is now possible to perform mitogenome assembly on a Mac or Windows with 16GB or 32GB RAM theoretically; 
(2) and actually only the `mitoz all` and `mitoz assemble` commands need much memory, all the other commands (`mitoz filter/findmitoscaf/annotate/visualize or `mitoz-tools`) need very little memory and thus can run on normal Mac or Windows (e.g. with 8GB RAM), and sometimes for these analyses, you do not want to upload the data to a Linux server.


## 2.2 Download the MitoZ image

You can download a pre-built Singularity (https://sylabs.io/) image from https://www.dropbox.com/sh/mqjqn656x41q2wb/AAD02t_kUCjNHbBgCeYpEM88a?dl=0 (**only for the latest version** due to limited Dropbox space); https://pan.baidu.com/s/1YIULJ9H3BeWKcIZdMZpcuw?pwd=7r9d (提取码:7r9d) 
 (**MitoZ version 3.2 and newer versions**); or from [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6581914.svg)](https://doi.org/10.5281/zenodo.6581914) (**MitoZ version 3.2 only**).

- After downloading MitoZ, you still need to install the **etetoolkit (NCBI Taxonomy) database**, especially when the automatic installation does not work for you. See [**6. The Etetoolkit database**](#6-the-etetoolkit-database) section below.


Within the Singularity image, 
MitoZ is installed at `/app/anaconda/bin/mitoz` and `/app/anaconda/lib/python3.9/site-packages/mitoz`. MitoZ's annotation database is at `/app/anaconda/lib/python3.9/site-packages/mitoz/profiles`.

## 2.3 Usage 1
```bash
$ /path/to/MitoZ_v3.4.sif -h
# for example, to use the `all` subcommand:
$ /path/to/MitoZ_v3.4.sif all -h

# or 
$ singularity run /path/to/MitoZ_v3.4.sif -h
$ singularity run /path/to/MitoZ_v3.4.sif all -h

```

However, if you want to use the `mitoz-tools` command, you need to do it this way:
```bash
$ singularity exec /path/to/MitoZ_v3.4.sif mitoz
$ singularity exec /path/to/MitoZ_v3.4.sif mitoz-tools

# to use the `all` command of MitoZ with the `exec` command, do this:
$ singularity exec /path/to/MitoZ_v3.4.sif mitoz all -h

```

Or, you can also 'shell' into the container, as shown by **Usage 2** below. 


## 2.4 Usage 2 
```bash
$ mkdir -p /my/workdir/projectID
$ cd /my/workdir/projectID

# the below command assume your fastq files are located under the `/my/workdir/projectID` directory,
# so within Singularity' shell, you can access these fastq files directly.
$ singularity shell /path/to/MitoZ_v3.4.sif
# after login the container, it is just like you are in anther Linux machine,
# so you can use the `mitoz` command directly:
Singularity> which mitoz
/app/anaconda/bin/mitoz
Singularity> mitoz -h
Singularity> mitoz-tools -h
#
# After you finish the analysis, use the `exit` command to exit the container:
Singularity> exit

# However, if your fastq files are located at another places, say `/my/workdir/raw_data_dir`,
# to let MitoZ Singularity container can access them, you need to mirror the path into the container using the `--bind` option:
$ cd /my/workdir/projectID
$ singularity shell --bind /my/workdir/raw_data_dir  /path/to/MitoZ_v3.4.sif

```

## 2.5 Copy MitoZ's database out of the Singularity container

Do this only if you want to **customize your PCG annotation database**, see https://github.com/linzhi2013/MitoZ/wiki/Extending-MitoZ%27s-database for more details.

```bash
$ mkidr ~/mitoz_custom_db/
$ singularity shell /path/to/MitoZ_v3.4.sif
Singularity> cp -r /app/anaconda/lib/python3.9/site-packages/mitoz/profiles ~/mitoz_custom_db
Singularity> exit
# Now the `profiles` are under the `~/mitoz_custom_db` of your host machine and you can modify them to create your custom database. 
```

# 3. Conda-Pack

The installation of MitoZ via [Conda](#4-conda) often has the missing Perl module problems. If you cannot use the [Singularity images](#1-singularity) nor [Docker images](#2-docker) methods, you can try this Conda-Pack version. This method is also useful if your server is not able to access the Internet (But you still need to install the Etetoolkit taxonomy database by yourself if there is no Internet).

Here I packaged the whole conda environment (including all files) into a file named `mitoz3.4.tar.gz` using the Conda-Pack tool (https://conda.github.io/conda-pack/).

I created this environment on a Linux machine, thus it should also work on another Linux machine.

Firstly, can download the `mitoz3.4.tar.gz` file from Dropbox (https://www.dropbox.com/sh/lcrz1cufew8ormx/AAA6sR23IxajzpdEapMCtHALa?dl=0).

Next, 

```bash

# Choose a directory on your machine for the installation of MitoZ, e.g. '~/soft/mitoz3.4'
$ mkdir -p ~/soft/mitoz3.4

# then unpack mitoz3.4.tar.gz into this target directory 
$ tar -xzf /path/to/downloaded/mitoz3.4.tar.gz -C ~/soft/mitoz3.4

# Activate the environment
$ source ~/soft/mitoz3.4/bin/activate

# Cleanup prefixes from in the active environment.
# Note that this command can also be run without activating the environment
# as long as some version of python is already installed on the machine.
(mitoz3.4) $ conda-unpack

# At this point the environment is exactly as if you installed it here
# using conda directly. All scripts should work fine.
(mitoz3.4) $ mitoz -h
(mitoz3.4) $ mitoz-tools -h

# Deactivate the environment to remove it from your path when your finish MitoZ analysis
(mitoz3.4) $ source ~/soft/mitoz3.4/bin/deactivate
```

Please refer to https://conda.github.io/conda-pack/ for more details.

Now you can go to install the [Etetoolkit database](#6-the-etetoolkit-database)



# 4. Conda

## 4.1 Installation of Conda and Mamba
Firstly, install Miniconda (https://docs.conda.io/en/latest/miniconda.html) (recommended) or Anaconda (https://www.anaconda.com/products/distribution#Downloads) :

```bash
$ curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ sh Miniconda3-latest-Linux-x86_64.sh
# setup channels
$ conda config --add channels defaults
$ conda config --add channels bioconda
$ conda config --add channels conda-forge

$ conda install mamba -n base -c conda-forge  # "mamba" is much much faster than the "conda" command!

```

## 4.2 Installation of MitoZ

Conda version of MitoZ currently only fully functionally runs on Linux.

```bash
$ mamba clean -y -a # in case something unknown interferes with our installation
$ conda clean -y -a # in case something unknown interferes with our installation

# It is a good idea to install MitoZ into an independent environment, i.e. 'mitozEnv' here!
$ mamba create -n mitozEnv -c bioconda mitoz=3.4 # It's recommended to specify the version you want to install!

# Note: 
# 1. You can use any other name instead of 'mitozEnv' as the environment name, e.g. 'mitoz3.4', 
#.   so you can do 'mamba create -n mitoz3.4 -c bioconda mitoz=3.4'
#
# 2. You can also install MitoZ to a specific path, 
#    like 'mamba create -p /share/pool/guanliang/soft/mitoz3.4 -c bioconda mitoz=3.4',
#    and then use 'source activate /share/pool/guanliang/soft/mitoz3.4' to activate the environment. 


$ source activate mitozEnv   # or use "mamba" or "conda" instead of "source" the command here.
# If your MitoZ environment name is 'mitoz3.4', then you should do 'source activate mitoz3.4'!

$ circos --modules # check if all Perl modules required by circos are installed. Some modules could still be missing (don't know why conda did not fix them automatically). Similar problems can be seen at https://github.com/bioconda/bioconda-recipes/issues/9830

# Now we are ready to go:
$ mitoz # all subcommands are within this command now!

$ mitoz-tools # some useful tools for mitochondrial genome analysis
```

Now you can go to install the [Etetoolkit database](#6-the-etetoolkit-database)


## 4.3 Location of the installation

If you want to find the path where MitoZ is installed, execute:
```bash
$ conda env list
# conda environments:
#
base                  *  /home/guanliang/soft/miniconda3
mitozEnv                 /home/guanliang/soft/miniconda3/envs/mitozEnv
```
The exact path for me is: `/home/guanliang/soft/miniconda3/envs/mitozEnv/lib/python3.7/site-packages/mitoz`. For example, this is the path for MitoZ's database:
```bash
$ ll /home/guanliang/soft/miniconda3/envs/mitozEnv/lib/python3.7/site-packages/mitoz/profiles
total 16K
-rw-rw-r-- 2 guanliang    0 May 12 06:47 __init__.py
drwxrwxr-x 2 guanliang 4.0K May 24 16:06 CDS_HMM
drwxrwxr-x 2 guanliang 4.0K May 24 16:06 rRNA_CM
drwxrwxr-x 2 guanliang 4.0K May 24 16:06 __pycache__
drwxrwxr-x 2 guanliang 4.0K May 24 17:36 MT_database
```
See also  [**Extending MitoZ's database**](https://github.com/linzhi2013/MitoZ/wiki/Extending-MitoZ%27s-database).



## 4.4 Problems

**Make sure that you are the owner of the `conda`/`mamba` commands**, it happened to me that when I used another user's `conda` command I got a lot of trouble. In this case, you can follow the very beginning instruction of this page and install your own Miniconda/Anaconda.

### 4.4.1 mamba installation
```bash
$ conda install mamba -n base -c conda-forge
Collecting package metadata (current_repodata.json): done
Solving environment: / 
The environment is inconsistent, please check the package plan carefully
The following packages are causing the inconsistency:

  - defaults/linux-64::python-language-server==0.34.1=py38_0
  - defaults/noarch::python-jsonrpc-server==0.3.4=py_1                                                  
\ failed with initial frozen solve. Retrying with flexible solve.
```

Possible solutions:

```bash
$ mamba clean -y -a # in case something unknown interferes with our installation
$ conda clean -y -a # in case something unknown interferes with our installation
```
and try again.

Or, install the `mamba` into a separate environment:

```bash
$ conda create -n mambaEnv -c conda-forge mamba
# and then use the `mamba` command within this env:
$ conda activate mambaEnv
```

**Finally, you can try to install a new Miniconda** (https://docs.conda.io/en/latest/miniconda.html) (recommended) or Anaconda (https://www.anaconda.com/products/distribution#Downloads) **at a totally different place**.

Please use Google to find the solution that works for you. 

**Or, you can simply keep using the `conda` command (just replace the `mamba` with `conda`) to install MitoZ, which might cost you extra time though.**

### 4.4.2 Circos Missing Perl Modules

After the `mamba create -n mitozEnv -c bioconda mitoz` command, you should check if there are some missing Perl modules required by Circos, sometimes they are missing, and I do not know the exact reason.
 
```bash
$ source activate mitozEnv
$ circos --modules # check if all Perl modules required by circos are installed. Some modules could still be missing (don't know why conda did not fix them automatically). Similar problems can be seen at https://github.com/bioconda/bioconda-recipes/issues/9830

# For me, the Perl modules "GD" and "GD::Polyline" were missing (although conda said they have been installed already when I ran 'conda install perl-gd'), I fixed them by running the following three commands:
$ mamba install -c conda-forge pkg-config
$ mamba install -c anaconda gcc_linux-64
$ cpanm install GD
# I will try to fix the circos' problem in bioconda's MitoZ recipe file, but for the moment, please use the above solution, or try the "mitozEnv.yaml" solution below.
# You can use the ”cpanm“ command to install other missing Perl modules if necessary. 
```

- A user proposed this solution: https://github.com/linzhi2013/MitoZ/issues/152 for the missing GD module problem.
You can test this solution and see if it works, and then leave some comments at https://github.com/linzhi2013/MitoZ/issues/152, so the other users and I can know if this is a universal solution. Thanks a lot for helping to improve the software!

- Another user reported that after she installed a new Conda, the problem got solved without changing the Perl version (i.e. using 5.26).


_**Finally, if the above methods don't work for you, then don't waste your time on them, try to use the [Singularity images](#2-singularity) or [Docker images](#3-docker) instead.**_

# 5. Source codes

## 5.1 Installation of Conda and Mamba
Firstly, install Miniconda (https://docs.conda.io/en/latest/miniconda.html) (recommended) or Anaconda (https://www.anaconda.com/products/distribution#Downloads) :

```bash
$ curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ sh Miniconda3-latest-Linux-x86_64.sh
# setup channels
$ conda config --add channels defaults
$ conda config --add channels bioconda
$ conda config --add channels conda-forge

$ conda install mamba -n base -c conda-forge  # "mamba" is much much faster than the "conda" command!

```

## 5.2 Installation of dependencies

```bash
$ mamba clean -y -a # in case something unknown interferes with our installation
$ conda clean -y -a # in case something unknown interferes with our installation
$ mamba env create -n mitozEnv -f https://github.com/linzhi2013/MitoZ/releases/download/3.4/mitozEnv.yml
$ conda activate mitozEnv

# Note: 
# 1. You can use any other name instead of 'mitozEnv' as the environment name, e.g. 'mitoz3.4', 
#.   so you can do 'mamba env create -n mitoz3.4 -f https://github.com/linzhi2013/MitoZ/releases/download/3.4/mitozEnv.yml'
#
# 2. You can also install MitoZ to a specific path, 
#    like 'mamba env create -p /share/pool/guanliang/soft/mitoz3.4 -f https://github.com/linzhi2013/MitoZ/releases/download/3.4/mitozEnv.yml',
#    and then use 'source activate /share/pool/guanliang/soft/mitoz3.4' to activate the environment. 
```

## 5.3 Installation of MitoZ

```bash
# Next, please download the newest version of MitoZ source code from https://github.com/linzhi2013/MitoZ/releases/
$ pip install ./mitoz-3.4.tar.gz 
# or 
$ tar -zxvf mitoz-3.4.tar.gz
$ cd mitoz-3.4
$ python3 setup.py install

# finally, check
$ circos --modules # check if all Perl modules required by circos are installed. Some modules could still be missing (don't know why conda did not fix them automatically). Similar problems can be seen at https://github.com/bioconda/bioconda-recipes/issues/9830
```

Now you can go to install the [Etetoolkit database](#6-the-etetoolkit-database)

## 5.4 Location of the installation

If you want to find the path where MitoZ is installed, execute:
```bash
$ conda env list
# conda environments:
#
base                  *  /home/guanliang/soft/miniconda3
mitozEnv                 /home/guanliang/soft/miniconda3/envs/mitozEnv
```
The exact path for me is: `/home/guanliang/soft/miniconda3/envs/mitozEnv/lib/python3.7/site-packages/mitoz`.


**The newest version may not always be available on bioconda, because it takes time for the bioconda team to incorporate a new version of software into the bioconda channel, besides, the [Bioconda website](https://bioconda.github.io/recipes/mitoz/README.html) often does not show the latest available versions or builds. Thus, in this case, you may want to check the https://anaconda.org/bioconda/mitoz/files or use the second method for installation.**


# 6. The Etetoolkit database

- After the installation of MitoZ, you still need to install the **etetoolkit (NCBI Taxonomy) database**, especially when the automatic installation does not work for you. 

- **Warning: it is reported that a broken etetoolkit (NCBI Taxonomy) database would result in some PCGs not annotated (https://github.com/linzhi2013/MitoZ/issues/89), or MitoZ getting "Error" during the run (e.g. during the `findmitoscaf` step). Thus, please make sure this database works well before running MitoZ.** 

- It is **recommended to run the test dataset before** applying MitoZ to your own samples, just to make sure your installation is okay. See [**8. Running the test dataset**](#8-running-the-test-dataset).


## 6.1 Installation of the etetoolkit database

Unless you install MitoZ via the Docker method, otherwise you always need further to install the etetoolkit database. 

Firstly try:
```bash
$ conda activate mitozEnv

$ python3
Python 3.9.7 (default, Sep 16 2021, 08:50:36)
[Clang 10.0.0 ] :: Anaconda, Inc. on darwin
Type "help", "copyright", "credits" or "license" for more information.
>>> from ete3 import NCBITaxa
>>> ncbi = NCBITaxa()
>>> exit()
```

If you are using the Singularity image, you need to shell into the container first:

```bash
$ singularity shell /path/to/MitoZ_v3.4.sif
Singularity> python3
Singularity> Python 3.9.7 (default, Sep 16 2021, 08:50:36)
Singularity> [Clang 10.0.0 ] :: Anaconda, Inc. on darwin
Singularity> Type "help", "copyright", "credits" or "license" for more information.
Singularity> >>> from ete3 import NCBITaxa
Singularity> >>> ncbi = NCBITaxa()
Singularity> >>> exit()

Singularity> exit
```

**Now verify the database:**
```bash
$ conda activate mitozEnv
# or shell into the singularity container:
# $ singularity shell /path/to/MitoZ_v3.4.sif

$ python3
>>> from ete3 import NCBITaxa
>>> a = NCBITaxa()
>>> a.get_name_translator(["Arthropoda"])
{'Arthropoda': [6656]}
```

If the above works for you, then you are finished and can go to [**6. Running the test dataset**](#6-running-the-test-dataset). Otherwise, please read the below instructions.


If you have trouble downloading and installing the Etetoolkit database, you can download the `taxdump.tar.gz` file or my pre-built database from https://www.dropbox.com/sh/mqjqn656x41q2wb/AAD02t_kUCjNHbBgCeYpEM88a?dl=0
 or from https://pan.baidu.com/s/1YIULJ9H3BeWKcIZdMZpcuw?pwd=7r9d (提取码:7r9d).

Then execute:
```bash
$ conda activate mitozEnv

$ python3
Python 3.9.7 (default, Sep 16 2021, 08:50:36)
[Clang 10.0.0 ] :: Anaconda, Inc. on darwin
Type "help", "copyright", "credits" or "license" for more information.
>>> from ete3 import NCBITaxa
>>> ncbi = NCBITaxa(taxdump_file='/path/to/downloaded/taxdump.tar.gz')
Loading node names...
2424313 names loaded.
277227 synonyms loaded.
Loading nodes...
2424313 nodes loaded.
Linking nodes...
Tree is loaded.
Updating database: /Users/gmeng/.etetoolkit/taxa.sqlite ...
 2424000 generating entries...
Uploading to /Users/gmeng/.etetoolkit/taxa.sqlite

Inserting synonyms:      275000
Inserting taxid merges:  65000
Inserting taxids:       2420000
>>> exit()

$ ls -lhrt ~/.etetoolkit/
total 1171272
-rw-r--r--  1 gmeng  staff    12M Jun  2 11:35 taxa.sqlite.traverse.pkl
-rw-r--r--  1 gmeng  staff   558M Jun  2 11:36 taxa.sqlite
```


**However, if you got something like this:**
```bash
>>> from ete3 import NCBITaxa
>>> ncbi = NCBITaxa(taxdump_file='taxdump.tar.gz')
Loading node names...
2424313 names loaded.
277265 synonyms loaded.
Loading nodes...
2424313 nodes loaded.
Linking nodes...
Tree is loaded.
Updating database: /home/guanliang/.etetoolkit/taxa.sqlite ...
 2424000 generating entries...
Uploading to /home/guanliang/.etetoolkit/taxa.sqlite

Inserting synonyms:      75000 Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/app/anaconda/lib/python3.9/site-packages/ete3/ncbi_taxonomy/ncbiquery.py", line 106, in __init__
    self.update_taxonomy_database(taxdump_file)
  File "/app/anaconda/lib/python3.9/site-packages/ete3/ncbi_taxonomy/ncbiquery.py", line 131, in update_taxonomy_database
    update_db(self.dbfile, taxdump_file)
  File "/app/anaconda/lib/python3.9/site-packages/ete3/ncbi_taxonomy/ncbiquery.py", line 760, in update_db
    upload_data(dbfile)
  File "/app/anaconda/lib/python3.9/site-packages/ete3/ncbi_taxonomy/ncbiquery.py", line 802, in upload_data
    db.execute("INSERT INTO synonym (taxid, spname) VALUES (?, ?);", (taxid, spname))
sqlite3.IntegrityError: UNIQUE constraint failed: synonym.spname, synonym.taxid
```

There are some bugs in the ETE 3.1.1 package, you got this problem because you installed MitoZ (**<= 3.3**) via `conda` or `mamba` commands, and unfortunately, at the early builds on BioConda, I wrongly specified `ete3=3.1.1`, I should use `ete3>=3.1.2` instead. 

In this case, you can download my **pre-build version of the etetoolkit database** (filename: `etetoolkit.tgz`) from https://www.dropbox.com/sh/mqjqn656x41q2wb/AAD02t_kUCjNHbBgCeYpEM88a?dl=0
 or from https://pan.baidu.com/s/1YIULJ9H3BeWKcIZdMZpcuw?pwd=7r9d (提取码:7r9d), and then:
```bash
$ mv /path/to/etetoolkit.tgz ~
$ cd ~
$ rm -rf ~/.etetoolkit
$ tar -zxvf etetoolkit.tgz
```
**OR, you can upgrade MitoZ** 
- via the `mamba env create -n mitozEnv -f mitozEnv.yaml` method (see the beginning) 
- via the `mamba create -n mitozEnv -c bioconda mitoz=3.4` command (see the beginning) 
- get out of the `mitozEnv` environment (`conda deactivate mitozEnv`), then install an ete3 in this 'base' environment via `mamba install -c conda-forge ete3>=3.1.2`. And then, use the Python and ete3 in this 'base' environment to create the etetoolkit database by following the beginning part of [**4-the-etetoolkit-database**](#4-the-etetoolkit-database).


## 6.2 Upgrading the local database
When you want to upgrade the etetoolkit database you can check this. See http://etetoolkit.org/docs/latest/tutorial/tutorial_ncbitaxonomy.html#upgrading-the-local-database

```bash
$ python3
Python 3.9.7 (default, Sep 16 2021, 08:50:36)
[Clang 10.0.0 ] :: Anaconda, Inc. on darwin
Type "help", "copyright", "credits" or "license" for more information.
>>> from ete3 import NCBITaxa
>>> ncbi = NCBITaxa()
>>> ncbi.update_taxonomy_database()
>>> exit()
```

Or you can also download the latest file from `https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz`, 
```bash
$ wget -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

$ python3
Python 3.9.7 (default, Sep 16 2021, 08:50:36)
[Clang 10.0.0 ] :: Anaconda, Inc. on darwin
Type "help", "copyright", "credits" or "license" for more information.
>>> from ete3 import NCBITaxa
>>> ncbi = NCBITaxa(taxdump_file='/path/to/downloaded/taxdump.tar.gz')
>>> exit()
```

# 7. About the thread number and data size

- For MitoAssemble (`--assembler mitoassemble`), using 8 to 16 threads + 2 to 8 G bp fastq data is good enough, for example `--thread_number 8`, or `--thread_number 12`. A bigger thread could take a lot of RAM (e.g. 150 GB) for the assembly step.
    - More data does not necessarily mean better mitogenome
    - Too many threads do not necessarily mean faster.

- For Megahit (`--assembler megahit`), I tested with 16 threads and 15 G bp fastq data, and set `--memory 50`, and it took around 50 GB RAM, so you can use more data and threads with Megahit. 
    - While memory usage with more data (e.g. 15G bp) seems not to be a big problem for Megahit, using more data does take more time, so it is recommended to use fewer data to save time.
    - You can also increase `--memory` usage to save time if your servers have enough RAM.

- For Spades (`--assembler spades`), I did not record the RAM usage, may be similar to Megahit?



**How to check how much RAM MitoZ uses?**

You can use the `top` or `htop` (recommended; https://anaconda.org/conda-forge/htop) to check how many resources MitoZ uses if you are running MitoZ on your server; or you can use the `qstat` command if you are using an SGE cluster.



# 8. Running the test dataset

Before applying MitoZ to your own samples, it is important to run MitoZ on the test dataset.

```bash
$ mkdir ~/test
$ cd ~/test

$ wget -c https://raw.githubusercontent.com/linzhi2013/MitoZ/master/test/test.1.fq.gz 
$ wget -c https://raw.githubusercontent.com/linzhi2013/MitoZ/master/test/test.2.fq.gz

$ conda activate mitozEnv
$ mitoz all --fq1 test.1.fq.gz  --fq2 test.2.fq.gz --outprefix ttt --clade Arthropoda --data_size_for_mt_assembly 0 --assembler megahit --requiring_taxa Arthropoda
```

You can then analyze your samples by following the [Tutorial](https://github.com/linzhi2013/MitoZ/wiki/Tutorial)