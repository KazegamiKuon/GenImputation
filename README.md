# GenImputation

This is template repository. Structure and automated workflow for AI project.

Programming language is **python**. However for easy test, demo and debug, I usualy code in by [jupyter notebook](https://jupyter.org/install) on [anaconda](https://docs.anaconda.com/anaconda/install/)

You should install anconda before jupyter

Anyway for begin, we start with python using terminal and git

## Who are you?

If  you configed your git then pass.

Check your config:

```script
git config --list --show-origin
```

If your **email** and **name** already existed then you configed or you can change your configuration too

```script
git config --global user.name your_user_name
git config --global user.email your_email
```

## Create new repository

You create repository on github and then follow the directions on github.

**Remember**, you should create your folder with name same your repository name and cd to that folder before run script.

## Create workflow

For easy quality control and merge, I highly recomend you should follow one workflow. I follow this [workflow](https://nvie.com/posts/a-successful-git-branching-model/) also add in something. Check my workflow to see more.

You must create **dev** *branch* from **master** or **main** *branch* before starting. Highly recomend you create it on github. Then change your branch.

```script
git checkout dev
```

Every time when you want to create a branch, please create from develop branch. After you have created, change branch again.

## Install/Uninstall your environment

I make file bash for easy do it.

**setup-library.sh** will install some necessary reponsitory and setup conda environment.

```script
source setup-library.sh
```

If the environment already exists, **samtools** and **bcftools** should installed and avaible.
Run below script to check it available or not.

```script
bcftools
samtools
```

If you didn't see library document and shown some error, Use below script to fix it or see **Samtools, bcftools other way install**.

```script
conda install -c bioconda samtools=1.9 --force-reinstall
conda install -c bioconda bcftools=1.9 --force-reinstall
```

For some reason, If bcftools still can't activate, you should run this script:

```script
conda install openssl=1.0
```

## Samtools, bcftools other way install

Activate environment and setup required [samtools](http://www.htslib.org/download/). Download **htslib**, **samtools** and **bcftools**.Unless have it, then run script below to setup.

```script
wget https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2
tar xjf htslib-1.12.tar.bz2
cd ./htslib-1.12/
sudo bash ./configure
make
sudo make install
cd ../
wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2
tar xif samtools-1.12.tar.bz2
cd ./samtools-1.12/
sudo bash ./configure
make
sudo make install
cd ../
wget https://github.com/samtools/bcftools/releases/download/1.12/bcftools-1.12.tar.bz2
tar xif bcftools-1.12.tar.bz2
cd bcftools-1.12/
sudo bash ./configure
make
sudo make install
```

## (coming soon)
