#!/usr/bin/env bash

# Entering in the library folder
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd $DIR
# Install FAST
sudo perl -MCPAN -e 'install FAST'
# Install pigz
sudo apt-get install zlib1g zlib1g-dev libpthread-stubs0-dev 
wget https://zlib.net/pigz/pigz-2.4.tar.gz
tar -zxvf pigz-2.4.tar.gz
cd pigz-2.4
make
mv pigz ../
mv unpigz ../
cd ../
rm -rf pigz-2.4/ pigz-2.4.tar.gz
# Install pandaseq
sudo apt-add-repository ppa:neufeldlab/ppa && sudo apt-get update && sudo apt-get install pandaseq
# Install Trimmomatic
sudo apt-get install openjdk-11-jdk-headless git
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
rm -rf Trimmomatic-0.39.zip
# Install ORFM
wget https://github.com/wwood/OrfM/releases/download/v0.7.1/orfm-0.7.1_Linux_x86_64.tar.gz
tar -zxf orfm-0.7.1_Linux_x86_64.tar.gz
mv orfm-0.7.1_Linux_x86_64/orfm-0.7.1_Linux_x86_64 ./orfm
rm -rf orfm-0.7.1_Linux_x86_64 orfm-0.7.1_Linux_x86_64.tar.gz
# Install CTDDClass
## Note you need to have installed python3 and the following packages:
## sys, os, shutil, scipy, argparse, collections, platform, math, re, numpy (1.13.1), sklearn (0.19.1), matplotlib (2.1.0), and pandas (0.20.1)
git clone https://github.com/Superzchen/iLearn
rsync -avzP iLearn/descproteins/CTDDClass.py ./
rsync -avzP iLearn/descproteins/saveCode.py* ./
rsync -avzP iLearn/descproteins/readFasta.py* ./
rm -rf iLearn/
# Install PLASS
 wget https://mmseqs.com/plass/plass-static_sse41.tar.gz; tar xvfz plass-static_sse41.tar.gz; export PATH=$(pwd)/plass/bin/:$PATH
# Install PALADIN
git clone https://github.com/twestbrookunh/paladin.git
cd paladin/
make
# Install samtools
conda install -c bioconda samtools
# Install MMSeqs2
 wget https://mmseqs.com/latest/mmseqs-linux-sse41.tar.gz; tar xvfz mmseqs-linux-sse41.tar.gz; export PATH=$(pwd)/mmseqs2/bin/:$PATH
# Install eXpress
wget https://pachterlab.github.io/eXpress/downloads/express-1.5.1/express-1.5.1-linux_x86_64.tgz
tar -zxf express-1.5.1-linux_x86_64.tgz
# Making sure executability
chmod +x *.py *.sh *.R
