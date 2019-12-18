#!/usr/bin/env bash

echo "
#########################################################################
######## FACS pipeline - Fast AMP Clustering System   ###################
#########################################################################
########                    Authors:                  ###################
########                                              ###################
######## Célio Dias Santos Júnior, Luis Pedro Coelho  ###################
#########################################################################
######## ISTBI - FUDAN University / Shanghai - China  ###################
#########################################################################
"

if ! which conda > /dev/null; then
    echo "[ Conda not found. Please install miniconda and add 'conda' to the PATH: "
    echo "    curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    echo "    sh Miniconda3-latest-Linux-x86_64.sh"
    exit 1
fi

eval "$(conda shell.bash hook)"

echo "[ ## 1.] Installing routine linux softwares"
echo "# Creating new environment for FACS"
mkdir -p envs
conda create -p ./envs/FACS_env python=3.7
conda activate ./envs/FACS_env
conda config --env --add channels r
conda config --env --add channels defaults
conda config --env --add channels bioconda
conda config --env --add channels conda-forge

conda config --env --set auto_activate_base false

echo "# Installing Bioinformatics packages"
conda install -y sqlite
conda install -y trimmomatic
conda install -y megahit
conda install -y paladin
conda install -y samtools
conda install -y eXpress
conda install -y pigz
conda install -y parallel

echo "# Installing Python packages"
conda install -y matplotlib
conda install -y scikit-learn
conda install -y pandas

echo "# Installing R packages"
conda install -y r-essentials
conda install -y r-base
conda install -y r-caret
conda install -y r-randomforest
conda install -y r-dplyr
conda install -y r-data.table

echo "# Installing non-conda R packages"
echo "
##########################################################################
install.packages(\"Peptides\", repos = \"http://cran.us.r-project.org\", dependencies=TRUE)
install.packages(\"doParallel\", repos = \"http://cran.us.r-project.org\", dependencies=TRUE)
install.packages(\"obliqueRF\", repos = \"http://cran.us.r-project.org\", dependencies=TRUE)
##########################################################################
" > inst.R
R --vanilla --slave < inst.R
rm -rf inst.R

conda deactivate

echo "[ ## 2.] Installing prodigal_modified"

(cd prodigal_modified && make && mv prodigal ../envs/FACS_env/bin/prodigal_sm )

#echo "[ ## 3.] Getting python scripts"
#echo "# Downloading"
#curl -O https://raw.githubusercontent.com/Superzchen/iLearn/master/descproteins/CTDDClass.py
#curl -O https://raw.githubusercontent.com/Superzchen/iLearn/master/descproteins/readFasta.py
#curl -O https://raw.githubusercontent.com/Superzchen/iLearn/master/descproteins/saveCode.py

#chmod +x CTDDClass.py readFasta.py saveCode.py

echo "############ Installation procedures finished
****** Thank you for installing FACS ********
--- If any bug appears report to:
celio.diasjunior@gmail.com
"
