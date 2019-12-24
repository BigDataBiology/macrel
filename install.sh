#!/usr/bin/env bash

set -e

echo "
#########################################################################
######## FACS - Fast AMP Classification System        ###################
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

echo "# Creating new environment for FACS"
mkdir -p envs
conda create --yes -p ./envs/FACS_env python=3.7
source activate ./envs/FACS_env
conda config --env --add channels r
conda config --env --add channels defaults
conda config --env --add channels bioconda
conda config --env --add channels conda-forge

echo "# Installing conda packages"
conda install -y \
		sqlite \
        trimmomatic \
        megahit \
        paladin \
        samtools \
        eXpress \
        pigz \
        parallel \
        matplotlib \
        scikit-learn \
        pandas \
        r-essentials \
        r-base \
        r-caret \
        r-randomforest \
        r-dplyr \
        r-data.table \
		--quiet

echo "# Installing non-conda R packages"
echo "
##########################################################################
install.packages(\"Peptides\", repos = \"http://cran.us.r-project.org\", dependencies=TRUE)
install.packages(\"doParallel\", repos = \"http://cran.us.r-project.org\", dependencies=TRUE)
install.packages(\"obliqueRF\", repos = \"http://cran.us.r-project.org\", dependencies=TRUE)
##########################################################################
" > inst.R
R --vanilla --slave < inst.R --quiet
rm -rf inst.R

echo "[ ## 2.] Installing prodigal_modified"

cd prodigal_modified
make --quiet # conda will add $GCC to environment
mv prodigal ../envs/FACS_env/bin/prodigal_sm

source deactivate

cd ..
chmod +x CTDDClass.py

echo "############ Installation procedures finished
****** Thank you for installing FACS ********
--- Please submit bugreports/comments to
celio.diasjunior@gmail.com
"
