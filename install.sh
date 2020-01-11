#!/usr/bin/env bash

set -e

echo "
# Macrel - (Meta)genomic AMP Classification and Retrieval

AUTHORS: Célio Dias Santos Júnior, Luis Pedro Coelho

"

if ! which conda > /dev/null; then
    echo "[ Conda not found. Please install miniconda and add 'conda' to the PATH: "
    echo "    curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    echo "    sh Miniconda3-latest-Linux-x86_64.sh"
    exit 1
fi

eval "$(conda shell.bash hook)"

BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

echo "# Creating new environment for Macrel"
mkdir -p envs
conda create --yes -p $BASEDIR/envs/Macrel_env python=3.7
source activate $BASEDIR/envs/Macrel_env
conda config --env --add channels defaults
conda config --env --add channels bioconda
conda config --env --add channels conda-forge

echo "# Installing conda packages"

# tzlocal is necessary for rpy2-pandas integration
# rpy2 must be 3.0.0 or later!
conda install -y \
        ngless \
        megahit \
        paladin \
        pandas \
        "rpy2 > 3" \
        tzlocal \
        r-base \
        r-essentials \
        r-peptides \
        r-randomforest \
        --quiet

echo "# Installing modified version of prodigal (prodigal_sm)"
cd prodigal_modified
make CC=$GCC --quiet # conda will add $GCC to environment
cp -pir prodigal ../envs/Macrel_env/bin/prodigal_sm
cd ..

python setup.py install


echo "############ Installation procedures finished
****** Thank you for installing Macrel ********
--- Please submit bugreports/comments to
celio.diasjunior@gmail.com

or

https://github.com/BigDataBiology/macrel/issues
"
