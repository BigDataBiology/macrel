#!/usr/bin/env bash

set -e

# This is set in travis
if [[ x$PYTHON_VERSION == x ]]; then
    PYTHON_VERSION=3.7
fi


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
conda create --yes -p $BASEDIR/envs/Macrel_env python=$PYTHON_VERSION
conda config --env --add channels defaults
conda config --env --add channels bioconda
conda config --env --add channels conda-forge

if ! which mamba > /dev/null; then
    CONDA_INSTALL_CMD=conda
else
    echo "# Installing packages with MAMBA"
    CONDA_INSTALL_CMD=mamba
fi

${CONDA_INSTALL_CMD} install -y \
          --prefix $BASEDIR/envs/Macrel_env \
          ngless \
          pyrodigal>=0.7.3 \
          megahit \
          paladin \
          pandas \
          scikit-learn \
          atomicwrites \
          tzlocal \
          make \
          c-compiler

# only activate AFTER install all the packages to get the right environmental variables
source activate $BASEDIR/envs/Macrel_env

# env forwards the current environment and sets some extra variables needed by
# ./recipe/build.sh
env \
        PYTHON=$(which python) \
        SRC_DIR=$PWD \
        PREFIX=$CONDA_PREFIX \
        ./recipe/build.sh


echo "############ Installation procedures finished
****** Thank you for installing Macrel ********
--- Please submit bugreports/comments to
celio.diasjunior@gmail.com
or
https://github.com/BigDataBiology/macrel/issues
"
