#!/usr/bin/env bash

export PATH=$HOME/miniconda3/bin:$PATH

if test -e $HOME/miniconda3/bin; then
    echo "miniconda already installed."
else
    rm -rf $HOME/miniconda3

    wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda3.sh
    chmod +x miniconda3.sh
    ./miniconda3.sh -b -p $HOME/miniconda3
    conda update --yes conda --quiet
fi

# For debugging:
conda info -a

