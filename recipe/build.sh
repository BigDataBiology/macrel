#!/bin/bash

set -e

mkdir -p $PREFIX/bin

# compile c
cd prodigal_modified
make
chmod +x prodigal
cd ..

cp $SRC_DIR/prodigal_modified/prodigal $PREFIX/bin/prodigal_sm
python setup.py install

