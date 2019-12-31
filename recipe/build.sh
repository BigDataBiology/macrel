#!/bin/bash

outdir=$PREFIX/share/$PKG_NAME-$PKG_VERSION-$PKG_BUILDNUM
mkdir -p $outdir
mkdir -p $PREFIX/bin
mkdir -p $outdir/envs/FACS_env/bin

# compile c
cd prodigal_modified
make
chmod +x prodigal
cd ..

cp $SRC_DIR/prodigal_modified/prodigal $outdir/envs/FACS_env/bin/prodigal_sm
cp $SRC_DIR/CTDDClass.py $outdir/CTDDClass.py
cp $SRC_DIR/FACS.sh $outdir/FACS.sh
cp $SRC_DIR/LICENSE $outdir/LICENSE
cp $SRC_DIR/Predict_130819.R $outdir/Predict_130819.R
cp $SRC_DIR/features_130819.R $outdir/features_130819.R
cp $SRC_DIR/r22_largeTraining.rds $outdir/r22_largeTraining.rds
cp $SRC_DIR/rf_dataset1.rds $outdir/rf_dataset1.rds


ln -s $outdir/* $PREFIX/bin/

