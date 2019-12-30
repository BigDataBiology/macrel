#!/bin/bash

conda config --env --add channels defaults
conda config --env --add channels r
conda config --env --add channels bioconda
conda config --env --add channels conda-forge

outdir=$PREFIX/share/$PKG_NAME-$PKG_VERSION-$PKG_BUILDNUM
mkdir -p $outdir
mkdir -p $PREFIX/bin
mkdir -p $outdir/envs/FACS_env/bin

# compile c
cd prodigal_modified
make
chmod +x prodigal
cd ..

# save source code
cp -r $SRC_DIR/example_seqs/ $outdir/example_seqs/
cp $SRC_DIR/prodigal_modified/prodigal $outdir/envs/FACS_env/bin/prodigal_sm
cp $SRC_DIR/CTDDClass.py $outdir/CTDDClass.py
cp $SRC_DIR/FACS.sh $outdir/FACS.sh
cp $SRC_DIR/LICENSE $outdir/LICENSE
cp $SRC_DIR/Predict_130819.R $outdir/Predict_130819.R
cp $SRC_DIR/features_130819.R $outdir/features_130819.R
cp $SRC_DIR/r22_largeTraining.rds $outdir/r22_largeTraining.rds
cp $SRC_DIR/rf_dataset1.rds $outdir/rf_dataset1.rds


# soft connection
ln -s $outdir/* $PREFIX/bin/

# start script
#cp $RECIPE_DIR/FACS.py $outdir/FACS
#chmod +x $outdir/FACS

#soft connection
#ln -s $outdir/FACS $PREFIX/bin