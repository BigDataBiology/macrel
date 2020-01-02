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
for f in \
        FACS.sh \
        LICENSE \
        Predict_130819.R \
        AMP-features.py \
        r22_largeTraining.rds \
        rf_dataset1.rds \
        trim.pe.ngl \
        trim.se.ngl \
        ; do
    cp $SRC_DIR/$f $outdir/$f
done

ln -s $outdir/* $PREFIX/bin/

