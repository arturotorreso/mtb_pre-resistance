#!/bin/bash

#$ -S /bin/bash
#$ -N bactdating_raxml
#$ -l h_vmem=10G
#$ -l tmem=10G
#$ -l h_rt=5000:0:0
#$ -j y
#$ -t 1-3
#$ -l tscratch=3G
#$ -o /

date
hostname

lin=l4
OUT_DIR=/phylogenetics/bactDating/$lin.wg;
IN_DIR=/phylogenetics/$lin.wg;
IN_FILE=$IN_DIR/$lin.aln.dates.tree;

LENGTH=4411532

MODEL=mixedgamma
ITER=1e7
RATE=0.5

PREFIX=$lin.$MODEL.$ITER.$RATE.wg${SGE_TASK_ID}


Rscript --vanilla /SAN/ballouxlab/tb_lima/dev/runBactDat.R -i $IN_FILE -p $PREFIX -d $OUT_DIR -l $LENGTH -m $MODEL -n $ITER -s TRUE -M $RATE > $OUT_DIR/${PREFIX}_simplified.txt

date
hostname
