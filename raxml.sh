#!/bin/bash

#$ -S /bin/bash
#$ -N tb_raxml
#$ -pe mpi 25
#$ -R y
#$ -l h_vmem=5G
#$ -l tmem=5G
#$ -l h_rt=9000:0:0
#$ -j y
#$ -t 1-20
#$ -o /

RAXML=/raxml-ng-0.9.0/bin/raxml-ng-mpi;
MPIRUN=/openmpi-3.1.1/bin/mpirun;

lin=l4

IN_DIR=/alignment;
OUT_DIR=/phylogenetics/$lin.wg;

FASTA=$IN_DIR/$lin.aln.fa;
$RAXML --parse --msa $FASTA --model GTR+G --prefix $OUT_DIR/$lin.aln --seed $RANDOM --precision 12 --blmin 0.000000001

IN_FILE=$OUT_DIR/$lin.aln.raxml.rba;
PREFIX=$lin.aln.${SGE_TASK_ID};

if [[ $SGE_TASK_ID -le 10 ]]; then
  $MPIRUN --map-by node -np $NSLOTS --mca btl tcp,self $RAXML --msa $IN_FILE --model GTR+G --prefix $OUT_DIR/$PREFIX --tree rand{1} --threads 1 --seed $RANDOM --precision 12 --blmin 0.000000001;
else
  $MPIRUN --map-by node -np $NSLOTS --mca btl tcp,self $RAXML --msa $IN_FILE --model GTR+G --prefix $OUT_DIR/$PREFIX --tree pars{1} --threads 1 --seed $RANDOM --precision 12 --blmin 0.000000001;
fi;
