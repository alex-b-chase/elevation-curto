#!/bin/bash
#$ -N concatCORE
#$ -m a
#$ -pe openmp 32
#$ -q bio

module load RAxML/8.2.10  

BASE=/bio/abchase/refDB/coregenes

ALIGNMENT=concat.aligned.coregenes.fa

cd $BASE

raxmlHPC-PTHREADS-AVX -s $ALIGNMENT -m PROTGAMMAWAG -n coregenes -x 100 -# 100 -p 4321 -f a -T 32

