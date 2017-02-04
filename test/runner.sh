#!/bin/bash
time=20
ppn=12
cpt=2
iter=1
vol=30000

mkdir -p runs.$ppn

#for nodes in 1 2 4 8 16 32
for nodes in 8 
#for nodes in 1 2 4 8
do
  echo running $nodes
  srun -t $time:0 -c $cpt --ntasks-per-node=$ppn -N $nodes ./scaling -v $vol -i $iter > runs.$ppn/$((ppn * nodes)).$((vol * iter))K.8bytes.out
done
