#!/bin/bash
time=20
ppn=12
cpt=2
vol=10000

mkdir -p runs.norm

for iter in 1 2 5 10
do
  echo "$((iter * 10)) K runs"
  for nodes in 1 2 4 8 16 32
  #for nodes in 1 2 4 8
  do
    echo running $nodes
    srun -t $time:0 -c $cpt --ntasks-per-node=$ppn -N $nodes ./scaling -v $vol -i $iter > runs.norm/$((ppn * nodes)).$((vol * iter))K.8bytes.out
  done
done
