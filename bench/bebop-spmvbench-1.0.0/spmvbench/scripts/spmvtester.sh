#!/bin/sh

echo "m,n,r,c,num_nonzero_blocks,num_trials,median_time,min_time,max_time,mflops,num_loads,num_stores" > ../benchdata_itanium2_small.csv

for i in 1 2 3 4 6 8 12
do
  for j in 1 2 3 4 6 8 12
  do
    ../lowlevel_tester 407 407 $i $j .3 5 ../benchdata_itanium2_small.csv
  done
done
