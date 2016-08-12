#!/bin/sh

echo "m,n,r,c,num_nonzero_blocks,num_trials,median_time,min_time,max_time,mflops,num_loads,num_stores" > ../nnz_per_row_itanium2.csv

for (( j = 1; j <= 200; j++ ))
do
    ../lowlevel_tester 500 500 1 1 -$j 10 ../nnz_per_row_itanium2.csv
done
