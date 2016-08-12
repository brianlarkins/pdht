#!/bin/sh

for (( j = 1; j <= 10; j++ ))
do
  echo "m,n,r,c,num_nonzero_blocks,num_trials,median_time,min_time,max_time,mflops,num_loads,num_stores" > ../sizedata/sizedata_itanium2_$j.csv

  for (( i = 100; i <= 2000; i += 100 ))
  do
      ../lowlevel_tester $i $i 1 1 -$j 5 ../sizedata/sizedata_itanium2_$j.csv
  done
done
