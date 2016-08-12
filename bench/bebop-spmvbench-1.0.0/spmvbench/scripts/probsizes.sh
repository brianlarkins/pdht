#!/bin/sh


for ((j = 1; j <= 25; j++ ))
do

echo "m,n,r,c,num_nonzero_blocks,num_trials,median_time,min_time,max_time,mflops,num_loads,num_stores" > ../sizedata_per_row/sizedata_opteron_$j.csv

for (( i = 100; i <= 1000; i += 100 ))
do
    ../lowlevel_tester $i $i 1 1 -$j 10 ../sizedata_per_row/sizedata_opteron_$j.csv
done

for (( i = 2000; i <= 9000; i += 1000 ))
do
    ../lowlevel_tester $i $i 1 1 -$j 10 ../sizedata_per_row/sizedata_opteron_$j.csv
done

for (( i = 10000; i <= 100000; i += 10000 ))
do
    ../lowlevel_tester $i $i 1 1 -$j 10 ../sizedata_per_row/sizedata_opteron_$j.csv
done

#for (( i = 110000; i <= 200000; i += 10000 ))
#do
#   ../lowlevel_tester $i $i 1 1 $j 10 ../sizedata/sizedata_opteron_$j.csv
#done

done
