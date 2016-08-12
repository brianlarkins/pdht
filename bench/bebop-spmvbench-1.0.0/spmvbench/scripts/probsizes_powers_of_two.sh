#!/bin/sh

for r in 2 3 4 6 8
do

for c in 2 3 4 6 8
do

for ((j = 50; j <= 60; j++ ))
do

echo "m,n,r,c,num_nonzero_blocks,num_trials,median_time,min_time,max_time,mflops,num_loads,num_stores" > ../sizedata_powers_of_two/sizedata_opteron_"$r"x"$c"_$j.csv

for (( i = 1024; i <= 5000000; i *= 2 ))
do
    ../lowlevel_tester $i $i $r $c -$j 10 ../sizedata_powers_of_two/sizedata_opteron_"$r"x"$c"_$j.csv
done

done

done

done
