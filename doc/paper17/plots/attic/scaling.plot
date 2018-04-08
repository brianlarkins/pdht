#TITLE = "Read Times vs. Matchlist Length"
#set title TITLE offset char 0, char -1
set style data linespoints
#set term pdf size 800, 600
set terminal postscript eps enhanced color 'Helvetica,16' lw 2

set output "throughput.eps"

set xlabel "Nodes"
set logscale x
set xrange [1:76]
set xtics (1,2,4,8,16,32,64,72)
set ylabel "Throughput (reads/sec)"
#set logscale y
#set yrange [0:60]
set key top left

#  procs * LV == total data volume
#  time / 1000 = time in seconds
# throughput = volume / seconds in reads/sec

plot 'scaling.dat' using 1:($1*10000) /($2/1000) title "LV 10K", \
     'scaling.dat' using 1:($1*20000) /($3/1000) title "LV 20K", \
     'scaling.dat' using 1:($1*50000) /($4/1000) title "LV 50K", \
     'scaling.dat' using 1:($1*100000)/($5/1000) title "LV 100K"        
