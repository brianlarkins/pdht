#TITLE = "Read Times vs. Matchlist Length"
#set title TITLE offset char 0, char -1
set style data linespoints
#set term pdf size 800, 600
set terminal postscript eps enhanced color 'Helvetica,16' lw 2

set output "throughput.eps"

set xlabel "Processors"
set logscale x
set xrange [12:384]
set xtics (12,24,48,96,192,384)
set ylabel "Throughput (reads/sec)"
set logscale y
#set yrange [0:60]
set key top left

#  procs * LV == total data volume
#  time / 1000 = time in seconds
# throughput = volume / seconds in reads/sec

plot 'scaling.dat' using 1:($1*10000) /($2/1000) title "10K local elems", \
     'scaling.dat' using 1:($1*20000) /($3/1000) title "20K local elems", \
     'scaling.dat' using 1:($1*30000) /($4/1000) title "30K local elems", \
     'scaling.dat' using 1:($1*40000)/($5/1000) title "40K local elems", \
     'scaling.dat' using 1:($1*40000)/($7/1000) title "50K local elems"        
