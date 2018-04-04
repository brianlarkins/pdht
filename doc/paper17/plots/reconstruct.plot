#TITLE = "Read Times vs. Matchlist Length"
#set title TITLE offset char 0, char -1
set style data linespoints
#set term pdf size 800, 600
set terminal postscript eps enhanced color 'Helvetica,24' lw 2

set output "reconstruct.eps"

set xlabel "Processes"
set logscale x
#set xrange [100:100000]
set xtics (24,  40, 60, 120, 160, 240)
set ylabel "Execution Time (sec)"
set yrange [0:10]
set ytics (1,5,10)
set key top left

plot 'madness.dat' index 0 using 1:4 title "PDHT", \
     'madness.dat' index 1 using 1:4 title "MVAPICH2"
