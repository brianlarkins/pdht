#TITLE = "Read Times vs. Matchlist Length"
#set title TITLE offset char 0, char -1
set style data linespoints
#set term pdf size 800, 600
set terminal postscript eps enhanced color 'Helvetica,24' lw 2

set output "meraculous.eps"

set xlabel "Processes"
set logscale x
#set xrange [100:100000]
set xtics (10, 20, 40, 60, 80, 120)
set ylabel "Execution Time (sec)"
set yrange [0:175]
set ytics (0, 25, 50, 75 ,100, 125, 150, 175)
set key top left

plot 'shepard.meraculous.dat' index 0 using 1:3 title "GCC-UPC"
#     'shepard.meraculous.dat' index 0 using 1:3 title "MVAPICH2"
