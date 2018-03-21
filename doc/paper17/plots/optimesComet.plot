#TITLE = "Read Times vs. Matchlist Length"
#set title TITLE offset char 0, char -1
set style data linespoints
#set term pdf size 800, 600
set terminal postscript eps enhanced color 'Helvetica,24' lw 2

set output "mpilatency.eps"

set xlabel "# Entries"
set logscale x
#set xrange [100:100000]
set xtics (1,10,100,1000,10000,50000)
set ylabel "Read Latency (usec)"
#set logscale y
set yrange [.1:13.5]
set key top left

plot 'optimesComet.dat' index 0 using 1:2 title "PDHT local Ordered", \
     'optimesComet.dat' index 1 using 1:2 title "PDHT remote Ordered", \
     'optimesComet.dat' index 2 using 1:2 title "PDHT local Unordered", \
     'optimesComet.dat' index 3 using 1:2 title "PDHT remote Unordered", \
     'optimesComet.dat' index 4 using 1:2 title "PDHT local searching Unordered", \
     'optimesComet.dat' index 5 using 1:2 title "PDHT MPI local", \
     'optimesComet.dat' index 6 using 1:2 title "PDHT MPI remote"

