#TITLE = "Read Times vs. Matchlist Length"
#set title TITLE offset char 0, char -1
set style data linespoints
#set term pdf size 800, 600
set terminal postscript eps enhanced color 'Helvetica,24' lw 2

set output "mpilatency.eps"

set xlabel "# Entries"
set logscale x
#set xrange [100:100000]
set xtics (1,10,100,"1k" 1000,"10k" 10000,"25k" 25000, "100k" 100000, "" 150000)
set ylabel "Read Latency (usec)"
#set logscale y
set yrange [.1:20]
set key top left

plot 'latency.dat' index 3 using 1:3 title "OpenMPI/Portals", \
     'latency.dat' index 5 using 1:3 title "MPICH/Portals", \
     'latency.dat' index 2 using 1:3 title "PDHT remote", \
     'latency.dat' index 2 using 1:2 title "PDHT local", \
     'latency.dat' index 4 using 1:3 title "MVAPICH2 MPI"
