#TITLE = "Read Times vs. Matchlist Length"
#set title TITLE offset char 0, char -1
set style data linespoints
#set term pdf size 800, 600
set terminal postscript eps enhanced color 'Helvetica,24' lw 2

set output "scaling64.eps"

set xlabel "Processes"
#set logscale x
set xrange [10:512]
set xtics (12, "" 24, "" 48, 60, 120, 180, 240, 288, 360, 480)
set ylabel "Throughput (reads/sec)"
#set logscale y
#set yrange [.1:20]
set key top left

#  procs * LV == total data volume
#  time / 1000 = time in seconds
# throughput = volume / seconds in reads/sec

# plotting #procs vs. (#procs * 1000 reads each) / (time-in-ms / 1000)
plot 'comet.scaling.dat' index 0 using 1:($1*1000)/($2/1000) title "MVAPICH2", \
     'comet.scaling.dat' index 0 using 1:($1*1000)/($4/1000) title "PDHT"
