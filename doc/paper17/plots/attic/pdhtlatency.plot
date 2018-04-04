#TITLE = "Read Times vs. Matchlist Length"
#set title TITLE offset char 0, char -1
set style data linespoints
#set term pdf size 800, 600
set terminal postscript eps enhanced color 'Helvetica,24' lw 2

set output "pdhtlatency.eps"

set xlabel "# Entries"
set logscale x
#set xrange [100:100000]
set xtics (1,10,100,1000,10000,50000)
set ylabel "Read Latency (usec)"
set logscale y
set yrange [.5:850]
set ytics (1,3,6,10,50,100,400)
set key top left

plot 'pdhtlatency.dat' index 0 using 1:2 title "local", \
     'pdhtlatency.dat' index 0 using 1:3 title "unordered", \
     'pdhtlatency.dat' index 1 using 1:3 title "multi PTE (5)", \
     'pdhtlatency.dat' index 2 using 1:3 title "matchlist"
