#TITLE = "Read Times vs. Matchlist Length"
#set title TITLE offset char 0, char -1
set style data linespoints
#set term pdf size 800, 600
set terminal postscript eps enhanced color 'Helvetica,24' lw 2
set pointsize 1.5

set output "pdhtlatency.eps"

set xlabel "Number of Entries"
set logscale x
#set xrange [100:100000]
#set xtics (1,10,100,"1k" 1000,"10k" 10000, "25k" 25000, "100k" 100000, "" 150000)
set ylabel "Read Latency (usec)"
set logscale y
#set yrange [.3:850]
#set ytics (1,5,10,50,100,500)
set key top left

plot 'latency.dat' index 0 using 1:3 title "matchlist" lw 2, \
     'latency.dat' index 1 using 1:3 title "multi PTE (5)" lw 2, \
     'latency.dat' index 2 using 1:3 title "unordered" lw 2, \
     'latency.dat' index 2 using 1:2 title "local" lw 2
