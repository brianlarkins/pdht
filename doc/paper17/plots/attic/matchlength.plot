#TITLE = "Read Times vs. Matchlist Length"
#set title TITLE offset char 0, char -1
set style data linespoints
#set term pdf size 800, 600
set terminal postscript eps enhanced color 'Helvetica,24' lw 2

set output "mlen.eps"

set xlabel "Matchlist Length"
set logscale x
#set xrange [100:100000]
set xtics (1,10,100,1000,10000,100000)
set ylabel "Read Latency (usec)"
set logscale y
set yrange [5:5000]
set key top left

#plot 'matchlength.dat' using 1:2 title "8 bytes local", \
#     'matchlength.dat' using 1:3 title "8 bytes shared", \

plot 'matchlength.dat' using 1:4 title "8 list", \
     'matchlength.dat' using 1:6 title "1K list", \
     'matchlength.dat' using 1:7 title "8 hashed", \
     'matchlength.dat' using 1:9 title "1K hashed"