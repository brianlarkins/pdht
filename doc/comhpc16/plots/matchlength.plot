TITLE = "Read Times vs. Matchlist Length"
set title TITLE offset char 0, char -1
set style data linespoints
#set term pdf size 800, 600
set terminal postscript eps enhanced color 'Helvetica,24' lw 2

set output "mlen.eps"

set xlabel "Matchlist Length"
set logscale x
set xrange [100:100000]
set ylabel "Time (ms)"
set logscale y
set key top left

plot 'matchlength.dat' using 1:2 title "8 bytes local", \
     'matchlength.dat' using 1:3 title "8 bytes remote", \
     'matchlength.dat' using 1:4 title "128 bytes remote", \
     'matchlength.dat' using 1:5 title "1K bytes remote"        
