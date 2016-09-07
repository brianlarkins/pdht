TITLE = "Read Times vs. Matchlist Length"
set title TITLE offset char 0, char -1
set style data linespoints
set term png size 800, 600
set output "mlen.png"

set xlabel "Matchlist Length"
set logscale x
set xrange [100:100000]
set ylabel "Time (ms)"

plot 'matchlength.dat' using 1:2 title "8 bytes local", \
     'matchlength.dat' using 1:3 title "8 bytes remote", \
     'matchlength.dat' using 1:4 title "128 bytes remote", \
     'matchlength.dat' using 1:5 title "1K bytes remote"        
