# Gnuplot script to plot impulse response and exact impulse from "impulse.txt"
# Pascal Getreuer 2013

set multiplot layout 2,1
set xlabel "n"

# Impulse Response
plot "impulse.txt" u 1:2 w linespoints linestyle 3 ps 0.5 title "output", \
    "impulse.txt" u 1:3 w linespoints linestyle 8 ps 0.5 title "exact"

# Impulse Response Error
plot "impulse.txt" u 1:($2 - $3) w linespoints linestyle 4 ps 0.5 title "error"

unset multiplot
pause -1
