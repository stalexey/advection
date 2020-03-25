#!/usr/bin/gnuplot
set style data lines
set xrange [0:10]
set yrange [-0.5:1.5]
plot 'interpolationLinear.dat' lw 2
replot 'interpolationCatmullRom.dat' lw 2
replot 'interpolationMonotonicCubic.dat' lw 2
pause -1