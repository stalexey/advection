#!/usr/bin/gnuplot
set style data lines
set xrange [0:10]
set yrange [-0.1:3.1]
plot 'interpolationLinear.dat' lw 1
replot 'interpolationCatmullRom.dat' lw 1
replot 'interpolationMonotonicCubicFedkiw.dat' lw 1
replot 'interpolationMonotonicCubicFritschCarlson.dat' lw 1
pause -1