#!/usr/bin/gnuplot
set style data lines
set xrange [0:10]
set yrange [-0.1:3.1]
plot 'interpolationLinear.dat' lw 1 lt rgb "orange"
replot 'interpolationCatmullRom.dat' lw 1 lt rgb "blue"
replot 'interpolationMonotonicCubicFedkiw.dat' lw 1 lt rgb "pink"
replot 'interpolationMonotonicCubicFritschCarlson.dat' lw 1 lt rgb "violet"
pause -1