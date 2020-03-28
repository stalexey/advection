#!/usr/bin/gnuplot
set style data lines
set xrange [0:10]
set yrange [-0.1:3.1]
plot 'interpolation_Linear.dat' lw 1 lt rgb "orange"
replot 'interpolation_CatmullRom.dat' lw 1 lt rgb "blue"
replot 'interpolation_MonotonicCubicFedkiw.dat' lw 1 lt rgb "pink"
replot 'interpolation_MonotonicCubicFritschCarlson.dat' lw 1 lt rgb "violet"
pause -1