#!/usr/bin/gnuplot
set style data lines
set xrange [0:10]
set yrange [-0.5:1.5]
plot 'advectionLinear.dat' lw 2
replot 'advectionCatmullRom.dat' lw 2
replot 'advectionMonotonicCubicFedkiw.dat' lw 2
replot 'advectionMonotonicCubicFritschCarlson.dat' lw 2
pause -1