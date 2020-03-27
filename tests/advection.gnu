#!/usr/bin/gnuplot
set style data lines
set xrange [0:10]
set yrange [-0.5:2]
plot 'advectionReference.dat' lw 1 lt rgb "gray"
replot 'advectionLinear.dat' lw 1 lt rgb "orange"
replot 'advectionCatmullRom.dat' lw 1 lt rgb "blue"
replot 'advectionMonotonicCubicFedkiw.dat' lw 1 lt rgb "pink"
replot 'advectionMonotonicCubicFritschCarlson.dat' lw 1 lt rgb "violet"
replot 'advectionLaxWendroffCDF.dat' lw 1 lt rgb "cyan"
pause -1