#!/usr/bin/gnuplot
set style data lines
set xrange [0:10]
set yrange [-0.5:2]
plot 'advection_Reference.dat' lw 1
replot 'advection_SemiLagrangian_Linear.dat' lw 1
replot 'advection_SemiLagrangian_CatmullRom.dat' lw 1
#replot 'advection_SemiLagrangian_MonotonicCubicFedkiw.dat' lw 1
replot 'advection_SemiLagrangian_MonotonicCubicFritschCarlson.dat' lw 1
replot 'advection_MacCormack_Linear.dat' lw 1
#replot 'advection_MacCormack_CatmullRom.dat' lw 1
#replot 'advection_MacCormack_MonotonicCubicFedkiw.dat' lw 1
#replot 'advection_MacCormack_MonotonicCubicFritschCarlson.dat' lw 1
replot 'advection_BFECC_Linear.dat' lw 1
#replot 'advection_BFECC_CatmullRom.dat' lw 1
#replot 'advection_BFECC_MonotonicCubicFedkiw.dat' lw 1
#replot 'advection_BFECC_MonotonicCubicFritschCarlson.dat' lw 1
replot 'advection_LaxWendroffCDS.dat' lw 1
pause -1