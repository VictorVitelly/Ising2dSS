set terminal qt size 1400,600
set grid x,y
set lmargin 15
set bmargin 6
set title 'Critical temperature' font ',22'
set xlabel 'q' font ',20'
set ylabel 'T_c' font ',20' offset -3.0 rotate by 0
set xtics font ',16'
set ytics font ',16'
set xrange [0.88:1.12]
set yrange [2:2.6]
set key font ',14'
set key outside

array point[1]

p 'results.dat' u 1:2:3 w errorbars lw 2 title 'Simulation', point us (1):(2/log(1+sqrt(2))) pt 7 lw 2 title 'Theoretical q=1'

pause -1
