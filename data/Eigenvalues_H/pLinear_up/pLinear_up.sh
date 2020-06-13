#!/bin/gnuplot
set terminal pdfcairo enhanced font "DejaVuSerif"
set output "b_1.000_l1_N_1500_E1.000e+007.pdf"
set xlabel '{/Symbol x}'
set ylabel '{/Symbol l_{1}}'
unset key
set mytics
set grid xtics ytics mytics
set grid back ls 6
set format y "%1.4e"
#set ytics (0.1,0.2,0.3,0.4,0.5,0.6,0.7)
#set xtics (0.1,1,10,100,1000,10000,100000,1000000)
set lmargin 15
plot 'b_1.000_eigen_1500_E1.000e+007.dat' u 1:2 w p lt 9 pt 13 ps 0.5  
