#!/bin/gnuplot
set terminal pdfcairo enhanced font "DejaVuSerif"
set output "b_1.00_l3_N_1500_E1.000e+007.pdf"
set xlabel "ξ"
set ylabel "λ_{1}"
unset key
set mytics
set grid xtics ytics mytics
set format y "%1.2e"
#set ytics (0.1,0.2,0.3,0.4,0.5,0.6,0.7)
#set xtics (0.1,1,10,100,1000,10000,100000,1000000)
set lmargin 15
plot 'b_1.00_eigen_1500_E1.000e+007.dat' u 1:3 w p lt 9 pt 13 ps 0.5  
