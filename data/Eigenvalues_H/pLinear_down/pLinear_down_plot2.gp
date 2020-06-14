set terminal pdfcairo enhanced font "DejaVuSerif"
set pointsize 0.1
array E[7]
array N[7]
array meow[20]

E[1] = 1
E[2] = 2
E[3] = 3
E[4] = 4
E[5] = 5
E[6] = 6
E[7] = 7

N[1] = 500
N[2] = 1000
N[3] = 1500
N[4] = 3000
N[5] = 5500
N[6] = 7000
N[7] = 10000

meow[1]  = "ξ"
meow[2]  = "λ_{1}"
meow[3]  = "λ_{2}"
meow[4]  = "λ_{3}"
meow[5]  = "Δλ_{21}"
meow[6]  = "Δλ_{31}"
meow[7]  = "φ_{1}"
meow[8]  = "φ_{2}"
meow[9]  = "φ_{3}"
meow[10] = "Re(ψ_1)"
meow[11] = "Re(ψ_2)"
meow[12] = "Re(ψ_3)"
meow[13] = "Im(ψ_1)"
meow[14] = "Im(ψ_2)"
meow[15] = "Im(ψ_3)"
meow[16] = "Δλ_{21}/ λ_{1}"
meow[17] = "Δλ_{31}/ λ_{1}"
meow[18] = "ψ_1"
meow[19] = "ψ_2"
meow[20] = "ψ_3"

dfile(f,a,b) = sprintf("b_%.2f_eigen_%d_E1.000e+00%d.dat", f, a, b)
ofile(f,a,b) = sprintf("b_%.3f_eigen_%d_%d.pdf", f, a, b)
Efile(f,a,b) = sprintf("b+%.3f_eigen_%d_%d.pdf", f, a , b)


# #############################################################################
# ### Column: see 'meow' array.
# ### L,N,C :: UpperLimit: 0.2/1.0 ; NumberOfPoints: 1–7 ; Column: 2–17 ; energy dependence;
varE(a,b,c) = sprintf("set xlabel \"%s\" ; set ylabel \"%s\"; plot for [i=1:7] dfile(%.3f, N[%d], E[i]) u 1:%d t sprintf(\"10^{%s}\", E[i])", "ξ", meow[c], a, b, c, "%d")
# ### L,E,C :: UpperLimit; 0.2/1.0 ; Energy: 1–7 ; Column: 2–17 ; umber of points dependence;
varN(a,b,c) = sprintf("set xlabel \"%s\" ; set ylabel \"%s\";  plot for [i=1:7] dfile(%.3f, N[i], E[%d]) u 1:%d t sprintf(\"%s\", N[i])", "ξ", meow[c], a, b, c, "%d")
# #############################################################################
# ### 3D plot, phase picture.
# ### L,N,C :: UpperLimit: 0.2/1.0 ; NumberOfPoints: 1–7 ; Component: 1–3 ; phase picture of psi_{1,2,3} ,energy dependence;
phPvE(a,b,c) = sprintf("set xlabel \"%s\" ; set ylabel \"%s\" ; set zlabel \"%s\"; splot for [i=1:7] dfile(%.3f, N[%d], E[i]) u %d:%d:1 t sprintf(\"10^{%s}\", E[i])", meow[c+9], meow[c+12], "ξ", a, b, c+9, c+12, "%d")
# ### L,E,C :: UpperLimit: 0.2/1.0 ; Energy: 1–7 ; Component: 1–3 ; phase picture of psi_{1,2,3}, number of points dependence; 
phPvN(a,b,c) = sprintf("set xlabel \"%s\" ; set ylabel \"%s\" ; set zlabel \"%s\"; splot for [i=1:7] dfile(%.3f, N[i], E[%d]) u %d:%d:1 t sprintf(\"%s\", N[i])", meow[c+9], meow[c+12], "ξ", a, b, c+9, c+12, "%d")
# #############################################################################
# ### Phase Picture, single energy, single points set.
# ### L,N,C,E :: UpperLimit: 0.2/1.0 ; NumberOfPoints: 1–7 ; Component: 1–3 ; Energy (power): 1–7 ;
phP(a,b,c,p) = sprintf("set xlabel \"%s\" ; set ylabel \"%s\" ; set zlabel \"%s\" ; set title sprintf(\"Phase picture for %s at E = 10^{%s} and N = %s\", E[%d], N[%d]) ; splot dfile(%.3f, N[%d], E[%d]) u %d:%d:1 t ''", meow[c+9], meow[c+12], "ξ", meow[c+17], "%d", "%d", p, b, a, b, p, c+9, c+12)

set output Efile(1.0,2,9)
eval varE(1.0,2,9) 
