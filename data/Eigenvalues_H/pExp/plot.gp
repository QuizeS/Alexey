array E[5]
array N[7]
array meow[8]

E[1] = 1
E[2] = 2
E[3] = 3
E[4] = 4
E[5] = 5

N[1] = 500
N[2] = 1000
N[3] = 1500
N[4] = 3000
N[5] = 5500
N[6] = 7000
N[7] = 10000

meow[1] = "λ_{1}"
meow[2] = "λ_{2}"
meow[3] = "λ_{3}"
meow[4] = "Δλ_{21}"
meow[5] = "Δλ_{31}"
meow[6] = "ф_{1}"
meow[7] = "φ_{2}"
meow[8] = "φ_{3}"

dfile(a,b) = sprintf("eigen_%d_E1.000e+00%d.dat", a, b)

# ### N,C :: NumberOfPoints: 1–7 ; Column: 2–8;
varE(b,c) = sprintf("plot for [i=1:5] dfile(N[%d], E[i]) u 1:%d t sprintf(\"%s, N=%s\", meow[%d], N[%d])", b, c, "%s", "%d", c-1, b)
# ### E,C :: Energy: 1–5 ; Column: 2–8;
varN(b,c) = sprintf("plot for [i=1:7] dfile(N[i], E[%d]) u 1:%d t sprintf(\"%s, E=10^{%s}\", meow[%d], E[%d])", b, c, "%s", "%d", c-1, b)
# plt_w(a,b,c) = "if(a == 0) { varE(b,c) } else { varN(b,c) }"
