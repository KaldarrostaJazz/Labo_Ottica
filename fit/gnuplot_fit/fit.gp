I = 70
k = 2 * pi / (632.8e-9)
D = 0.00030
d = 0.00010
c = 0.0602872
L = 1.4325
O = 1
f(x) = I * (cos(0.5*k*D*sin(atan((x-c)/L))))**2 * (sin(0.5*k*d*sin(atan((x-c)/L)))/(0.5*k*d*sin(atan((x-c)/L))))**2 + O
set fit logfile 'fit.log'
fit f(x) '../../osservazioni/lab3/12_130.dat' u 1:2:($3+0.3) yerrors via I, D, d, c, L, O

set multiplot title "Diffrazione da doppia fendìtura"
#set colorsequence default
set grid
set key
set xrange [0.050:0.070]
unset title

set origin 0, 0.34
set size 0.96, 0.64
unset xlabel
set ylabel "Intensità (ua)"
set format x ""
plot '../../osservazioni/lab3/12_130.dat' u 1:2:3 w yerrorbars title "Dati", [0.050:0.070] f(x) w l title "Fit"

set origin 0.045, 0
set size 0.915, 0.38
unset format x
unset ylabel
set xlabel "Posizione (m)"
#set ytics ("-0.2" -0.2, "-0.1" -0.1, "0" 0, "0.1" 0.1, "0.2", 0.2)
#set yrange [-0.25:0.25]
plot '../../osservazioni/lab3/12_130.dat' u 1:(($2-f($1))/$3) w l linetype 2 title "Residui"
