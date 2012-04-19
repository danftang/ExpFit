splot [][-11:0][] "computation.out" using 1:2:3
c(x,y) = A*x + B*y + C*x*y
fit [][-11:0][] c(x,y) "computation.out" using 1:2:3:(1) via A,B,C
replot c(x,y)

# interpolation using fastest FFT
f(x) = 34*(3*x + 2)/9;
replot f(x)

