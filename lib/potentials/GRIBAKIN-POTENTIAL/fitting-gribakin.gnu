set terminal pdf transparent size 18.3 cm,14.6 cm lw 3 enhanced font "Nimbus,11"

outputName ="fitting.pdf"
set output outputName

inputdata = "data-gribakin4.13"

mean(x)= m
fit mean(x) inputdata u 1:2 via m
SST = FIT_WSSR/(FIT_NDF+1)
r = 4.13
#f(x) = -1*(1-exp(-1*(x**6)/(2.0)**6  ))/(2*x**4)
g(x) = c0*exp(-0.001*(x-r)**2) + c1*exp(-0.002*(x-r)**2) + c2*exp(-0.004*(x-r)**2) + c3*exp(-0.008*(x-r)**2) + c4*exp(-0.016*(x-r)**2) + c5*exp(-0.032*(x-r)**2) + c6*exp(-0.064*(x-r)**2) + c7*exp(-0.128*(x-r)**2) + c8*exp(-0.256*(x-r)**2) + c9*exp(-0.512*(x-r)**2) + c10*exp(-1*(x-r)**2) + c11*exp(-2*(x-r)**2) + c12*exp(-3*(x-r)**2) + c13*exp(-4*(x-r)**2) + c14*exp(-5*(x-r)**2) + c15*exp(-6*(x-r)**2) + c16*exp(-7*(x-r)**2) + c17*exp(-8*(x-r)**2) + c18*exp(-9*(x-r)**2) + c19*exp(-10*(x-r)**2) + c20*exp(-20*(x-r)**2) + c21*exp(-30*(x-r)**2) + c22*exp(-40*(x-r)**2) + c23*exp(-50*(x-r)**2) + c24*exp(-100*(x-r)**2) 
fit g(x) inputdata u 1:2 via c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17, c18, c19, c20, c21, c22, c23, c24
SSE=FIT_WSSR/(FIT_NDF)
SSR=SST-SSE
R2=SSR/SST
print R2
set  xrange[0:50]
#set  yrange[-0.03:0]
#plot inputdata u 1:2, g(x)
#plot inputdata u 1:2
plot g(x)

print c0,  c1,  c2,  c3,  c4,  c5,  c6,  c7,  c8,  c9,  c10,  c11,  c12,  c13,  c14,  c15,  c16,  c17,  c18,  c19,  c20,  c21,  c22,  c23,  c24


