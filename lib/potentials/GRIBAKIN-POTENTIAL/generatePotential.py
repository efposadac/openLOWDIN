#!/usr/bin/python
import math
import sys
import os
import subprocess

potentialName = sys.argv[1]
outFile = open (potentialName,"w")

angstromToBohr = 1.88973

# "input"...
species="E+"
origin = ((0.000, 0.000, 0.000),(0.000, 0.000, 2.001),(0.000, 0.000, 4.130)) # a.u
polarizabilities =(0.387, 1.283, 0.956) # angstrom

# parameters

exponents = (0.001,0.002,0.004,0.008,0.016,0.032,0.064,0.128,0.256,0.512,1,2,3,4,5,6,7,8,9,10,20,30,40,50,100)
cutoff = 2.00 # a.u
maxR = 100 # a.u 
step = 0.04 # a.u 

#######

polarizabilities = [(angstromToBohr**3)*x for x in polarizabilities] # a.u.
maxR = int(maxR/step)
print maxR

angularMoment = list()
for i in range(0,len(exponents)) :
    angularMoment.append(0)

auxcoef = ""
for i in range(0,len(exponents)) :
    auxcoef = auxcoef + "c"+str(i)+", "
auxcoef = auxcoef[:-2]

######

outFile.write("O-"+species+"\n")
outFile.write("#\n")
outFile.write("%i\n" % (len(exponents)*len(origin)))


ii = 0
for atom in range(0,len(origin)):
    print origin[atom]
    realPotentialFileName = potentialName + "."+ str(origin[atom][2]) +".data"

    alpha = polarizabilities[atom]
    r = origin[atom][2] 
    realPotentialFile = open (realPotentialFileName, "w")

    for i in xrange(1,maxR+1,1):
        i = i*step
        realPotentialFile.write(str(i) + " "+ str( -1.0*alpha/(2.0*(i-r)**4) * (1 - math.exp(-((i-r)**6)/(cutoff**6))) ) + "\n" )
    realPotentialFile.close()

    auxfunction = "g(x) = "
    for i in range(0,len(exponents)) :
        auxfunction = auxfunction + "c"+str(i)+"*exp(-"+str(exponents[i])+"*(x-"+str(r)+")**2) + "
    auxfunction = auxfunction[:-2]
    #print auxfunction

    gnuplotFileName = potentialName + "."+ str(origin[atom][2]) +".gnp"
    gnuplotFile = open(gnuplotFileName,"w")
    gnuplotFile.write("""
set terminal pdf transparent size 18.3 cm,14.6 cm lw 3 enhanced font "Nimbus,11"

inputdata = '"""+realPotentialFileName+"""'

mean(x)= m
fit mean(x) inputdata u 1:2 via m
SST = FIT_WSSR/(FIT_NDF+1)
r = 4.13
"""+auxfunction+"""
fit g(x) inputdata u 1:2 via """+auxcoef+"""
SSE=FIT_WSSR/(FIT_NDF)
SSR=SST-SSE
R2=SSR/SST
print R2
set  xrange[0:50]
outputName = '"""+gnuplotFileName+'.pdf'+"""'
set output outputName


plot g(x)
print "BEGIN RESULTS"
print """+auxcoef+"""
""")
    gnuplotFile.close()

    os.system("gnuplot "+gnuplotFileName + " 2> "+gnuplotFileName+".log")

    gnuplotOutput = open(gnuplotFileName+".log","r")
    gnuplotOut = gnuplotOutput.readlines()

    coeff = list()
    coefficients = list()

    for j in range(0,len(gnuplotOut)):
        line = gnuplotOut[j]

        if "BEGIN RESULTS" in line:
            coeff = (gnuplotOut[j+1].split())

    gnuplotOutput.close()

    if len(coeff) == 0:
        print "Fitting error, results not found"
        continue

    coefficients = [float(k) for k in coeff]

    for i in range(0,len(exponents)) :
        outFile.write( str(i+ii+1) + " "+ str(angularMoment[i])+"\n" )
        outFile.write( "%.8f %.8e\n" % (exponents[i], coefficients[i]) )
        aux = ""
        for k in range(0,3):
            aux = aux + str(origin[atom][k]) + " "
        outFile.write( aux + "\n")
    ii = ii + i + 1

outFile.close()






