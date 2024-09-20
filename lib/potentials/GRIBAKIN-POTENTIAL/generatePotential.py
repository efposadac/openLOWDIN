#!/usr/bin/python
import math
import sys
import os
import subprocess

if len(sys.argv) > 1:
    potentialName = sys.argv[1]
else:
    print ("ERROR Please provide a name for the potential file, and run like 'python generatePotential.py POTNAME'")
    exit()

outFile = open (potentialName,"w")

# units
angstromToBohr = 1.88973

# -----------------------------
# Declaring the system

# species name
species="E+"
# molecule = atom, x, y, z, alpha_a, cutoff
#            [string, ang, ang, ang, a.u.**3, a.u.]
# benzene example
#molecule = """
#C  0.0000   1.4066  0.0000 10.32 2.20 
#C  1.2182   0.7033  0.0000 10.32 2.20 
#C  1.2182  -0.7033  0.0000 10.32 2.20 
#C  0.0000  -1.4066  0.0000 10.32 2.20 
#C -1.2182  -0.7033  0.0000 10.32 2.20 
#C -1.2182   0.7033  0.0000 10.32 2.20 
#H  0.0000   2.4998  0.0000 3.37 1.57 
#H  2.1649   1.2499  0.0000 3.37 1.57 
#H  2.1649  -1.2499  0.0000 3.37 1.57
#H  0.0000  -2.4998  0.0000 3.37 1.57 
#H -2.1649  -1.2499  0.0000 3.37 1.57 
#H -2.1649   1.2499  0.0000 3.37 1.57 
#"""
#molecule = """
#C  0.0000   1.4066  0.0000 8.268 2.30
#C  1.2182   0.7033  0.0000 8.268 2.30
#C  1.2182  -0.7033  0.0000 8.268 2.30
#C  0.0000  -1.4066  0.0000 8.268 2.30
#C -1.2182  -0.7033  0.0000 8.268 2.30
#C -1.2182   0.7033  0.0000 8.268 2.30
#H  0.0000   2.4998  0.0000 1.383 1.86
#H  2.1649   1.2499  0.0000 1.383 1.86
#H  2.1649  -1.2499  0.0000 1.383 1.86
#H  0.0000  -2.4998  0.0000 1.383 1.86
#H -2.1649  -1.2499  0.0000 1.383 1.86
#H -2.1649   1.2499  0.0000 1.383 1.86
#"""

# HCN example
#            [string, ang, ang, ang, ang**3, a.u.]
#molecule = """
#H   0.000   0.000   0.000 0.387 2.00
#C   0.000   0.000   1.059 1.283 2.00
#N   0.000   0.000   2.186 0.956 2.00
#"""

molecule = """
Be  0.000   0.000   0.000   38  2.686
"""

molecule = """
Be  0.000   0.000   -1.226801  38  2.686
Be  0.000   0.000   1.2268015  38  2.686
"""



# -----------------------------

# reshaping to an array format
molecule = molecule.split("\n")
molecule = molecule[1:] #removing the first jump (just for nice formatting) 
molecule = molecule[:-1] #removing the last jump 
n_atoms = len(molecule)
for i in range(0,n_atoms) :
    molecule[i] = molecule[i].split()

# convert units to a.u. (if needed)
for i in range(0,n_atoms):
    molecule[i][1] = float(molecule[i][1])*angstromToBohr
    molecule[i][2] = float(molecule[i][2])*angstromToBohr
    molecule[i][3] = float(molecule[i][3])*angstromToBohr
#    molecule[i][4] = float(molecule[i][4])*angstromToBohr**3
    molecule[i][4] = float(molecule[i][4])
    molecule[i][5] = float(molecule[i][5])
    #molecule[i][4] = molecule[i][4]*angstromToBohr**3

# -----------------------------
# parameters

exponents = (0.001,0.002,0.004,0.008,0.016,0.032,0.064,0.128,0.256,0.512,1,2,3,4,5,6,7,8,9,10,20,30,40,50,100)
maxR = 100 # a.u 
step = 0.04 # a.u 

# -----------------------------
# Temporal arrays
maxR = int(maxR/step)

angularMoment = list()
for i in range(0,len(exponents)) :
    angularMoment.append(0)

auxcoef = ""
for i in range(0,len(exponents)) :
    auxcoef = auxcoef + "c"+str(i)+", "
auxcoef = auxcoef[:-2]


# -----------------------------
# Potential headers
outFile.write("O-"+species+"\n")
outFile.write("#\n")
outFile.write("%i\n" % (len(exponents)*n_atoms))


# -----------------------------
# Potential fitting
ii = 0
print ("Fitting potential for: atom, x, y, z, alpha, cutoff")
print (" [string, ang, ang, ang, a.u.**3, a.u.]")

for atom in range(0,n_atoms):
    print ( "%s %.4f %.4f %.4f %.4f %.4f " % ( molecule[atom][0] , molecule[atom][1],  molecule[atom][2] , molecule[atom][3], molecule[atom][4], molecule[atom][5] ) ) 
    realPotentialFileName = potentialName + "."+ str(molecule[atom][0]) +".data"

    alpha = molecule[atom][4]
    # all potential will be centered on zero
    r = 0.0
    realPotentialFile = open (realPotentialFileName, "w")

    # calculating the potential in a line
    for i in xrange(1,maxR+1,1):
        i = i*step
        cutoff = molecule[atom][5]
        realPotentialFile.write(str(i) + " "+ str( -1.0*alpha/(2.0*(i-r)**4) * (1 - math.exp(-((i-r)**6)/(cutoff **6))) ) + "\n" )
    realPotentialFile.close()

    # saving in a string the fitting potential
    auxfunction = "g(x) = "
    for i in range(0,len(exponents)) :
        auxfunction = auxfunction + "c"+str(i)+"*exp(-"+str(exponents[i])+"*(x-"+str(r)+")**2) + "
    auxfunction = auxfunction[:-2]

    gnuplotFileName = potentialName + "."+ str(molecule[atom][0]) +".gnp"
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

    # calling gnuplot to fit the potential
    os.system("gnuplot "+gnuplotFileName + " 2> "+gnuplotFileName+".log")

    gnuplotOutput = open(gnuplotFileName+".log","r")
    gnuplotOut = gnuplotOutput.readlines()

    coeff = list()
    coefficients = list()

    # extract the fitting parameters from gnuplot log
    for j in range(0,len(gnuplotOut)):
        line = gnuplotOut[j]

        if "BEGIN RESULTS" in line:
            coeff = (gnuplotOut[j+1].split())

    gnuplotOutput.close()

    if len(coeff) == 0:
        print("Fitting error, results not found")
        continue

    # writting the fitted potential to Lowdin format
    coefficients = [float(k) for k in coeff]

    for i in range(0,len(exponents)) :
        outFile.write( str(i+ii+1) + " "+ str(angularMoment[i])+"\n" )
        outFile.write( "%.8f %.8e\n" % (exponents[i], coefficients[i]) )
        aux = ""
        for k in range(0,3):
            aux = aux + str(molecule[atom][k+1]) + " "
        outFile.write( aux + "\n")
    ii = ii + i + 1

outFile.close()
print ("The new fitted potential is stored in: %s " % potentialName)
