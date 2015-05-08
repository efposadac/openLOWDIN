#!/usr/bin/env python
import os
from colorstring import *

testName = "H2O.APMO.UP2.IS"
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refTotalEnergy = -75.930658814297
refOrb1h1_P2 = -16.9153

# Run calculation

os.system ("lowdin2 -i " + inputName)

output = open ( outputName, "r") 
outputRead = output.readlines()

# Values

count = 0
Orb5 = False
Orb1H1 = False
for line in outputRead:
	if "TOTAL ENERGY =" in line :
		totalEnergy = float(line.split()[3]) 

	if "Results for spin-orbital: 1 of species: H-A_1" in line :
		Orb1H1 = True

	if " Optimized second order pole:" in line and Orb1H1 :
		Orb1h1_P2 = float(outputRead[count].split()[4]) 
		Orb1H1 = False

	count = count + 1

diffTotalEnergy = abs(refTotalEnergy - totalEnergy)
diffOrb1h1_P2 = abs(refOrb1h1_P2 - Orb1h1_P2 )

if (diffTotalEnergy <= 1E-12 and diffOrb1h1_P2 == 0 ) :
	print testName + str_green(" ... OK" )
else :
	print testName + str_red(" ... NOT OK")

	print "\tDifference HF: " + str( diffTotalEnergy )
	print "\tDifference orbital 1 H-A_1 P2: " + str( diffOrb1h1_P2 ) 	


output.close()
