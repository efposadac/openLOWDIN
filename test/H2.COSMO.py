#!/usr/bin/env python
import os
from colorstring import *

testName = "H2.COSMO"
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refTotalEnergy = -1.117447274199

# Run calculation

os.system ("lowdin2 -i " + inputName)

output = open ( outputName, "r") 
outputRead = output.readlines()

# Values

for line in outputRead:
	if "TOTAL ENERGY =" in line :
		totalEnergy = float(line.split()[3]) 

diffTotalEnergy = abs(refTotalEnergy - totalEnergy)

if (diffTotalEnergy <= 1E-10 ) :
	print testName + str_green(" ... OK" )
else :
	print testName + str_red(" ... NOT OK")

	print "\tDifference HF: " + str( diffTotalEnergy )

output.close()
