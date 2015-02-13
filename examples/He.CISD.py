#!/usr/bin/env python
import os
from colorstring import *

testName = "He.CISD"
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refTotalEnergy = -2.859895424516
refCISDEnergy = -2.876418360249

# Run calculation

os.system ("lowdin2 -i " + inputName)

output = open ( outputName, "r") 
outputRead = output.readlines()

# Values

for line in outputRead:
	if "TOTAL ENERGY =" in line :
		totalEnergy = float(line.split()[3]) 
	if "GROUND-STATE ENERGY =" in line :
		CISDEnergy = float(line.split()[3]) 

diffTotalEnergy = abs(refTotalEnergy - totalEnergy)
diffCISDEnergy = abs(refCISDEnergy - CISDEnergy)

if (diffTotalEnergy <= 1E-12 and CISDEnergy <=1E-12 ) :
	print testName + str_green(" ... OK" )
else :
	print testName + str_red(" ... NOT OK")

	print "\tDifference HF: " + str( diffTotalEnergy )
	print "\tDifference CISD: " + str( diffCISDEnergy )


output.close()


