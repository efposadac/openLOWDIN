#!/usr/bin/env python
import os
from colorstring import *

testName = "HF.APMO"
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refTotalEnergy = -100.018384990309
refMP2Energy = -1.00318839332310560E+02

# Run calculation

os.system ("lowdin2 -i " + inputName)

output = open ( outputName, "r") 
outputRead = output.readlines()

# Values

for line in outputRead:
	if "TOTAL ENERGY =" in line :
		totalEnergy = float(line.split()[3]) 
	if "E(MP2)=" in line :
		MP2Energy = float(line.split()[1]) 

diffTotalEnergy = abs(refTotalEnergy - totalEnergy)
diffMP2Energy = abs(refMP2Energy - MP2Energy)

if (diffTotalEnergy <= 1E-12 and diffMP2Energy <= 1E-12 ) :
	print testName + str_green(" ... OK" )
else :
	print testName + str_red(" ... NOT OK")

	print "\tDifference HF: " + str( diffTotalEnergy )
	print "\tDifference MP2: " + str( diffMP2Energy )

output.close()
