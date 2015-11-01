#!/usr/bin/env python
import os
from colorstring import *

testName = "H2O.CCSD"
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refTotalEnergy = -75.984002559 # GAMESS
refCCSDEnergy = -76.120771451 #GAMESS

# Run calculation

os.system ("lowdin2 -i " + inputName)

output = open ( outputName, "r") 
outputRead = output.readlines()

# Values

for line in outputRead:
	if "TOTAL ENERGY =" in line :
		totalEnergy = float(line.split()[3]) 
	if "Total Energy (HF+CCSD)" in line :
		CCSDEnergy = float(line.split()[4]) 

diffTotalEnergy = abs(refTotalEnergy - totalEnergy)
diffCCSDEnergy = abs(refCCSDEnergy - CCSDEnergy)

if (diffTotalEnergy <= 1E-08 and CCSDEnergy <= 5E-04 ) : # 5E-04 ...
	print testName + str_green(" ... OK" )
else :
	print testName + str_red(" ... NOT OK")

	print "\tDifference HF: " + str( diffTotalEnergy )
	print "\tDifference CCSD: " + str( diffCCSDEnergy )


output.close()


