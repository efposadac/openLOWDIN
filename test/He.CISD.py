#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from colorstring import *

testName = "He.CISD"
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refTotalEnergy = -2.859895424516
refCISDEnergy = -2.876418360249

# Run calculation

status = os.system("lowdin2 -i " + inputName)

if status:
    print(testName + str_red(" ... NOT OK"))
    sys.exit(1)

output = open(outputName, "r")
outputRead = output.readlines()

# Values

for line in outputRead:
    if "TOTAL ENERGY =" in line:
        totalEnergy = float(line.split()[3])
    if "STATE:   1 ENERGY =" in line:
        CISDEnergy = float(line.split()[4])

diffTotalEnergy = abs(refTotalEnergy - totalEnergy)
diffCISDEnergy = abs(refCISDEnergy - CISDEnergy)

if (diffTotalEnergy <= 1E-10 and CISDEnergy <= 1E-10):
    print(testName + str_green(" ... OK"))
else:
    print(testName + str_red(" ... NOT OK"))
    print("Difference HF: " + str(diffTotalEnergy))
    print("Difference CISD: " + str(diffCISDEnergy))
    sys.exit(1)

output.close()
