#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from colorstring import *

testName = "Li2.UCISD"
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refTotalEnergy = -14.781739628142
refCISDEnergy = -14.821294812311
refHFCoefficient = 0.926869142476

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
    if "HF COEFFICIENT =" in line:
        HFCoefficient = float(line.split()[3])


diffTotalEnergy = abs(refTotalEnergy - totalEnergy)
diffCISDEnergy = abs(refCISDEnergy - CISDEnergy)
diffHFCoefficient = abs(refHFCoefficient - abs(HFCoefficient))

if (diffTotalEnergy <= 1E-10 and diffCISDEnergy <= 1E-8 and diffHFCoefficient <= 1E-5 ):
    print(testName + str_green(" ... OK"))
else:
    print(testName + str_red(" ... NOT OK"))
    print("Difference HF: " + str(diffTotalEnergy))
    print("Difference CISD: " + str(diffCISDEnergy))
    print("Difference HF Coefficient: " + str(diffHFCoefficient))
#    sys.exit(1)

output.close()
