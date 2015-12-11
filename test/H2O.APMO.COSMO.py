#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from colorstring import *

testName = "H2O.APMO.COSMO"
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refTotalEnergy = -75.956252723159

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

diffTotalEnergy = abs(refTotalEnergy - totalEnergy)

if (diffTotalEnergy <= 1E-10):
    print(testName + str_green(" ... OK"))
else:
    print(testName + str_red(" ... NOT OK"))
    print("Difference HF: " + str(diffTotalEnergy))
    sys.exit(1)

output.close()
