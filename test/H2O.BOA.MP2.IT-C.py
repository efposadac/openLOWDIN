#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from colorstring import *

testName = "H2O.BOA.MP2.IT-C"
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refTotalEnergy = -76.008843007735
refMP2Energy = -7.61654101128488605E+01

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
    if "E(MP2)=" in line:
        MP2Energy = float(line.split()[1])

diffTotalEnergy = abs(refTotalEnergy - totalEnergy)
diffMP2Energy = abs(refMP2Energy - MP2Energy)

if (diffTotalEnergy <= 1E-10 and diffMP2Energy <= 1E-10):
    print(testName + str_green(" ... OK"))
else:
    print(testName + str_red(" ... NOT OK"))
    print("Difference HF: " + str(diffTotalEnergy))
    print("Difference MP2: " + str(diffMP2Energy))
    sys.exit(1)


output.close()
