#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from colorstring import *

if len(sys.argv)==2:
    lowdinbin = sys.argv[1]
else:
    lowdinbin = "lowdin2"

testName = "HF.APMO"
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refTotalEnergy = -100.018384990309
refMP2Energy = -1.00318839332310560E+02

# Run calculation

status = os.system(lowdinbin + " -i " + inputName)

if status:
    print(testName + str_red(" ... NOT OK"))
    sys.exit(1)

output = open(outputName, "r")
outputRead = output.readlines()

# Values

for line in outputRead:
    if "TOTAL ENERGY =" in line:
        totalEnergy = float(line.split()[3])
    if "E(MP2) =" in line:
        MP2Energy = float(line.split()[2])

diffTotalEnergy = abs(refTotalEnergy - totalEnergy)
diffMP2Energy = abs(refMP2Energy - MP2Energy)

if (diffTotalEnergy <= 1E-7 and diffMP2Energy <= 1E-5):
    print(testName + str_green(" ... OK"))
else:
    print(testName + str_red(" ... NOT OK"))
    print("Difference HF: " + str(diffTotalEnergy))
    print("Difference MP2: " + str(diffMP2Energy))
    sys.exit(1)

output.close()
