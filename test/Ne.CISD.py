#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from colorstring import *

if len(sys.argv)==2:
    lowdinbin = sys.argv[1]
else:
    lowdinbin = "lowdin2"

testName = "Ne.CISD"
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refTotalEnergy = -128.522668749513
refFCIEnergy = -128.678103255522

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
    if "STATE:   1 ENERGY =" in line:
        FCIEnergy = float(line.split()[4])

diffTotalEnergy = abs(refTotalEnergy - totalEnergy)
diffFCIEnergy = abs(refFCIEnergy - FCIEnergy)

if (diffTotalEnergy <= 1E-8 and diffFCIEnergy <= 1E-6):
    print(testName + str_green(" ... OK"))
else:
    print(testName + str_red(" ... NOT OK"))
    print("Difference HF: " + str(diffTotalEnergy))
    print("Difference FCI: " + str(diffFCIEnergy))
    sys.exit(1)

output.close()
