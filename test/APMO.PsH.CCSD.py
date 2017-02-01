#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from colorstring import *

# PsH e-: SHARON-E-6S2P and e+: SHARON-E+6S2P

testName = "APMO.PsH.CCSD"
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refTotalEnergy = -0.6667830618
refmp2Energy = -0.7157914862
refccsdEnergy = -0.7468665505

omt= "OMP_NUM_THREADS=1"

# Run calculation

status = os.system("export " + omt + " && " + "lowdin2 -i " + inputName)

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
	mp2Energy = float(line.split()[1])
    if "Total CCSD energy:" in line:
	ccsdEnergy = float(line.split()[3])

diffTotalEnergy = abs(refTotalEnergy - totalEnergy)
diffmp2Energy = abs(refmp2Energy - mp2Energy)
diffccsdEnergy = abs(refccsdEnergy - ccsdEnergy)

if (diffTotalEnergy <= 1E-10):
    print(testName + str_green(" ... OK"))
else:
    print(testName + str_red(" ... NOT OK"))
    print("Difference HF: " + str(diffTotalEnergy))
    sys.exit(1)

if (diffmp2Energy <= 1E-10):
    print(testName + str_green(" ... OK"))
else:
    print(testName + str_red(" ... NOT OK"))
    print("Difference MP2: " + str(diffTotalEnergy))
    sys.exit(1)

if (diffccsdEnergy <= 1E-9):
    print(testName + str_green(" ... OK"))
else:
    print(testName + str_red(" ... NOT OK"))
    print("Difference CCSD: " + str(diffTotalEnergy))
    sys.exit(1)

output.close()
