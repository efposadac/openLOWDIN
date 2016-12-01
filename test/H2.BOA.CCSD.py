#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from colorstring import *

# H2 6-31G

testName = "H2.BOA.CCSD"
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refTotalEnergy = -1.127274894913
refmp2Energy = -1.1448669169067
refccsdEnergy = -1.1520073923782

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

# H2O - STO-3G

testName = "H2O.BOA.CCSD"
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refTotalEnergy = -74.963681389320683
refmp2Energy = -74.999610505497799
refccsdEnergy = -75.013656085609256

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

if (diffmp2Energy <= 1E-9):
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
