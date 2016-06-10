#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from colorstring import *

testName = "H2O.BOA.P2"
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refTotalEnergy = -76.008843007735
refOrb5_P2 = -10.9105
refOrb6_P2 = 3.5273

# Run calculation

status = os.system("lowdin2 -i " + inputName)

if status:
    print(testName + str_red(" ... NOT OK"))
    sys.exit(1)

output = open(outputName, "r")
outputRead = output.readlines()

# Values

orbital5 = False
orbital6 = False

for line in outputRead:
    if "TOTAL ENERGY =" in line:
        totalEnergy = float(line.split()[3])
    if "Results for spin-orbital: 5 of species: E-" in line:
        orbital5 = True
    if "Optimized second order pole:" in line and orbital5 == True:
        Orb5_P2 = float(line.split()[4])
        orbital5 = False
    if "Results for spin-orbital: 6 of species: E-" in line:
        orbital6 = True
    if "Optimized second order pole:" in line and orbital6 == True:
        Orb6_P2 = float(line.split()[4])
        orbital6 = True

diffTotalEnergy = abs(refTotalEnergy - totalEnergy)
diffOrb5_P2 = abs(refOrb5_P2 - Orb5_P2)
diffOrb6_P2 = abs(refOrb6_P2 - Orb6_P2)

if (diffTotalEnergy <= 1E-10 and diffOrb5_P2 == 0 and diffOrb6_P2 == 0):
    print(testName + str_green(" ... OK"))
else:
    print(testName + str_red(" ... NOT OK"))
    print("Difference HF: " + str(diffTotalEnergy))
    print("Difference orbital 5 P2: " + str(diffOrb5_P2))
    print("Difference orbital 6 P2: " + str(diffOrb6_P2))
    sys.exit(1)

output.close()
