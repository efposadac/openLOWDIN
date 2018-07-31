#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from colorstring import *

testName = "LiH.TOP2"
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refTotalEnergy = -7.851847922443
refOrb2alpha_P2 = -7.5768

# Run calculation

status = os.system("lowdin2 -i " + inputName)
if status:
    print(testName + str_red(" ... NOT OK"))
    sys.exit(1)

output = open(outputName, "r")
outputRead = output.readlines()

# Values

count = 0
Orb2 = False
for line in outputRead:
    if "TOTAL ENERGY =" in line:
        totalEnergy = float(line.split()[3])

    if "Results for spin-orbital: 2 of species: E-ALPHA" in line:
        Orb2 = True
    if "FactorOS: 1.00000 FactorSS: 1.00000" in line and Orb2 is True:
        Orb2alpha_P2 = float(outputRead[count + 1].split()[4])

    count = count + 1

diffTotalEnergy = abs(refTotalEnergy - totalEnergy)
diffOrb2alpha_P2 = abs(refOrb2alpha_P2 - Orb2alpha_P2)

if (diffTotalEnergy <= 1E-10 and diffOrb2alpha_P2 == 0 ) :
    print(testName + str_green(" ... OK"))
else:
    print(testName + str_red(" ... NOT OK"))
    print("Difference HF: " + str(diffTotalEnergy))
    print("Difference orbital 2 P2: " + str(diffOrb2alpha_P2))
    sys.exit(1)

output.close()
