#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from colorstring import *

testName = "H2O.APMO.UP2"
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refTotalEnergy = -75.930658814297
refOrb5alpha_P2 = -10.4052
refOrb5alphascs_P2 = -10.7531
refOrb5alphasos_P2 = -10.9292
refOrb1h1_P2 = -16.9153

# Run calculation

status = os.system("lowdin2 -i " + inputName)
if status:
    print(testName + str_red(" ... NOT OK"))
    sys.exit(1)

output = open(outputName, "r")
outputRead = output.readlines()

# Values

count = 0
Orb5 = False
Orb1H1 = False
for line in outputRead:
    if "TOTAL ENERGY =" in line:
        totalEnergy = float(line.split()[3])

    if "Results for spin-orbital: 5 of species: E-ALPHA" in line:
        Orb5 = True
    if "FactorOS: 1.00000 FactorSS: 1.00000" in line and Orb5 is True:
        Orb5alpha_P2 = float(outputRead[count + 1].split()[4])

    if "FactorOS: 1.20000 FactorSS: 0.33333" in line and Orb5 is True:
        Orb5alphascs_P2 = float(outputRead[count + 1].split()[4])

    if "FactorOS: 1.30000 FactorSS: 0.00000" in line and Orb5 is True:
        Orb5alphasos_P2 = float(outputRead[count + 1].split()[4])
        Orb5 = False

    if "Results for spin-orbital: 1 of species: H-A_1" in line:
        Orb1H1 = True

    if " Optimized second order pole:" in line and Orb1H1:
        Orb1h1_P2 = float(outputRead[count].split()[4])
        Orb1H1 = False

    count = count + 1

diffTotalEnergy = abs(refTotalEnergy - totalEnergy)
diffOrb5alpha_P2 = abs(refOrb5alpha_P2 - Orb5alpha_P2)
diffOrb5alphascs_P2 = abs(refOrb5alphascs_P2 - Orb5alphascs_P2)
diffOrb5alphasos_P2 = abs(refOrb5alphasos_P2 - Orb5alphasos_P2)
diffOrb1h1_P2 = abs(refOrb1h1_P2 - Orb1h1_P2)

if (diffTotalEnergy <= 1E-10 and diffOrb5alpha_P2 == 0 and diffOrb5alphascs_P2 == 0 and diffOrb5alphasos_P2 == 0 and
        diffOrb1h1_P2 == 0):
    print(testName + str_green(" ... OK"))
else:
    print(testName + str_red(" ... NOT OK"))
    print("Difference HF: " + str(diffTotalEnergy))
    print("Difference orbital 5 P2: " + str(diffOrb5alpha_P2))
    print("Difference orbital 5 SCS-P2: " + str(diffOrb5alphascs_P2))
    print("Difference orbital 5 SOS-P2: " + str(diffOrb5alphasos_P2))
    print("Difference orbital 1 H-A_1 P2: " + str(diffOrb1h1_P2))
    sys.exit(1)

output.close()
