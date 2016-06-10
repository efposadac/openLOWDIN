#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from colorstring import *

testName = "PsH.P3"
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refTotalEnergy = -0.662179820165
refOrb1Positron_KT = -5.3910
refOrb1Positron_EP2 = -5.8791
refOrb1Positron_P3 = -6.1426
refOrb1Positron_EP3 = -6.1731
refOrb1Positron_OVGF_A = -6.3559
refOrb1Positron_OVGF_B = -6.3685
refOrb1Positron_OVGF_C = -6.3470
refOrb1Positron_RENP3 = -6.2039

# Run calculation

status = os.system("lowdin2 -i " + inputName)

if status:
    print(testName + str_red(" ... NOT OK"))
    sys.exit(1)

output = open(outputName, "r")
outputRead = output.readlines()

# Values

count = 0
for line in outputRead:
    if "TOTAL ENERGY =" in line:
        totalEnergy = float(line.split()[3])

    if "SUMMARY OF PROPAGATOR RESULTS FOR THE SPIN-ORBITAL: 1  OF SPECIES:POSITRON" in line:

        orb1Positron_KT = float(outputRead[count + 4].split()[1])
        orb1Positron_EP2 = float(outputRead[count + 5].split()[1])
        orb1Positron_P3 = float(outputRead[count + 6].split()[1])
        orb1Positron_EP3 = float(outputRead[count + 7].split()[1])
        orb1Positron_OVGF_A = float(outputRead[count + 8].split()[2])
        orb1Positron_OVGF_B = float(outputRead[count + 9].split()[2])
        orb1Positron_OVGF_C = float(outputRead[count + 10].split()[2])
        orb1Positron_RENP3 = float(outputRead[count + 11].split()[1])

    count = count + 1

diffTotalEnergy = abs(refTotalEnergy - totalEnergy)
diffOrb1Positron_KT = abs(refOrb1Positron_KT - orb1Positron_KT)
diffOrb1Positron_EP2 = abs(refOrb1Positron_EP2 - orb1Positron_EP2)
diffOrb1Positron_P3 = abs(refOrb1Positron_P3 - orb1Positron_P3)
diffOrb1Positron_EP3 = abs(refOrb1Positron_EP3 - orb1Positron_EP3)
diffOrb1Positron_OVGF_A = abs(refOrb1Positron_OVGF_A - orb1Positron_OVGF_A)
diffOrb1Positron_OVGF_B = abs(refOrb1Positron_OVGF_B - orb1Positron_OVGF_B)
diffOrb1Positron_OVGF_C = abs(refOrb1Positron_OVGF_C - orb1Positron_OVGF_C)
diffOrb1Positron_RENP3 = abs(refOrb1Positron_RENP3 - orb1Positron_RENP3)

errorInP3 = 0.0001

if (diffTotalEnergy <= 1E-10 and
        diffOrb1Positron_KT <= errorInP3 and
        diffOrb1Positron_EP2 <= errorInP3 and
        diffOrb1Positron_P3 <= errorInP3 and
        diffOrb1Positron_EP3 <= errorInP3 and
        diffOrb1Positron_OVGF_A <= errorInP3 and
        diffOrb1Positron_OVGF_B <= errorInP3 and
        diffOrb1Positron_OVGF_C <= errorInP3 and
        diffOrb1Positron_RENP3 <= errorInP3):

    print(testName + str_green(" ... OK"))
else:
    print(testName + str_red(" ... NOT OK"))

    print("\tDifference HF: " + str(diffTotalEnergy))
    print("\tDifference in e+ values")
    print("\tDifference KT     " + str(diffOrb1Positron_KT))
    print("\tDifference EP2    " + str(diffOrb1Positron_EP2))
    print("\tDifference P3     " + str(diffOrb1Positron_P3))
    print("\tDifference EP3    " + str(diffOrb1Positron_EP3))
    print("\tDifference OVGF_A " + str(diffOrb1Positron_OVGF_A))
    print("\tDifference OVGF_B " + str(diffOrb1Positron_OVGF_B))
    print("\tDifference OVGF_C " + str(diffOrb1Positron_OVGF_C))
    print("\tDifference RENP3  " + str(diffOrb1Positron_RENP3))
    sys.exit(1)

output.close()
