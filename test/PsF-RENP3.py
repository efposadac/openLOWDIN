#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from colorstring import *

if len(sys.argv)==2:
    lowdinbin = sys.argv[1]
else:
    lowdinbin = "lowdin2"
    
testName = "PsF-RENP3"
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refTotalEnergy = -99.639838220375
refOrb1Positron_KT = -5.045419
refOrb1Positron_EP2 = -5.653732
refOrb1Positron_P3 = -5.915010
refOrb1Positron_RENP3 = -5.952235

# Run calculation

status = os.system(lowdinbin+" -i " + inputName)

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
        orb1Positron_RENP3 = float(outputRead[count + 11].split()[1])

    count = count + 1

diffTotalEnergy = abs(refTotalEnergy - totalEnergy)
diffOrb1Positron_KT = abs(refOrb1Positron_KT - orb1Positron_KT)
diffOrb1Positron_EP2 = abs(refOrb1Positron_EP2 - orb1Positron_EP2)
diffOrb1Positron_P3 = abs(refOrb1Positron_P3 - orb1Positron_P3)
diffOrb1Positron_RENP3 = abs(refOrb1Positron_RENP3 - orb1Positron_RENP3)

errorInP3 = 0.001

if (diffTotalEnergy <= 1E-8 and
        diffOrb1Positron_KT <= errorInP3 and
        diffOrb1Positron_EP2 <= errorInP3 and
        diffOrb1Positron_P3 <= errorInP3 and
        diffOrb1Positron_RENP3 <= errorInP3):

    print(testName + str_green(" ... OK"))
else:
    print(testName + str_red(" ... NOT OK"))

    print("\tDifference HF: " + str(diffTotalEnergy))
    print("\tDifference in e+ values")
    print("\tDifference KT     " + str(diffOrb1Positron_KT))
    print("\tDifference EP2    " + str(diffOrb1Positron_EP2))
    print("\tDifference P3     " + str(diffOrb1Positron_P3))
    print("\tDifference RENP3  " + str(diffOrb1Positron_RENP3))
    sys.exit(1)

output.close()
