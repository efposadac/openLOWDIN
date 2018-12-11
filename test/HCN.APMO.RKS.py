#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from colorstring import *

testName = "HCN.APMO.RKS"
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refExchangeCorrelationEnergy=-10.15758267
refNuclearElectronCorrelationEnergy=-0.03005334
refTotalEnergy=-93.367393799928
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
    if "Exchange correlation energy with the final grid" in line:
        exchangeCorrelationEnergy = float(line.split()[7])
    if "E-/H_1 Correlation contribution:" in line:
        nuclearElectronCorrelationEnergy = float(line.split()[3])


diffTotalEnergy = abs(refTotalEnergy - totalEnergy)
diffExchangeCorrelationEnergy = abs(refExchangeCorrelationEnergy - exchangeCorrelationEnergy)
diffNuclearElectronCorrelationEnergy = abs(refNuclearElectronCorrelationEnergy - nuclearElectronCorrelationEnergy)

if (diffTotalEnergy <= 1E-5 and diffExchangeCorrelationEnergy <= 1E-5 and diffNuclearElectronCorrelationEnergy <= 1E-5 ):
    print(testName + str_green(" ... OK"))
else:
    print(testName + str_red(" ... NOT OK"))
    print("Difference RKS energy: " + str(diffTotalEnergy))
    print("Difference ExcCor energy: " + str(diffExchangeCorrelationEnergy))
    print("Difference E- number: " + str(diffRefNumberOfE))
    sys.exit(1)

output.close()
