#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from colorstring import *

testName = "FHF-.HF+EPC17"
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refExchangeCorrelationEnergy=-0.035742864487
refNumberOfE=19.99997642
refNumberOfP=0.99999993
refTotalEnergy=-199.524726631186
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
    if "Total Exchange Correlation energy" in line:
        exchangeCorrelationEnergy = float(line.split()[5])
    if "Number of E- particles in the final grid" in line:
        numberOfE = float(line.split()[8])
    if "Number of H_1 particles in the final grid" in line:
        numberOfP = float(line.split()[8])


diffTotalEnergy = abs(refTotalEnergy - totalEnergy)
diffExchangeCorrelationEnergy = abs(refExchangeCorrelationEnergy - exchangeCorrelationEnergy)
diffRefNumberOfE = abs(refNumberOfE - numberOfE)
diffRefNumberOfP = abs(refNumberOfP - numberOfP)

if (diffTotalEnergy <= 1E-6 and diffExchangeCorrelationEnergy <= 1E-3 and diffRefNumberOfE <= 1E-4 and diffRefNumberOfP <= 1E-4 ):
    print(testName + str_green(" ... OK"))
else:
    print(testName + str_red(" ... NOT OK"))
    print("Difference RKS energy: " + str(diffTotalEnergy))
    print("Difference ExcCor energy: " + str(diffExchangeCorrelationEnergy))
    print("Difference E- number: " + str(diffRefNumberOfE))
    print("Difference H_1 number: " + str(diffRefNumberOfP))
    sys.exit(1)

output.close()
