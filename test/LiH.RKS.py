#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from colorstring import *

testName = "LiH.RKS"
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refExchangeCorrelationEnergy=-2.220930439475
refNumberOfE=3.99998704
refTotalEnergy=-8.070168253764
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


diffTotalEnergy = abs(refTotalEnergy - totalEnergy)
diffExchangeCorrelationEnergy = abs(refExchangeCorrelationEnergy - exchangeCorrelationEnergy)
diffRefNumberOfE = abs(refNumberOfE - numberOfE)

if (diffTotalEnergy <= 1E-6 and diffExchangeCorrelationEnergy <= 1E-3 and diffRefNumberOfE <= 1E-4 ):
    print(testName + str_green(" ... OK"))
else:
    print(testName + str_red(" ... NOT OK"))
    print("Difference RKS energy: " + str(diffTotalEnergy))
    print("Difference ExcCor energy: " + str(diffExchangeCorrelationEnergy))
    print("Difference E- number: " + str(diffRefNumberOfE))
    sys.exit(1)

output.close()
