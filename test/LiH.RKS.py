#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from colorstring import *

testName = "LiH.RKS"
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refExchangeCorrelationEnergy=-2.21707300
refNumberOfE=3.99991937
refTotalEnergy=-8.072145588372
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
    if "Number of E- particles in the final grid" in line:
        numberOfE = float(line.split()[8])


diffTotalEnergy = abs(refTotalEnergy - totalEnergy)
diffExchangeCorrelationEnergy = abs(refExchangeCorrelationEnergy - exchangeCorrelationEnergy)
diffRefNumberOfE = abs(refNumberOfE - numberOfE)

if (diffTotalEnergy <= 1E-5 and diffExchangeCorrelationEnergy <= 1E-5 and diffRefNumberOfE <= 1E-5 ):
    print(testName + str_green(" ... OK"))
else:
    print(testName + str_red(" ... NOT OK"))
    print("Difference RKS energy: " + str(diffTotalEnergy))
    print("Difference ExcCor energy: " + str(diffExchangeCorrelationEnergy))
    print("Difference E- number: " + str(diffRefNumberOfE))
    sys.exit(1)

output.close()
