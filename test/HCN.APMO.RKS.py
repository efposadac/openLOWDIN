#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from colorstring import *

if len(sys.argv)==2:
    lowdinbin = sys.argv[1]
else:
    lowdinbin = "lowdin2"

testName = "HCN.APMO.RKS"
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refExchangeCorrelationEnergy=-10.125588571062
refNuclearElectronCorrelationEnergy=-0.033320942512
refTotalEnergy=-93.367985322582
# Run calculation

status = os.system(lowdinbin + " -i " + inputName)

if status:
    print(testName + str_red(" ... NOT OK"))
    sys.exit(1)

output = open(outputName, "r")
outputRead = output.readlines()

# Values

for line in outputRead:
    if "TOTAL ENERGY =" in line:
        totalEnergy = float(line.split()[3])
    if "E- Exc.Corr. energy" in line:
        exchangeCorrelationEnergy = float(line.split()[4])
    if "E-/H_1 Corr. energy" in line:
        nuclearElectronCorrelationEnergy = float(line.split()[4])


diffTotalEnergy = abs(refTotalEnergy - totalEnergy)
diffExchangeCorrelationEnergy = abs(refExchangeCorrelationEnergy - exchangeCorrelationEnergy)
diffNuclearElectronCorrelationEnergy = abs(refNuclearElectronCorrelationEnergy - nuclearElectronCorrelationEnergy)

if (diffTotalEnergy <= 1E-6 and diffExchangeCorrelationEnergy <= 1E-3 and diffNuclearElectronCorrelationEnergy <= 1E-3 ):
    print(testName + str_green(" ... OK"))
else:
    print(testName + str_red(" ... NOT OK"))
    print("Difference RKS energy: " + str(diffTotalEnergy))
    print("Difference ExcCor energy: " + str(diffExchangeCorrelationEnergy))
    print("Difference NucElCorr energy: " + str(diffNuclearElectronCorrelationEnergy))
    sys.exit(1)

output.close()
