#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from colorstring import *

testName = sys.argv[0][:-3]
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refValues = {
"HF energy" : [-76.008843007653,1E-8],
"Orb5alpha_P2" : [-10.9103,1E-3],
"Orb6alpha_P2" : [3.5274,1E-3],
"Orb5beta_P2" : [-10.9103,1E-3],
"Orb6beta_P2" : [3.5274,1E-3]
}

testValues = dict(refValues) #copy 
for value in testValues: #reset
    testValues[value] = 0 #reset
    
# Run calculation

status = os.system("lowdin2 -i " + inputName)

if status:
    print(testName + str_red(" ... NOT OK"))
    sys.exit(1)

output = open(outputName, "r")
outputRead = output.readlines()

# Values
for i in range(0,len(outputRead)):
    line = outputRead[i]
    if "TOTAL ENERGY =" in line:
        testValues["HF energy"] = float(line.split()[3])
    if "Results for spin-orbital: 5 of species: E-ALPHA" in line:
        for j in range(i,len(outputRead)): 
            if "Optimized second order pole:" in outputRead[j] :
                testValues["Orb5alpha_P2"] = float(outputRead[j].split()[4])
                break
    if "Results for spin-orbital: 6 of species: E-ALPHA" in line:
        for j in range(i,len(outputRead)): 
            if "Optimized second order pole:" in outputRead[j] :
                testValues["Orb6alpha_P2"] = float(outputRead[j].split()[4])
                break
    if "Results for spin-orbital: 5 of species: E-BETA" in line:
        for j in range(i,len(outputRead)): 
            if "Optimized second order pole:" in outputRead[j] :
                testValues["Orb5beta_P2"] = float(outputRead[j].split()[4])
                break
    if "Results for spin-orbital: 6 of species: E-BETA" in line:
        for j in range(i,len(outputRead)): 
            if "Optimized second order pole:" in outputRead[j] :
                testValues["Orb6beta_P2"] = float(outputRead[j].split()[4])
                break

passTest = True

for value in refValues:
    diffValue = abs(refValues[value][0] - testValues[value]) 
    if ( diffValue <= refValues[value][1] ):
        passTest = passTest * True
    else :
        passTest = passTest * False
        print("%s %.8f %.8f %.2e" % ( value, refValues[value][0], testValues[value], diffValue))

if passTest :
    print(testName + str_green(" ... OK"))
else:
    print(testName + str_red(" ... NOT OK"))
    sys.exit(1)

output.close()
