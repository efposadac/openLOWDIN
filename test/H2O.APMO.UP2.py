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
"HF energy" : [-75.930658814297,1E-8],
"Orb5alpha_P2" : [-10.4052,1E-3],
"Orb5alpha_SCS_P2" : [-10.7531,1E-3],
"Orb5alpha_SOS_P2" : [-10.9292,1E-3],
"Orb6alpha_P2" : [3.4416,1E-3],
"Orb1H1a_P2" : [-16.9153,1E-3],
"Orb2H1a_P2" : [-50.7107,1E-2],
"Orb1H1b_P2" : [-16.9153,1E-3],
"Orb2H1b_P2" : [-50.7107,1E-2]
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
    if "Results for spin-orbital: 5 of species: E-ALPHA" in line:
        for j in range(i,len(outputRead)): 
            if "FactorOS: 1.20000 FactorSS: 0.33333" in outputRead[j] :
                testValues["Orb5alpha_SCS_P2"] = float(outputRead[j+1].split()[4])
                break
    if "Results for spin-orbital: 5 of species: E-ALPHA" in line:
        for j in range(i,len(outputRead)): 
            if "FactorOS: 1.30000 FactorSS: 0.00000" in outputRead[j] :
                testValues["Orb5alpha_SOS_P2"] = float(outputRead[j+1].split()[4])
                break
    if "Results for spin-orbital: 1 of species: H-A_1" in line:
        for j in range(i,len(outputRead)): 
            if "Optimized second order pole:" in outputRead[j] :
                testValues["Orb1H1a_P2"] = float(outputRead[j].split()[4])
                break
    if "Results for spin-orbital: 2 of species: H-A_1" in line:
        for j in range(i,len(outputRead)): 
            if "Optimized second order pole:" in outputRead[j] :
                testValues["Orb2H1a_P2"] = float(outputRead[j].split()[4])
                break
    if "Results for spin-orbital: 1 of species: H-B_1" in line:
        for j in range(i,len(outputRead)): 
            if "Optimized second order pole:" in outputRead[j] :
                testValues["Orb1H1b_P2"] = float(outputRead[j].split()[4])
                break
    if "Results for spin-orbital: 2 of species: H-B_1" in line:
        for j in range(i,len(outputRead)): 
            if "Optimized second order pole:" in outputRead[j] :
                testValues["Orb2H1b_P2"] = float(outputRead[j].split()[4])
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
