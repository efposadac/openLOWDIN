#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from colorstring import *

testName = "HCN.e+.gribakin.uhf"
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refValues = {
    "HF energy" : [-92.903413350633,1E-7],
    "KT 1" : [-0.0016701915,1E-6]
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
    if "Eigenvalues for: POSITRON" in line:
        for j in range(i,len(outputRead)):
            linej = outputRead[j]
            if "1" in linej:
                testValues["KT 1"] = float(linej.split()[1])
                break

passTest = True

for value in refValues:
    diffValue = abs(refValues[value][0] - testValues[value]) 
    if ( diffValue <= refValues[value][1] ):
        passTest = passTest * True
    else :
        passTest = passTest * False
        print(value + " " + str(refValues[value]) +" " +  str(testValues[value]) + " "+ str(diffValue))

if passTest :
    print(testName + str_green(" ... OK"))
else:
    print(testName + str_red(" ... NOT OK"))
    sys.exit(1)

output.close()
