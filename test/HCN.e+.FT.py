#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from colorstring import *

testName = "HCN.e+.FT"
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refValues = {
"HF energy" : -92.901807583131,
"KT 1" : -0.0000643677
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
    if "BEGIN EIGENVALUES" in line:
        for j in range(i,len(outputRead)):
            linej = outputRead[j]
            if "POSITRON" in linej:
                for k in range(j+2,len(outputRead)): #j+2 is important
                    linek = outputRead[k]
                    if linek.split()[0] == "1" :
                        testValues["KT 1"] = float(linek.split()[1])
                        break
                break
        break

passTest = True

for value in refValues:
    diffValue = abs(refValues[value] - testValues[value]) 
    if ( diffValue <= 1E-8 ):
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
