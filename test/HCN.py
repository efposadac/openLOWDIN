#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from colorstring import *

if len(sys.argv)==2:
    lowdinbin = sys.argv[1]
else:
    lowdinbin = "lowdin2"

testName = "HCN"
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refValues = {
"HF energy" : -92.901743215468
}

testValues = dict(refValues) #copy 
for value in testValues: #reset
    testValues[value] = 0 #reset
    
# Run calculation

status = os.system(lowdinbin + " -i " + inputName)

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
