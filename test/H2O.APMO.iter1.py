#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from colorstring import *

if len(sys.argv)==2:
    lowdinbin = sys.argv[1]
else:
    lowdinbin = "lowdin2"

testName = sys.argv[0][:-3]
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values and tolerance

refValues = {
"HF energy" : [-75.905888726266,1E-8],
"Iterations" : [1094,800],
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
    if " Total energy converged after" in line:
        testValues["Iterations"] = float(line.split()[4])

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

