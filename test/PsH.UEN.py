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
"HF energy" : [-6.660682661374E-01,1E-8],
"MP2 energy" : [-7.28307377384545984E-01,1E-8],
"NS-EN2" : [-7.55482317109419710E-01,1E-3],
"EN2" : [-7.55876101902173581E-01,1E-3]
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
    if "E(MP2) =" in line:
        testValues["MP2 energy"] = float(line.split()[2])
    if "E(NS-EN2) =" in line:
        testValues["NS-EN2"] = float(line.split()[2])
    if "E(EN2) =" in line:
        testValues["EN2"] = float(line.split()[2])

passTest = True

for value in refValues:
    diffValue = abs(refValues[value][0] - testValues[value]) 
    if ( diffValue <= refValues[value][1] ):
        passTest = passTest * True
        #print("%s %.8f %.8f %.2e" % ( value, refValues[value][0], testValues[value], diffValue))
    else :
        passTest = passTest * False
        print("%s %.8f %.8f %.2e" % ( value, refValues[value][0], testValues[value], diffValue))

if passTest :
    print(testName + str_green(" ... OK"))
else:
    print(testName + str_red(" ... NOT OK"))
    sys.exit(1)

output.close()
