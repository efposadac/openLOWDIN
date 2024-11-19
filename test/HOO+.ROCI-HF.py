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
"HF energy" : [-149.774990836259,1E-8],
"CI 1" : [-149.784025139570,1E-7],
"CI 2" : [-149.783893009711,1E-7],
"H_1 Kin 1" : [0.010857088555,1E-7],
"H_1 Kin 2" : [0.010941452494,1E-7],
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
stateFlag=0
for i in range(0,len(outputRead)):
    line = outputRead[i]
    if "TOTAL ENERGY =" in line:
        testValues["HF energy"] = float(line.split()[3])
    if "STATE:   1 ENERGY =" in line:
        testValues["CI 1"] = float(line.split()[4])
        stateFlag=1
    if "STATE:   2 ENERGY =" in line:
        testValues["CI 2"] = float(line.split()[4])
        stateFlag=2
    if "H_1 Kinetic energy =" in line:
        if stateFlag == 1:
            testValues["H_1 Kin 1"] = float(line.split()[4])
        elif stateFlag == 2:
            testValues["H_1 Kin 2"] = float(line.split()[4])
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
