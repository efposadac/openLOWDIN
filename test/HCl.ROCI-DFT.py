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
"KS energy" : [-460.740603426038,1E-6],
"CI 1" : [-460.763938669911,1E-6],
"CI 2" : [-460.763767549387,1E-6],
"H_1 Kin 1" : [0.004031326872,1E-6],
"H_1 Kin 2" : [0.004117680624,1E-6],
"E-/H_1 Corr 1" : [-0.026900321722,1E-4],
"E-/H_1 Corr 2" : [-0.026900318425,1E-4],
"scaled CI 1" : [-460.753068611851,1E-6],
"scaled CI 2" : [-460.752986688385,1E-6],
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
        testValues["KS energy"] = float(line.split()[3])
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
    if "E-/H_1 DFTcorrelation energy =" in line:
        if stateFlag == 1:
            testValues["E-/H_1 Corr 1"] = float(line.split()[4])
        elif stateFlag == 2:
            testValues["E-/H_1 Corr 2"] = float(line.split()[4])
    if "STATE:   3 ENERGY =" in line:
            stateFlag = 3
    if "STATE:   1 SCALED ENERGY =" in line:
        testValues["scaled CI 1"] = float(line.split()[5])
        stateFlag=1
    if "STATE:   2 SCALED ENERGY =" in line:
        testValues["scaled CI 2"] = float(line.split()[5])
        stateFlag=2

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
