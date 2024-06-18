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
"HF energy" : [-343.383191892820,1E-8],
"U-HOMO" : [-371.890049816287,1E-1],
"H_1-HOMO" : [-1.019360160964,1E-4],
"He_4-HOMO" : [-652.366876763392,1E-1],
"e-HOMO" : [-0.585414450602,1E-4],
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
checkArray=[0,0,0,0]
for i in range(0,len(outputRead)):
    line = outputRead[i]
    if "TOTAL ENERGY =" in line:
        testValues["HF energy"] = float(line.split()[3])
    if "Eigenvalues for:" in line:
        species=line.split()[2]
        if species == "E-":
            checkArray[0]=1
        elif species == "H_1":
            checkArray[1]=1
        elif species == "U-":
            checkArray[2]=1
        elif species == "HE_4":
            checkArray[3]=1
        
    if "1 " in line and checkArray[0]==1:
        checkArray[0]=0
        testValues["e-HOMO"] = float(line.split()[1])
    if "1 " in line and checkArray[1]==1:
        checkArray[1]=0
        testValues["H_1-HOMO"] = float(line.split()[1])
    if "1 " in line and checkArray[2]==1:
        checkArray[2]=0
        testValues["U-HOMO"] = float(line.split()[1])
    if "1 " in line and checkArray[3]==1:
        checkArray[3]=0
        testValues["He_4-HOMO"] = float(line.split()[1])

        
output.close()

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

