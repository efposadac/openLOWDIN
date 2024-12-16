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
"HF energy" : [-153.536301213424,1E-6],
"Embedded HF energy" : [-153.118946505217,1E-6],
"Embedded CISD" : [-153.310526831872,1E-6],
"Embedded MP2" : [-153.352174869549,1E-6],
"Embedded P2 H_1" : [-5.328000,1E-3]
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
flag=1
for i in range(0,len(outputRead)):
    line = outputRead[i]
    if "TOTAL ENERGY =" in line and flag==1:
        testValues["HF energy"] = float(line.split()[3])
        flag=2
        continue
    if "TOTAL ENERGY =" in line and flag==2:
        testValues["Embedded HF energy"] = float(line.split()[3])
        flag=3
        continue
    if "STATE:   1 ENERGY =" in line:
        testValues["Embedded CISD"] = float(line.split()[4])
    if "E(MP2) =" in line:
        testValues["Embedded MP2"] = float(line.split()[2])
    if "Optimized second order pole:" in line:
        testValues["Embedded P2 H_1"] = float(line.split()[4])

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
