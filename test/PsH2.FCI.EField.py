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
"HF energy" : [-1.330593658435,1E-8],
"CI 1" : [-1.472333332766,1E-8],
"HF dipole" : [0.35846814,1E-7],
"CI dipole" : [0.27622697,1E-4],
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
HF_prop = True 

# Values
for i in range(0,len(outputRead)):
    line = outputRead[i]
    if "TOTAL ENERGY =" in line:
        testValues["HF energy"] = float(line.split()[3])
    if "STATE:   1 ENERGY =" in line:
        testValues["CI 1"] = float(line.split()[4])

    if "DIPOLE: (A.U.)" in line and HF_prop:
        for j in range (i,len(outputRead)) :
            linej = outputRead[j]
            if "Total Dipole:"  in linej:
                testValues["HF dipole"] = float(linej.split()[5])
                HF_prop = False
                break

    if "DIPOLE: (A.U.)" in line and not HF_prop:
        for j in range (i,len(outputRead)) :
            linej = outputRead[j]
            if "Total Dipole:"  in linej:
                testValues["CI dipole"] = float(linej.split()[5])
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
