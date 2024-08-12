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
"HF energy" : [-112.737336707933,1E-8],
"CISD energy" : [-112.957030012578,1E-6],
"HF z-dipole" : [-0.33130728,1E-3],
"CISD z-dipole" :  [0.13640258,1E-3],
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
dipoleflag=False
ciflag=False
hfflag=True
for i in range(0,len(outputRead)):
    line = outputRead[i]
    if "TOTAL ENERGY =" in line:
        testValues["HF energy"] = float(line.split()[3])
    if "STATE:   1 ENERGY =" in line:
        testValues["CISD energy"] = float(line.split()[4])
    if "DIPOLE: (DEBYE)" in line:
        dipoleflag=True
    if "Total Dipole:" in line and dipoleflag and hfflag:
        testValues["HF z-dipole"] = float(line.split()[4])
        dipoleflag=False
        hfflag=False
    if "We are calculating properties for E-ALPHA in the CI ground state" in line:
        ciflag=True
    if "Total Dipole:" in line and dipoleflag and ciflag:
        testValues["CISD z-dipole"] = float(line.split()[4])
        dipoleflag=False

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

