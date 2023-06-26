#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from colorstring import *

testName = sys.argv[0][:-3]
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values and tolerance

refValues = {
"HF energy" : [-0.651377267238,1E-8],
"HA-TIP Ext Pot" : [0.005339817047,1E-6],
"HB-TIP Ext Pot" : [0.005265789260,1E-6],
"HA-TIP/HB-TIP Hartree" : [0.003393498016,1E-6],
"HA-TIP/Fixed interact." : [-0.025239485153,1E-6],
"HB-TIP/Fixed interact." : [-0.010156190638,1E-6]
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
    if "HA-TIP Ext Pot" in line:
        testValues["HA-TIP Ext Pot"] = float(line.split()[5])
    if "HB-TIP Ext Pot" in line:
        testValues["HB-TIP Ext Pot"] = float(line.split()[5])
    if "HA-TIP/HB-TIP Hartree" in line:
        testValues["HA-TIP/HB-TIP Hartree"] = float(line.split()[4])
    if "HA-TIP/Fixed interact." in line:
        testValues["HA-TIP/Fixed interact."] = float(line.split()[4])
    if "HB-TIP/Fixed interact." in line:
        testValues["HB-TIP/Fixed interact."] = float(line.split()[4])
           

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
