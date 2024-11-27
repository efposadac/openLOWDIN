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
    "Total mass" :                  [40383.2869, 5E-4],
    "E- Kinetic energy" :           [36.558549208752,1E-6],
    "MUON Kinetic energy" :           [4.823223401186,1E-6],
    "C_12 Kinetic energy" :           [0.020362665915,1E-6],
    "HE_4 Kinetic energy" :           [0.042798597334,1E-6],
    "H_1 Kinetic energy" :            [0.018810346458,1E-6],
    "H_2 Kinetic energy" :            [0.013654378944,1E-6],
    "H_3 Kinetic energy" :            [0.011128664979,1E-6],
    "HF energy" :                     [-74.168958869716,1E-8]
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
    if "MASS (m_e)" in line:
        testValues["Total mass"] = float(line.split()[3])
    if "E- Kinetic energy" in line:
        testValues["E- Kinetic energy"] = float(line.split()[4])
    if "MUON Kinetic energy" in line:
        testValues["MUON Kinetic energy"] = float(line.split()[4])
    if "C_12 Kinetic energy" in line:
        testValues["C_12 Kinetic energy"] = float(line.split()[4])
    if "HE_4 Kinetic energy" in line:
        testValues["HE_4 Kinetic energy"] = float(line.split()[4])
    if "H_1 Kinetic energy" in line:
        testValues["H_1 Kinetic energy"] = float(line.split()[4])
    if "H_2 Kinetic energy" in line:
        testValues["H_2 Kinetic energy"] = float(line.split()[4])
    if "H_3 Kinetic energy" in line:
        testValues["H_3 Kinetic energy"] = float(line.split()[4])
        
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

