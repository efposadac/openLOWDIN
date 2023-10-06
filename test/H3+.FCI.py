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
"HF energy" : [-1.255739154579,1E-8],
"FCI 1" : [-1.306431783235,1E-8],
"FCI 2" : [-1.291420173154,1E-8],
"FCI 3" : [-1.290426417940,1E-8],
"Natural Occ 1 H_1 1"  : [0.9919,1E-3],
"Natural Orb 1 H_1 1"  : [0.721928,1E-3],
"Natural Orb 1 H_1 2"  : [0.187901,1E-3],
"Natural Orb 1 H_1 3"  : [0.000000,1E-3],
"Natural Orb 1 H_1 4"  : [0.150257,1E-3],
"Natural Orb 1 H_1 5"  : [0.000000,1E-3],
"Natural Orb 1 H_1 6"  : [0.000000,1E-3],
"Natural Orb 1 H_1 7"  : [0.334924,1E-3],
"Natural Orb 1 H_1 8"  : [0.000000,1E-3],
"Natural Orb 1 H_1 9"  : [0.125434,1E-3],
"Natural Orb 1 H_1 10" : [0.000000,1E-3],
"Natural Orb 1 H_1 11" : [0.000000,1E-3],
"Natural Orb 1 H_1 12" : [0.116370,1E-3],
"Natural Orb 1 H_1 13" : [0.000000,1E-3],
"Natural Orb 1 H_1 14" : [0.117396,1E-3],
"Natural Orb 1 H_1 15" : [0.131109,1E-3],
"Natural Orb 1 H_1 16" : [0.000000,1E-3],
"Natural Orb 1 H_1 17" : [0.000000,1E-3],
"Natural Orb 1 H_1 18" : [0.127979,1E-3],
"Natural Orb 1 H_1 19" : [0.000000,1E-3],
"Natural Orb 1 H_1 20" : [0.161496,1E-3],
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
    if "STATE:   1 ENERGY =" in line:
        testValues["FCI 1"] = float(line.split()[4])
    if "STATE:   2 ENERGY =" in line:
        testValues["FCI 2"] = float(line.split()[4])
    if "STATE:   3 ENERGY =" in line:
        testValues["FCI 3"] = float(line.split()[4])

    if "  Natural Orbitals in state:            1  for: H_1" in line:
        testValues["Natural Occ 1 H_1 1"] = float(outputRead[i+2].split()[0])
        for j in range(1,20+1):
            linej = outputRead[i+3+j]
            testValues["Natural Orb 1 H_1 "+str(j)] = abs(float(linej.split()[3]))
            

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
