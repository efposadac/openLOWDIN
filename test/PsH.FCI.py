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
"HF energy" : [-0.666783062050,1E-8],
"FCI 1" : [-0.743335966767,1E-8],
"FCI 2" : [-0.595807154865,1E-8],
"FCI 3" : [-0.595807154727,1E-8],
"Natural Occ 1 e+ 1" : [0.9189,1E-4],
"Natural Orb 1 e+ 1" : [0.005431,1E-4],
"Natural Orb 1 e+ 2" : [0.033634,1E-4],
"Natural Orb 1 e+ 3" : [0.039679,1E-4],
"Natural Orb 1 e+ 4" : [0.163628,1E-4],
"Natural Orb 1 e+ 5" : [0.617836,1E-4],
"Natural Orb 1 e+ 6" : [0.303937,1E-4]
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
    if "STATE:   1 ENERGY =" in line:
        testValues["FCI 1"] = float(line.split()[4])
    if "STATE:   2 ENERGY =" in line:
        testValues["FCI 2"] = float(line.split()[4])
    if "STATE:   3 ENERGY =" in line:
        testValues["FCI 3"] = float(line.split()[4])

    if "  Natural Orbitals in state:            1  for: POSITRON" in line:
        testValues["Natural Occ 1 e+ 1"] = float(outputRead[i+2].split()[0])
        for j in xrange(1,6+1):
            linej = outputRead[i+3+j]
            testValues["Natural Orb 1 e+ "+str(j)] = abs(float(linej.split()[3]))
            

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
