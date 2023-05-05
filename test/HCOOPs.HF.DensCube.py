#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from colorstring import *

testName = sys.argv[0][:-3]
inputName = testName + ".lowdin"
outputName = testName + ".out"
cube1Name = testName + ".E-.dens.cub"                                                                      
cube2Name = testName + ".POSITRON.dens.cub"                                                                
# Reference values and tolerance

refValues = {
"HF energy" : [-188.362545831570,1E-8],
"Num e- in cube" : [23.98744049,1E-1],
"Num e+ in cube" : [0.96711581,1E-2],
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

output.close()

cube1 = open(cube1Name, "r")
cube1Read = cube1.readlines()
sumE=0
for i in range(0,len(cube1Read)):
    line = cube1Read[i]
    if i == 3: step=float(line.split()[1])
    if i > 10:
        values = line.split()
        for j in range(0,len(values)):
            sumE+=float(values[j])
testValues["Num e- in cube"]=sumE*step**3
cube1.close()

cube2 = open(cube2Name, "r")
cube2Read = cube2.readlines()
sumP=0
for i in range(0,len(cube2Read)):
    line = cube2Read[i]
    if i == 3: step=float(line.split()[1])
    if i > 10:
        values = line.split()
        for j in range(0,len(values)):
            sumP+=float(values[j])
testValues["Num e+ in cube"]=sumP*step**3
cube2.close()

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

