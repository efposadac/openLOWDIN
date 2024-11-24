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
densplot1Name = testName + ".E-.3D.dens"
densplot2Name = testName + ".E+.3D.dens"
orbplot1Name = testName + ".E-.3D.orb3"
orbplot2Name = testName + ".E+.3D.orb1"
# Reference values and tolerance

refValues = {
"HF energy" : [-75.587683637788,1E-8],
"Num e- in densplot" : [10.0,1E-0],
"Num e- in orbplot" :  [2.0,5E-1],
"Num e+ in densplot" : [1.0,5E-2],
"Num e+ in orbplot" :  [1.0,5E-2],
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

output.close()

densplot1 = open(densplot1Name, "r")
densplot1Read = densplot1.readlines()
sumE=0
for i in range(0,len(densplot1Read)):
    line = densplot1Read[i]
    if i > 1:
        values = line.split()
        if i == 3: x1=float(values[1])
        if i == 4: x2=float(values[1])
        if len(values) > 1: sumE+=(float(values[0])**2+float(values[1])**2)**(1.0/2.0)*float(values[2])
testValues["Num e- in densplot"]=2.0*sumE*(x2-x1)**2
densplot1.close()

densplot2 = open(densplot2Name, "r")
densplot2Read = densplot2.readlines()
sumP=0
for i in range(0,len(densplot2Read)):
    line = densplot2Read[i]
    if i > 1:
        values = line.split()
        if i == 3: x1=float(values[1])
        if i == 4: x2=float(values[1])
        if len(values) > 1: sumP+=(float(values[0])**2+float(values[1])**2)**(1.0/2.0)*float(values[2])
testValues["Num e+ in densplot"]=2.0*sumP*(x2-x1)**2
densplot2.close()

orbplot1 = open(orbplot1Name, "r")
orbplot1Read = orbplot1.readlines()
sumE=0
for i in range(0,len(orbplot1Read)):
    line = orbplot1Read[i]
    if i > 1:
        values = line.split()
        if i == 3: x1=float(values[1])
        if i == 4: x2=float(values[1])
        if len(values) > 1: sumE+=(float(values[0])**2+float(values[1])**2)**(1.0/2.0)*float(values[2])**2
testValues["Num e- in orbplot"]=2.0*sumE*(x2-x1)**2
orbplot1.close()

orbplot2 = open(orbplot2Name, "r")
orbplot2Read = orbplot2.readlines()
sumP=0
for i in range(0,len(orbplot2Read)):
    line = orbplot2Read[i]
    if i > 1:
        values = line.split()
        if i == 3: x1=float(values[1])
        if i == 4: x2=float(values[1])
        if len(values) > 1: sumP+=(float(values[0])**2+float(values[1])**2)**(1.0/2.0)*float(values[2])**2
testValues["Num e+ in orbplot"]=2.0*sumP*(x2-x1)**2
orbplot2.close()

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

