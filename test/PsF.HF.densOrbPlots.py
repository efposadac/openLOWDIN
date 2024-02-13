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
densplot1Name = testName + ".E-.2D.dens"                                                                      
densplot2Name = testName + ".POSITRON.2D.dens"                                                                
orbplot1Name = testName + ".E-.2D.orb2"                                                                      
orbplot2Name = testName + ".POSITRON.2D.orb1"                                                                
# Reference values and tolerance

refValues = {
"HF energy" : [-99.635031860198,1E-8],
"Num e- in densplot" : [9.99830528,1E-3],
"Num e- in orbplot" :  [0.99956206,1E-3],
"Num e+ in densplot" : [0.99999922,1E-5],
"Num e+ in orbplot" :  [0.99999921,1E-3],
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
        if i == 3: x1=float(values[0])
        if i == 4: x2=float(values[0])
        if len(values) > 1: sumE+=float(values[0])**2*float(values[1])
testValues["Num e- in densplot"]=2.0*3.14159265359*sumE*(x2-x1)
densplot1.close()

densplot2 = open(densplot2Name, "r")
densplot2Read = densplot2.readlines()
sumP=0
for i in range(0,len(densplot2Read)):
    line = densplot2Read[i]
    if i > 1:
        values = line.split()
        if i == 3: x1=float(values[0])
        if i == 4: x2=float(values[0])
        if len(values) > 1: sumP+=float(values[0])**2*float(values[1])
testValues["Num e+ in densplot"]=2.0*3.14159265359*sumP*(x2-x1)
densplot2.close()

orbplot1 = open(orbplot1Name, "r")
orbplot1Read = orbplot1.readlines()
sumE=0
for i in range(0,len(orbplot1Read)):
    line = orbplot1Read[i]
    if i > 1:
        values = line.split()
        if i == 3: x1=float(values[0])
        if i == 4: x2=float(values[0])
        if len(values) > 1: sumE+=float(values[0])**2*float(values[1])**2
testValues["Num e- in orbplot"]=2.0*3.14159265359*sumE*(x2-x1)
orbplot1.close()

orbplot2 = open(orbplot2Name, "r")
orbplot2Read = orbplot2.readlines()
sumP=0
for i in range(0,len(orbplot2Read)):
    line = orbplot2Read[i]
    if i > 1:
        values = line.split()
        if i == 3: x1=float(values[0])
        if i == 4: x2=float(values[0])
        if len(values) > 1: sumP+=float(values[0])**2*float(values[1])**2
testValues["Num e+ in orbplot"]=2.0*3.14159265359*sumP*(x2-x1)
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

