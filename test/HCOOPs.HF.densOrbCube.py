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
cubes=[testName + ".E-.dens.cub",
       testName + ".E+.dens.cub",                                                                
       testName + ".E+.orb1.cub",
       testName + ".E-.orb12.cub"]
plots=[testName + ".E+.2D.orb1",
       testName + ".E+.3D.orb1"
    ]

# Reference values and tolerance
refValues = {
"HF energy" : [-188.362545831570,1E-8],
"Num E- in density cube" : [24.0,1E-1],
"Num E+ in density cube" : [1.0,1E-2],
"Num E+ in orbital cube" : [1.0,1E-2],
"Num E- in orbital cube" : [1.0,1E-2],
"Num E+ in 2D orbital plot" : [0.00010605,1E-4],
"Num E+ in 3D orbital plot" : [0.03267345,3E-2],
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

index=0
cube = open(cubes[index], "r")
cubeRead = cube.readlines()
sumPart=0
for i in range(0,len(cubeRead)):
    line = cubeRead[i]
    if i == 3: step=float(line.split()[1])
    if i > 10:
        values = line.split()
        for j in range(0,len(values)):
            sumPart+=float(values[j])
testValues["Num E- in density cube"]=sumPart*step**3
cube.close()

index=1
cube = open(cubes[index], "r")
cubeRead = cube.readlines()
sumPart=0
for i in range(0,len(cubeRead)):
    line = cubeRead[i]
    if i == 3: step=float(line.split()[1])
    if i > 10:
        values = line.split()
        for j in range(0,len(values)):
            sumPart+=float(values[j])
testValues["Num E+ in density cube"]=sumPart*step**3
cube.close()

index=2
cube = open(cubes[index], "r")
cubeRead = cube.readlines()
sumPart=0
for i in range(0,len(cubeRead)):
    line = cubeRead[i]
    if i == 3: step=float(line.split()[1])
    if i > 10:
        values = line.split()
        for j in range(0,len(values)):
            sumPart+=float(values[j])**2
testValues["Num E+ in orbital cube"]=sumPart*step**3
cube.close()

index=3
cube = open(cubes[index], "r")
cubeRead = cube.readlines()
sumPart=0
for i in range(0,len(cubeRead)):
    line = cubeRead[i]
    if i == 3: step=float(line.split()[1])
    if i > 10:
        values = line.split()
        for j in range(0,len(values)):
            sumPart+=float(values[j])**2
testValues["Num E- in orbital cube"]=sumPart*step**3
cube.close()

index=0
orbplotName=plots[index]
orbplot = open(orbplotName, "r")
orbplotRead = orbplot.readlines()
sumPart=0
for i in range(0,len(orbplotRead)):
    line = orbplotRead[i]
    if i > 1:
        values = line.split()
        if i == 3: x1=float(values[0])
        if i == 4: x2=float(values[0])
        if len(values) > 1 and float(values[0]) > 0.0 : sumPart+=float(values[1])**2
testValues["Num E+ in 2D orbital plot"]=sumPart*(x2-x1)**3
orbplot.close()

index=1
orbplotName=plots[index]
orbplot = open(orbplotName, "r")
orbplotRead = orbplot.readlines()
sumPart=0
for i in range(0,len(orbplotRead)):
    line = orbplotRead[i]
    if i > 1:
        values = line.split()
        if i == 3: x1=float(values[1])
        if i == 4: x2=float(values[1])
        if len(values) > 1: sumPart+=float(values[2])**2
testValues["Num E+ in 3D orbital plot"]=sumPart*(x2-x1)**3
orbplot.close()


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

