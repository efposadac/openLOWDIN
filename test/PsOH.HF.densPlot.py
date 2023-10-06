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
plot1Name = testName + ".E-.3D.dens"                                                                      
plot2Name = testName + ".POSITRON.3D.dens"                                                                
# Reference values and tolerance

refValues = {
"HF energy" : [-75.587683637788,1E-8],
"Num e- in plot" : [10.45449875,1E-1],
"Num e+ in plot" : [0.99246436,1E-2],
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

plot1 = open(plot1Name, "r")
plot1Read = plot1.readlines()
sumE=0
for i in range(0,len(plot1Read)):
    line = plot1Read[i]
    if i > 1:
        values = line.split()
        if i == 3: x1=float(values[1])
        if i == 4: x2=float(values[1])
        if len(values) > 1: sumE+=(float(values[0])**2+float(values[1])**2)**(1.0/2.0)*float(values[2])
testValues["Num e- in plot"]=2.0*sumE*(x2-x1)**2
plot1.close()

plot2 = open(plot2Name, "r")
plot2Read = plot2.readlines()
sumP=0
for i in range(0,len(plot2Read)):
    line = plot2Read[i]
    if i > 1:
        values = line.split()
        if i == 3: x1=float(values[1])
        if i == 4: x2=float(values[1])
        if len(values) > 1: sumP+=(float(values[0])**2+float(values[1])**2)**(1.0/2.0)*float(values[2])
testValues["Num e+ in plot"]=2.0*sumP*(x2-x1)**2
plot2.close()

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

