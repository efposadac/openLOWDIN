#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from colorstring import *

testName = sys.argv[0][:-3]
inputName = testName + ".lowdin"
outputName = testName + ".out"
plot1Name = testName + ".E-.2D.dens"                                                                      
plot2Name = testName + ".POSITRON.2D.dens"                                                                
# Reference values and tolerance

refValues = {
"HF energy" : [-99.635031860198,1E-8],
"Num e- in plot" : [9.99830528,1E-3],
"Num e+ in plot" : [0.99999922,1E-5],
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

plot1 = open(plot1Name, "r")
plot1Read = plot1.readlines()
sumE=0
for i in range(0,len(plot1Read)):
    line = plot1Read[i]
    if i > 1:
        values = line.split()
        if i == 3: x1=float(values[0])
        if i == 4: x2=float(values[0])
        if len(values) > 1: sumE+=float(values[0])**2*float(values[1])
testValues["Num e- in plot"]=2.0*3.14159265359*sumE*(x2-x1)
plot1.close()

plot2 = open(plot2Name, "r")
plot2Read = plot2.readlines()
sumP=0
for i in range(0,len(plot2Read)):
    line = plot2Read[i]
    if i > 1:
        values = line.split()
        if i == 3: x1=float(values[0])
        if i == 4: x2=float(values[0])
        if len(values) > 1: sumP+=float(values[0])**2*float(values[1])
testValues["Num e+ in plot"]=2.0*3.14159265359*sumP*(x2-x1)
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

