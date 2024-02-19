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
molden1Name = testName + ".E-ALPHA.molden"                                                              
molden2Name = testName + ".E-BETA.molden"                                                             
molden3Name = testName + ".POSITRON.molden"                                                             
# Reference values and tolerance

refValues = {
"HF energy" : [-198.967220107085,1E-8],
"CISD+ energy" : [-198.977666606058,1E-6],
"Natural Occ 10 e-alpha 1" : [0.9987,1E-4],
"Natural Occ 11 e-alpha 1" : [0.0017,1E-4],
"Natural Occ 9 e-beta 1" : [0.9648,1E-4],
"Natural Occ 10 e-beta 1" : [0.0344,1E-4],
"Natural Occ 1 e+ 1" : [0.9653,1E-4],
"Natural Occ 2 e+ 1" : [0.0327,1E-4],
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
        testValues["CISD+ energy"] = float(line.split()[4])
output.close()

molden1 = open(molden1Name, "r")
molden1Read = molden1.readlines()
v=0
eigenv=[]
for i in range(0,len(molden1Read)):
    line = molden1Read[i]
    if "Occup=" in line:
        v+=1
        if v==10:
            testValues["Natural Occ 10 e-alpha 1"] = abs(float(line.split()[1]))
        if v==11:
            testValues["Natural Occ 11 e-alpha 1"] = abs(float(line.split()[1]))
        if "0.0" in line:
            break
molden1.close()

molden2 = open(molden2Name, "r")
molden2Read = molden2.readlines()
v=0
eigenv=[]
for i in range(0,len(molden2Read)):
    line = molden2Read[i]
    if "Occup=" in line:
        v+=1
        if v==9:
            testValues["Natural Occ 9 e-beta 1"] = abs(float(line.split()[1]))
        if v==10:
            testValues["Natural Occ 10 e-beta 1"] = abs(float(line.split()[1]))
        if "0.0" in line:
            break
molden2.close()

molden3 = open(molden3Name, "r")
molden3Read = molden3.readlines()
v=0
eigenv=[]
for i in range(0,len(molden3Read)):
    line = molden3Read[i]
    if "Occup=" in line:
        v+=1
        if v==1:
            testValues["Natural Occ 1 e+ 1"] = abs(float(line.split()[1]))
        if v==2:
            testValues["Natural Occ 2 e+ 1"] = abs(float(line.split()[1]))
        if "0.0" in line:
            break
molden3.close()

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
