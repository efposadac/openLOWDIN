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
molden1Name = testName + ".E-.molden"                                                                
molden2Name = testName + ".POSITRON.molden"                                                                      
# Reference values and tolerance

refValues = {
"HF energy" : [-281.379249259318,1E-8],
"e+HOMO" : [-3.48415E-03,1E-4],
"e-HOMO" : [-5.40431E-01,1E-4],
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

molden1 = open(molden1Name, "r")
molden1Read = molden1.readlines()
v=0
eigenv=[]
for i in range(0,len(molden1Read)):
    line = molden1Read[i]
    if "Ene=" in line:
        eigenv.append(float(line.split()[1]))
        v+=1
    if "Occup=" in line and "0.0" in line :
        testValues["e-HOMO"] = eigenv[v-2]
        break
molden1.close()

molden2 = open(molden2Name, "r")
molden2Read = molden2.readlines()
v=0
eigenv=[]
for i in range(0,len(molden2Read)):
    line = molden2Read[i]
    if "Ene=" in line:
        eigenv.append(float(line.split()[1]))
        v+=1
    if "Occup=" in line and "0.0" in line :
        testValues["e+HOMO"] = eigenv[v-2] 
        break
molden2.close()

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

