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
"HF energy" : [-75.895520937848,1E-8],
"HF dipole" : [1.06383820,1E-7],
"HF quadrupole xx" : [-7.21841981,1E-6],
"HF quadrupole yy" : [-4.02024678,1E-6],
"HF quadrupole zz" : [-6.20812552,1E-6],
"CI 1" : [-75.911063510564,1E-8],
"CI dipole" : [1.01277057,1E-7],
"CI quadrupole xx" : [-7.26705299,1E-6],
"CI quadrupole yy" : [-4.21742308,1E-6],
"CI quadrupole zz" : [-6.30621482,1E-6]
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
HF_dipo = True 
HF_quad = True 

# Values
for i in range(0,len(outputRead)):
    line = outputRead[i]
    if "TOTAL ENERGY =" in line:
        testValues["HF energy"] = float(line.split()[3])
    if "STATE:   1 ENERGY =" in line:
        testValues["CI 1"] = float(line.split()[4])

    if "DIPOLE: (A.U.)" in line and HF_dipo:
        for j in range (i,len(outputRead)) :
            linej = outputRead[j]
            if "Total Dipole:"  in linej:
                testValues["HF dipole"] = float(linej.split()[5])
                HF_dipo = False
                break

    if "DIPOLE: (A.U.)" in line and not HF_dipo:
        for j in range (i,len(outputRead)) :
            linej = outputRead[j]
            if "Total Dipole:"  in linej:
                testValues["CI dipole"] = float(linej.split()[5])
                break

    if "QUADRUPOLE NON-TRACELESS: (DEBYE ANGS)" in line and HF_quad:
        for j in range (i,len(outputRead)) :
            linej = outputRead[j]
            if "Total Quadrupole:"  in linej:
                testValues["HF quadrupole xx"] = float(linej.split()[2])
                testValues["HF quadrupole yy"] = float(linej.split()[3])
                testValues["HF quadrupole zz"] = float(linej.split()[4])
                HF_quad = False
                break

    if "QUADRUPOLE NON-TRACELESS: (DEBYE ANGS)" in line and not HF_quad:
        for j in range (i,len(outputRead)) :
            linej = outputRead[j]
            if "Total Quadrupole:"  in linej:
                testValues["CI quadrupole xx"] = float(linej.split()[2])
                testValues["CI quadrupole yy"] = float(linej.split()[3])
                testValues["CI quadrupole zz"] = float(linej.split()[4])
                break

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
