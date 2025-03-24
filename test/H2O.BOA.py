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
"HF energy" : [-76.010288769789,1E-8],
"HF dipole" : [1.01124369,1E-7],
"HF quadrupole xx" : [-7.25729057,1E-5],
"HF quadrupole yy" : [-4.07831587,1E-5],
"HF quadrupole zz" : [-6.25445146,1E-5],
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
HF_dip = True 
HF_quad = True 

# Values
for i in range(0,len(outputRead)):
    line = outputRead[i]
    if "TOTAL ENERGY =" in line:
        testValues["HF energy"] = float(line.split()[3])

    if "DIPOLE: (A.U.)" in line and HF_dip:
        for j in range (i,len(outputRead)) :
            linej = outputRead[j]
            if "Total Dipole:"  in linej:
                testValues["HF dipole"] = float(linej.split()[5])
                HF_dip = False
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
