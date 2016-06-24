#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from colorstring import *

testName = "H2O.APMO.Optmization"
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refdH1x = 0.000000 
refdH1y = 1.548552
refdH1z = 0.922960
refdHx = 1.548552
refdHy = 0.000000
refdHz = 0.949803
refdOx = 0.922960
refdOy = 0.949803
refdOz = 0.000000

# Run calculation

status = os.system("lowdin2 -i " + inputName)

if status:
    print(testName + str_red(" ... NOT OK"))
    sys.exit(1)

output = open(outputName, "r")
outputRead = output.readlines()

# Values

i = 0
for line in outputRead:
    i = i + 1
    if "END GEOMETRY OPTIMIZATION" in line:
        newline = outputRead[i + 21]
        dH1x = float(newline.split()[1])
        dH1y = float(newline.split()[2])
        dH1z = float(newline.split()[3])
        newline = outputRead[i + 22]
        dHx = float(newline.split()[1])
        dHy = float(newline.split()[2])
        dHz = float(newline.split()[3])
        newline = outputRead[i + 23]
        dOx = float(newline.split()[1])
        dOy = float(newline.split()[2])
        dOz = float(newline.split()[3])

diffHx = abs(refdH1x - dH1x)
diffHy = abs(refdH1y - dH1y)
diffHz = abs(refdH1z - dH1z)
diffDx = abs(refdHx  - dHx)
diffDy = abs(refdHy  - dHy)
diffDz = abs(refdHz  - dHz)
diffOx = abs(refdOx  - dOx)
diffOy = abs(refdOy  - dOy)
diffOz = abs(refdOz  - dOz)

if (diffHx <= 1E-6 and diffHy <= 1E-6 and diffHz <= 1E-6 and diffDx <= 1E-6 and diffDy <= 1E-6 and diffDz <= 1E-6 and diffOx <= 1E-6 and diffOy <= 1E-6 and diffOz <= 1E-6):
    print(testName + str_green(" ... OK"))
else:
    print(testName + str_red(" ... NOT OK"))
    print("Difference H1: " + str(diffHx) + str(diffHy) + str(diffHz))
    print("Difference H: " + str(diffDx) + str(diffDy) + str(diffDz))
    print("Difference O: " + str(diffOx) + str(diffOy) + str(diffOz))
    sys.exit(1)

output.close()
