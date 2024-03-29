#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from colorstring import *

if len(sys.argv)==2:
    lowdinbin = sys.argv[1]
else:
    lowdinbin = "lowdin2"

testName = "HDO.GRADIENTS"
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values

refGradHx = 0.000000000000
refGradHy = 0.061073535862
refGradHz = -0.041117203100
refGradDx = 0.000000000000
refGradDy = -0.055359213575
refGradDz = -0.037187905984
refGradOx = 0.000000000000
refGradOy = -0.005714322288
refGradOz = 0.078305109084

diffTol=1E-5
# Run calculation

status = os.system(lowdinbin + " -i " + inputName)

if status:
    print(testName + str_red(" ... NOT OK"))
    sys.exit(1)

output = open(outputName, "r")
outputRead = output.readlines()

# Values

i = 0
for line in outputRead:
    i = i + 1
    if "dE/dx            dE/dy            dE/dz" in line:
        newline = outputRead[i + 1]
        gradHx = float(newline.split()[0])
        gradHy = float(newline.split()[1])
        gradHz = float(newline.split()[2])
        newline = outputRead[i + 2]
        gradDx = float(newline.split()[0])
        gradDy = float(newline.split()[1])
        gradDz = float(newline.split()[2])
        newline = outputRead[i + 3]
        gradOx = float(newline.split()[0])
        gradOy = float(newline.split()[1])
        gradOz = float(newline.split()[2])

diffHx = abs(refGradHx - gradHx)
diffHy = abs(refGradHy - gradHy)
diffHz = abs(refGradHz - gradHz)
diffDx = abs(refGradDx - gradDx)
diffDy = abs(refGradDy - gradDy)
diffDz = abs(refGradDz - gradDz)
diffOx = abs(refGradOx - gradOx)
diffOy = abs(refGradOy - gradOy)
diffOz = abs(refGradOz - gradOz)

if (diffHx <= diffTol and diffHy <= diffTol and diffHz <= diffTol and diffDx <= diffTol and diffDy <= diffTol and diffDz <= diffTol and diffOx <= diffTol and diffOy <= diffTol and diffOz <= diffTol):
    print(testName + str_green(" ... OK"))
else:
    print(testName + str_red(" ... NOT OK"))
    print("Difference H: " + str(diffHx)+ " " + str(diffHy)+ " " + str(diffHz))
    print("Difference D: " + str(diffDx)+ " " + str(diffDy)+ " " + str(diffDz))
    print("Difference O: " + str(diffOx)+ " " + str(diffOy)+ " " + str(diffOz))
    sys.exit(1)

output.close()
