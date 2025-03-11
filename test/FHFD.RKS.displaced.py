
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
    "KS Total Energy" : [-200.802898216123,1E-6],
    "E- Exc.Corr. energy" : [-17.449677426723,1E-3],
    "E-/H_1 Corr. energy" : [-0.040554548788,1E-4],
    "E-/H_2 Corr. energy" : [-0.035414098071,1E-4],
    "Number of E-" : [20.00000726,1E-4],
    "Number of H_1" : [0.99999439,1E-4],
    "Number of H_2" : [0.99998816,1E-4]
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
        testValues["KS Total Energy"] = float(line.split()[3])
    if "E- Exc.Corr. energy" in line:
        testValues["E- Exc.Corr. energy"] = float(line.split()[4])
    if "E-/H_1 Corr. energy" in line:
        testValues["E-/H_1 Corr. energy"] = float(line.split()[4])
    if "E-/H_2 Corr. energy" in line:
        testValues["E-/H_2 Corr. energy"] = float(line.split()[4])
    if "Number of E- particles in the final grid" in line:
        testValues["Number of E-"] = float(line.split()[8])
    if "Number of H_1 particles in the final grid" in line:
        testValues["Number of H_1"] = float(line.split()[8])
    if "Number of H_2 particles in the final grid" in line:
        testValues["Number of H_2"] = float(line.split()[8])
            
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
