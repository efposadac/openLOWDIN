#!/usr/bin/env python
import sys
import lowdinTestFunctions as test

if len(sys.argv)==2:
    lowdinbin = sys.argv[1]
else:
    lowdinbin = test.getDefaultLowdinBin()

testName = sys.argv[0][:-3]
inputName = testName + ".lowdin"
outputName = testName + ".out"

# Reference values and tolerance
refValues = {
    "HF energy" : [0.00,1E-8],
    "QDO zero energy" : [75.00,1E-8],
    "KT 2" : [50.0,1E-7],
    "KT 5" : [100.0,1E-7],
    "RMS R_c" : [0.002139,1E-4]        
}                       

testValues = dict(refValues) #copy 
for value in testValues: #reset
    testValues[value] = 0 #reset
    
# Run calculation
test.runCalculation(lowdinbin,inputName)

output = open(outputName, "r")
outputRead = output.readlines()

# Values
testValues["HF energy"] = test.getHFTotalEnergy(outputRead)
testValues["KT 2"] = test.getHFeigenvalues(outputRead,"QP+",2)
testValues["KT 5"] = test.getHFeigenvalues(outputRead,"QP+",5)
for i in range(0,len(outputRead)):
    line = outputRead[i]
    if "TOTAL ENERGY plus QDO ZPE =" in line:
        testValues["QDO zero energy"] = float(line.split()[6])
    if "EXPECTED RMS RADIUS OF QUANTUM SPECIES" in line:
        for j in range(i,len(outputRead)):
            linej = outputRead[j]
            if "QP+" in linej:
                testValues["RMS R_c"] = float(linej.split()[1])
                break
output.close()

passTest=test.verifyResults(refValues,testValues)

test.exitPrintResults(passTest,testName)

