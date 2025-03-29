#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"HF energy" : [-0.651377267180,1E-8],
"NOCI 1" : [-0.651856477835,1E-6],
"NOCI 2" : [-0.646106175835,1E-6],
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["CI 1"] = test.getCIEnergy(testName,1)
    testValues["CI 2"] = test.getCIEnergy(testName,2)
    testValues["NOCI 1"] = test.getCIEnergy(testName,1)
    testValues["NOCI 2"] = test.getCIEnergy(testName,2)
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
