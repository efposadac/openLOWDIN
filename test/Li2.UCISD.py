#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues={
    "HF energy" : [-14.781739628142,1E-7],
    "CISD energy" : [-14.821294812311,1E-6],
    "HF coefficient" : [0.926869142476,1E-4]
}
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["CISD energy"] = test.getCIEnergy(testName,1)
    testValues["HF coefficient"] = test.getHFCoefficient(testName)
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
