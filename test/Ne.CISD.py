#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
    "HF Energy" : [-128.522668749513,1E-8],
    "FCI Energy" : [-128.678103255522,1E-6],
}
    return refValues

def getTestValues(testValues,testName):
    testValues["HF Energy"] = test.getSCFTotalEnergy(testName)
    testValues["FCI Energy"] = test.getCIEnergy(testName,1)
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
