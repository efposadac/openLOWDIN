#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
    "HF energy" : [-199.213740198536,1E-8],
    "eta e+" : [2.0,1E-8],
    "occupation e+" : [1.0,1E-8]
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["eta e+"] = test.getSpeciesEta(testName,"E+")
    testValues["occupation e+"] = test.getSpeciesOccupation(testName,"E+")
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
