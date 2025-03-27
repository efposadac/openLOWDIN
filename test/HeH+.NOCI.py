#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"HF energy" : [-2.828229313778,1E-8],
"NOCI 1" : [-2.839688559383,1E-7],
"NOCI 2" : [-2.839448377440,1E-7],
"NOCI 4" : [-2.839448372563,1E-7]
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["CI 1"] = test.getCIEnergy(testName,1)
    testValues["CI 2"] = test.getCIEnergy(testName,2)
    testValues["NOCI 1"] = test.getCIEnergy(testName,1)
    testValues["NOCI 2"] = test.getCIEnergy(testName,2)
    testValues["NOCI 4"] = test.getCIEnergy(testName,4)
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
