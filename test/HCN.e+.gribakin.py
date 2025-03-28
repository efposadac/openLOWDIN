#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
    "HF energy" : [-92.9034152787,1E-7],
    "KT 1" : [-0.0016720654,1E-6]
}
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["KT 1"] = test.getHFeigenvalues(testName,"E+",1)
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
