#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"Num E+ in density cube" : [1.0,1E-2],
"Num E+ in density cube #2" : [1.0,1E-2],
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["Num E+ in density cube"] = test.getParticlesInDensCube(testName,"E+")
    testValues["Num E+ in density cube #2"] = test.getParticlesInDensCube(testName,"E+",2)
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
