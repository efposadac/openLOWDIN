#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"HF energy" : [-282.73698376,1E-6],
"Iterations" : [2,1],
"KT e+ 1" : [-0.02189320,1E-4],
"KT e- 20" : [-0.51095700,1E-4],
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["KT e+ 1"] = test.getHFeigenvalues(testName,"E+",1)
    testValues["KT e- 20"] = test.getHFeigenvalues(testName,"E-",20)
    testValues["Iterations"] = test.getSCFiterations(testName)
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
