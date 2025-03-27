#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"HF energy" : [-281.379249167035,1E-6],
"Iterations" : [2,1],
"KT e+ 1" : [-3.48432E-03,1E-4],
"KT e- 19" : [-6.32003E-01,1E-4],
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["KT e+ 1"] = test.getHFeigenvalues(testName,"POSITRON",1)
    testValues["KT e- 19"] = test.getHFeigenvalues(testName,"E-",19)
    testValues["Iterations"] = test.getSCFiterations(testName)
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
