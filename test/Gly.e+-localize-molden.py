#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"HF energy" : [-282.73698370,1E-6],
"e+HOMO molden" : [-0.02189320,1E-4],
"e-HOMO molden" : [-0.60297000,1E-4],
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["e+HOMO molden"] = test.getHOMOmolden(testName,"E+")
    testValues["e-HOMO molden"] = test.getHOMOmolden(testName,"E-")
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
