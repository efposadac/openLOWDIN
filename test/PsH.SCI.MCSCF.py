#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"HF energy" : [-0.666783059430,1E-8],
"E_SCI+PT2" : [-0.743335409139,2E-5,2E-5],
-0.743333556464623
"E_MCSCF initial" : [-0.743334059515,2E-5],
}
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["E_SCI+PT2"] = test.getSCIPT2Energy(testName)
    testValues["E_MCSCF initial"] = test.getMCSCFInitialEnergy(testName)
    return

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
