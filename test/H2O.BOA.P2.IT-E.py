#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"HF energy" : [-76.008843007735,1E-8],
"Orb5e_P2" : [-10.9105,1E-3],
"Orb6e_P2" : [3.5273,1E-3]
}
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["Orb5e_P2"] = test.getP2orbEnergy(testName,"E-",5)
    testValues["Orb6e_P2"] = test.getP2orbEnergy(testName,"E-",6)
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
