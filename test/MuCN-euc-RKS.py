#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"KS energy" : [-93.071088537483,1E-6],
"U+/E- Corr energy" : [-0.107102067712,1E-4]
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["KS energy"] = test.getSCFTotalEnergy(testName)
    testValues["U+/E- Corr energy"] = test.getDFTCorrEnergy(testName,"E-","U+")
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
