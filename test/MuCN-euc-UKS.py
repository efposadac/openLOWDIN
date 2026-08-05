#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
        "KS energy" : [-93.071088372363,1E-6],
        "U+/E-A Corr energy" : [-0.053547483452,1E-4],
        "U+/E-B Corr energy" : [-0.053547483452,1E-4]
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["KS energy"] = test.getSCFTotalEnergy(testName)
    testValues["U+/E-A Corr energy"] = test.getDFTCorrEnergy(testName,"E-ALPHA","U+")
    testValues["U+/E-B Corr energy"] = test.getDFTCorrEnergy(testName,"E-BETA","U+")
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
