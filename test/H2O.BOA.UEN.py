#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"HF energy" : [-7.605254225730E+01,1E-8],
"MP2 energy" : [-7.62994459305299983E+01,1E-5],
"NS-EN2" : [-7.63476164937737138E+01,1E-3],
"EN2" : [-7.63476682780967622E+01,1E-3]
}
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["MP2 energy"] = test.getMP2Energy(testName)
    testValues["NS-EN2"] = test.getNSEN2Energy(testName)
    testValues["EN2"] = test.getEN2Energy(testName)
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
