#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"HF energy" : [-149.774990836259,1E-7],
"CI 1" : [-149.784025139570,1E-6],
"CI 2" : [-149.783893009711,1E-6],
"H_1 Kin 1" : [0.010857088555,1E-6],
"H_1 Kin 2" : [0.010941452494,1E-6],
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["CI 1"] = test.getCIEnergy(testName,1)
    testValues["CI 2"] = test.getCIEnergy(testName,2)
    testValues["H_1 Kin 1"] = test.getNOCIKineticEnergy(testName,"H_1",1)
    testValues["H_1 Kin 2"] = test.getNOCIKineticEnergy(testName,"H_1",2)
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
