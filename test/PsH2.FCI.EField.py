#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
    "HF energy" : [-1.330593658435,1E-8],
    "CI 1" : [-1.472333332766,1E-8],
    "HF dipole" : [0.35846814,1E-7],
    "CI dipole" : [0.27622697,1E-4],
    "HF Ext Pot" : [-0.006405591748,1E-6],               
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["CI 1"] = test.getCIEnergy(testName,1)
    testValues["HF dipole"] = test.getSCFDipole(testName,"total","A.U.")
    testValues["CI dipole"] = test.getCIDipole(testName,"total","A.U.")
    testValues["HF Ext Pot"] = test.getSCFTotalExtPotEnergy(testName)
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
