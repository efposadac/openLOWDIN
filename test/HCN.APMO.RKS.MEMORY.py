#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues={
    "E- Exc.Corr. energy":[-10.125588571062,1E-3],
    "E-/H_1 Corr. energy":[-0.033320942512,1E-3],
    "HF Energy":[-93.367985322582,1E-6]
}
    return refValues

def getTestValues(testValues,testName):
    testValues["HF Energy"] = test.getSCFTotalEnergy(testName)
    testValues["E- Exc.Corr. energy"] = test.getDFTExcCorrEnergy(testName,"E-")
    testValues["E-/H_1 Corr. energy"] = test.getDFTCorrEnergy(testName,"E-","H_1")
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
