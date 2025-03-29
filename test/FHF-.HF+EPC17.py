#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues={
    "Total Exc.Corr. energy":[-0.035742864487,1E-3],
    "Number of E-":[19.99997642,1E-4],
    "Number of H_1":[0.99999993,1E-4],
    "KS energy":[-199.524726631186,1E-6],
}
    return refValues

def getTestValues(testValues,testName):
    testValues["KS energy"] = test.getSCFTotalEnergy(testName)
    testValues["Number of E-"] = test.getParticlesInGrid(testName,"E-")
    testValues["Number of H_1"] = test.getParticlesInGrid(testName,"H_1")
    testValues["Total Exc.Corr. energy"] = test.getDFTTotalExcCorrEnergy(testName)
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
