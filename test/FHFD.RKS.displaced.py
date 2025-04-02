#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
    "KS energy" : [-200.802898216123,1E-6],
    "E- Exc.Corr. energy" : [-17.449677426723,1E-3],
    "E-/H_1 Corr. energy" : [-0.040554548788,1E-4],
    "E-/H_2 Corr. energy" : [-0.035414098071,1E-4],
    "Number of E-" : [20.00000726,1E-4],
    "Number of H_1" : [0.99999439,1E-4],
    "Number of H_2" : [0.99998816,1E-4],
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["KS energy"] = test.getSCFTotalEnergy(testName)
    testValues["E- Exc.Corr. energy"] = test.getDFTExcCorrEnergy(testName,"E-")
    testValues["E-/H_1 Corr. energy"] = test.getDFTCorrEnergy(testName,"E-","H_1")
    testValues["E-/H_2 Corr. energy"] = test.getDFTCorrEnergy(testName,"E-","H_2")
    testValues["Number of E-"] = test.getParticlesInGrid(testName,"E-")
    testValues["Number of H_1"] = test.getParticlesInGrid(testName,"H_1")
    testValues["Number of H_2"] = test.getParticlesInGrid(testName,"H_2")
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
