#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
    "HF Energy" : [-8.083046278111,1E-7],
    "ExchangeCorrelationEnergy" : [-2.239695236675,1E-4],
    "NumberOfE" : [3.99998722,1E-4]
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF Energy"] = test.getSCFTotalEnergy(testName)
    testValues["NumberOfE"] = test.getParticlesInGrid(testName,"E-")
    testValues["ExchangeCorrelationEnergy"] = test.getDFTTotalExcCorrEnergy(testName)
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
