#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
    "HF Energy" : [-7.916907112952,1E-7],
    "ExchangeCorrelationEnergy" : [-2.039879481164,1E-4],
    "NumberOfE" : [3.99998731,1E-4]
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
