#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues={
    "HF energy" : [-99.639838220375,1E-8],
    "Orb1Positron_KT" : [-5.045419,1E-3],
    "Orb1Positron_EP2" : [-5.653732,1E-3],
    "Orb1Positron_P3" : [-5.915010,1E-3],
    "Orb1Positron_RENP3" : [-5.952235,1E-3]
    }
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    P3values = test.getP3results(testName,"POSITRON",1)
    testValues["Orb1Positron_KT"] = P3values[0]
    testValues["Orb1Positron_EP2"] = P3values[1]
    testValues["Orb1Positron_P3"] = P3values[2]
    testValues["Orb1Positron_RENP3"] = P3values[7]
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
