#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"Orb5alpha_P2" : [-10.768693,1E-3],
"Orb5alpha_SCS_P2" : [-10.740363,1E-3],
"Orb5alpha_SOS_P2" : [-10.725240,1E-3],
}
    return refValues

def getTestValues(testValues,testName):
    testValues["Orb5alpha_P2"] = test.getP2orbEnergy(testName,"E-ALPHA",5)
    testValues["Orb5alpha_SCS_P2"] = test.getP2orbScaledEnergy(testName,"E-ALPHA",5,"SCS")
    testValues["Orb5alpha_SOS_P2"] = test.getP2orbScaledEnergy(testName,"E-ALPHA",5,"SOS")
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
