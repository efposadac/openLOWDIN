#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues={
    "HF energy" : [-0.662179820165,1E-8],
    "Orb1Positron_KT" : [-5.3910,1E-3],
    "Orb1Positron_EP2" : [-5.8791,1E-3],
    "Orb1Positron_P3" : [-6.1426,1E-3],
    "Orb1Positron_EP3" : [-6.1731,1E-3],
    "Orb1Positron_OVGF_A" : [-6.3559,1E-3],
    "Orb1Positron_OVGF_B" : [-6.3685,1E-3],
    "Orb1Positron_OVGF_C" : [-6.3470,1E-3],
    "Orb1Positron_RENP3" : [-6.2039,1E-3]
}
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    P3values = test.getP3results(testName,"POSITRON",1)
    testValues["Orb1Positron_KT"] = P3values[0]
    testValues["Orb1Positron_EP2"] = P3values[1]
    testValues["Orb1Positron_P3"] = P3values[2]
    testValues["Orb1Positron_EP3"] = P3values[3]
    testValues["Orb1Positron_OVGF_A"] = P3values[4]
    testValues["Orb1Positron_OVGF_B"] = P3values[5]
    testValues["Orb1Positron_OVGF_C"] = P3values[6]
    testValues["Orb1Positron_RENP3"] = P3values[7]
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
