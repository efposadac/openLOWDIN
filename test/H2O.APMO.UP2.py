#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"HF energy" : [-75.930658814297,1E-8],
"Orb5alpha_P2" : [-10.4052,1E-3],
"Orb5alpha_SCS_P2" : [-10.7531,1E-3],
"Orb5alpha_SOS_P2" : [-10.9292,1E-3],
"Orb6alpha_P2" : [3.4416,1E-3],
"Orb1H1a_P2" : [-16.9153,1E-3],
"Orb2H1a_P2" : [-50.7107,1E-2],
"Orb1H1b_P2" : [-16.9153,1E-3],
"Orb2H1b_P2" : [-50.7107,1E-2]
}
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["Orb1H1a_P2"] = test.getP2orbEnergy(testName,"H-A_1",1)
    testValues["Orb2H1a_P2"] = test.getP2orbEnergy(testName,"H-A_1",2)
    testValues["Orb1H1b_P2"] = test.getP2orbEnergy(testName,"H-B_1",1)
    testValues["Orb2H1b_P2"] = test.getP2orbEnergy(testName,"H-B_1",2)
    testValues["Orb5alpha_P2"] = test.getP2orbEnergy(testName,"E-ALPHA",5)
    testValues["Orb5alpha_SCS_P2"] = test.getP2orbScaledEnergy(testName,"E-ALPHA",5,"SCS")
    testValues["Orb5alpha_SOS_P2"] = test.getP2orbScaledEnergy(testName,"E-ALPHA",5,"SOS")
    testValues["Orb6alpha_P2"] = test.getP2orbEnergy(testName,"E-ALPHA",6)
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
