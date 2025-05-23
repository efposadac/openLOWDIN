#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"Orb5e-_P2" : [-10.140654,1E-3],
"Orb6e-_P2" : [4.925340,1E-3],
"Orb1H1a_P2" : [-17.007032,1E-3],
"Orb2H1a_P2" : [-57.734065,1E-2],
"Orb1H1b_P2" : [-17.007032,1E-3],
"Orb2H1b_P2" : [-57.734065,1E-2]
}
    return refValues

def getTestValues(testValues,testName):
    testValues["Orb1H1a_P2"] = test.getP2orbEnergy(testName,"H-A_1",1)
    testValues["Orb2H1a_P2"] = test.getP2orbEnergy(testName,"H-A_1",2)
    testValues["Orb1H1b_P2"] = test.getP2orbEnergy(testName,"H-B_1",1)
    testValues["Orb2H1b_P2"] = test.getP2orbEnergy(testName,"H-B_1",2)
    testValues["Orb5e-_P2"] = test.getP2orbEnergy(testName,"E-",5)
    testValues["Orb6e-_P2"] = test.getP2orbEnergy(testName,"E-",6)
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
