#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
        "HF energy" : [-199.396255393556,1E-8],
        "NOCI 1" : [-199.410508983518,1E-6],
        "NOCI 2" : [-199.396035449251,1E-6],
        "NOCI 4" : [-199.393067789142,1E-6],
        "Natural Occ 1 H_1 1" : [0.9954,1E-4],
        "Natural Orb 1 H_1 2" : [0.201157,1E-4],
        "Natural Orb 1 H_1 80" : [0.201157,1E-4],
        "Natural Occ 2 H_1 1" : [0.9892,1E-4],
        "Natural Orb 2 H_1 40" : [0.00000,1E-4],
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["NOCI 1"] = test.getCIEnergy(testName,1)
    testValues["NOCI 2"] = test.getCIEnergy(testName,2)
    testValues["NOCI 4"] = test.getCIEnergy(testName,4)
    testValues["Natural Occ 1 H_1 1"] = test.getNaturalOrbOcc(testName,"H_1",1,1)
    orbital = test.getNaturalOrb(testName,"H_1",1,1)
    testValues["Natural Orb 1 H_1 2"] = orbital[2-1] 
    testValues["Natural Orb 1 H_1 80"] = orbital[80-1] 
    testValues["Natural Occ 2 H_1 1"] = test.getNaturalOrbOcc(testName,"H_1",1,2)
    orbital = test.getNaturalOrb(testName,"H_1",1,2)
    testValues["Natural Orb 2 H_1 2"] = orbital[2-1] 
    testValues["Natural Orb 2 H_1 80"] = orbital[80-1] 
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
