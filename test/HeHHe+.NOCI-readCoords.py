#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
        "HF energy" : [-5.753361915874,1E-8],
        "NOCI 1" : [-5.757017190944,1E-6],
        "NOCI 2" : [-5.756828937998,1E-6],
        "Natural Occ 1 H_1 1" : [0.9900,1E-4],
        "Natural Orb 1 H_1 3" : [14.966751,1E-4],
        "Natural Orb 1 H_1 31" : [14.966751,1E-4],
        "Natural Occ 2 H_1 1" : [0.9899,1E-4],
        "Natural Orb 2 H_1 3" : [15.017286,1E-4],
        "Natural Orb 2 H_1 31" : [15.017286,1E-4],
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["NOCI 1"] = test.getCIEnergy(testName,1)
    testValues["NOCI 2"] = test.getCIEnergy(testName,2)
    testValues["Natural Occ 1 H_1 1"] = test.getNaturalOrbOcc(testName,"H_1",1,1)
    orbital = test.getNaturalOrb(testName,"H_1",1,1)
    testValues["Natural Orb 1 H_1 3"] = orbital[3-1] 
    testValues["Natural Orb 1 H_1 31"] = orbital[31-1] 
    testValues["Natural Occ 2 H_1 1"] = test.getNaturalOrbOcc(testName,"H_1",1,2)
    orbital = test.getNaturalOrb(testName,"H_1",1,2)
    testValues["Natural Orb 2 H_1 3"] = orbital[3-1] 
    testValues["Natural Orb 2 H_1 31"] = orbital[31-1] 
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
