#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"HF energy" : [-0.666783062050,1E-8],
"FCI 1" : [-0.743335966767,1E-8],
"FCI 2" : [-0.595807154865,1E-8],
"FCI 3" : [-0.595807154727,1E-8],
"Natural Occ 1 e+ 1" : [0.9189,1E-4],
"Natural Orb 1 e+ 1" : [0.005431,1E-4],
"Natural Orb 1 e+ 2" : [0.033634,1E-4],
"Natural Orb 1 e+ 3" : [0.039679,1E-4],
"Natural Orb 1 e+ 4" : [0.163628,1E-4],
"Natural Orb 1 e+ 5" : [0.617836,1E-4],
"Natural Orb 1 e+ 6" : [0.303937,1E-4]
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["CI 1"] = test.getCIEnergy(testName,1)
    testValues["CI 2"] = test.getCIEnergy(testName,2)
    testValues["FCI 1"] = test.getCIEnergy(testName,1)
    testValues["FCI 2"] = test.getCIEnergy(testName,2)
    testValues["FCI 3"] = test.getCIEnergy(testName,3)
    testValues["Natural Occ 1 e+ 1"] = test.getNaturalOrbOcc(testName,"E+",1)
    orbital = test.getNaturalOrb(testName,"E+",1)
    for j in range(1,6+1):
        testValues["Natural Orb 1 e+ "+str(j)] = orbital[j-1] 
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
