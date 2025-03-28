#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"HF energy" : [-1.255739154579,1E-8],
"FCI 1" : [-1.306431783235,1E-8],
"FCI 2" : [-1.291420173154,1E-8],
"FCI 3" : [-1.290426417940,1E-8],
"Natural Occ 1 H_1 1"  : [0.9919,1E-3],
"Natural Orb 1 H_1 1"  : [0.721928,1E-3],
"Natural Orb 1 H_1 2"  : [0.187901,1E-3],
"Natural Orb 1 H_1 3"  : [0.000000,1E-3],
"Natural Orb 1 H_1 4"  : [0.150257,1E-3],
"Natural Orb 1 H_1 5"  : [0.000000,1E-3],
"Natural Orb 1 H_1 6"  : [0.000000,1E-3],
"Natural Orb 1 H_1 7"  : [0.334924,1E-3],
"Natural Orb 1 H_1 8"  : [0.000000,1E-3],
"Natural Orb 1 H_1 9"  : [0.125434,1E-3],
"Natural Orb 1 H_1 10" : [0.000000,1E-3],
"Natural Orb 1 H_1 11" : [0.000000,1E-3],
"Natural Orb 1 H_1 12" : [0.116370,1E-3],
"Natural Orb 1 H_1 13" : [0.000000,1E-3],
"Natural Orb 1 H_1 14" : [0.117396,1E-3],
"Natural Orb 1 H_1 15" : [0.131109,1E-3],
"Natural Orb 1 H_1 16" : [0.000000,1E-3],
"Natural Orb 1 H_1 17" : [0.000000,1E-3],
"Natural Orb 1 H_1 18" : [0.127979,1E-3],
"Natural Orb 1 H_1 19" : [0.000000,1E-3],
"Natural Orb 1 H_1 20" : [0.161496,1E-3],
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["CI 1"] = test.getCIEnergy(testName,1)
    testValues["CI 2"] = test.getCIEnergy(testName,2)
    testValues["FCI 1"] = test.getCIEnergy(testName,1)
    testValues["FCI 2"] = test.getCIEnergy(testName,2)
    testValues["FCI 3"] = test.getCIEnergy(testName,3)
    testValues["Natural Occ 1 H_1 1"] = test.getNaturalOrbOcc(testName,"H_1",1)
    orbital = test.getNaturalOrb(testName,"H_1",1)
    for j in range(1,20+1):
        testValues["Natural Orb 1 H_1 "+str(j)] = orbital[j-1] 
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
