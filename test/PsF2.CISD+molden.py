#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"HF energy" : [-198.967220107085,1E-8],
"CISD+ energy" : [-198.977666606058,1E-6],
"Natural Occ 10 e-alpha 1 molden" : [0.9987,1E-4],
"Natural Occ 11 e-alpha 1 molden" : [0.0017,1E-4],
"Natural Occ 9 e-beta 1 molden" : [0.9648,1E-4],
"Natural Occ 10 e-beta 1 molden" : [0.0344,1E-4],
"Natural Occ 1 e+ 1 molden" : [0.9653,1E-4],
"Natural Occ 2 e+ 1 molden" : [0.0327,1E-4],
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["CISD+ energy"] = test.getCIEnergy(testName,1)
    testValues["Natural Occ 1 e+ 1"] = test.getNaturalOrbOcc(testName,"E+",1)
    testValues["Natural Occ 1 e+ 1 molden"] = test.getOccupMolden(testName,"E+",1)
    testValues["Natural Occ 2 e+ 1 molden"] = test.getOccupMolden(testName,"E+",2)
    testValues["Natural Occ 10 e-alpha 1 molden"] = test.getOccupMolden(testName,"E-ALPHA",10)
    testValues["Natural Occ 11 e-alpha 1 molden"] = test.getOccupMolden(testName,"E-ALPHA",11)
    testValues["Natural Occ 9 e-beta 1 molden"] = test.getOccupMolden(testName,"E-BETA",9)
    testValues["Natural Occ 10 e-beta 1 molden"] = test.getOccupMolden(testName,"E-BETA",10)
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
