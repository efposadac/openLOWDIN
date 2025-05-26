#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
        "FCI 1" : [-1.279372034617,1E-6],
        "FCI 2" : [-1.217524569872,1E-6],
        "Natural Occ 1 e+ 1 molden" : [0.9209583352,1E-4],
        "Natural Occ 2 e+ 1 molden" : [0.0246217069,1E-4],
        "Natural Occ 3 e+ 1 molden" : [0.0139684400,1E-4],
        "Natural Occ 1 e+ 2 molden" : [0.8617085763,1E-4],
        "Natural Occ 2 e+ 2 molden" : [0.0757075628,1E-4],
        "Natural Occ 3 e+ 2 molden" : [0.0108344857,1E-4],
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["FCI 1"] = test.getCIEnergy(testName,1)
    testValues["FCI 2"] = test.getCIEnergy(testName,2)
    testValues["Natural Occ 1 e+ 1 molden"] = test.getOccupMolden(testName,"E+",1)
    testValues["Natural Occ 2 e+ 1 molden"] = test.getOccupMolden(testName,"E+",2)
    testValues["Natural Occ 3 e+ 1 molden"] = test.getOccupMolden(testName,"E+",3)
    testValues["Natural Occ 1 e+ 2 molden"] = test.getOccupMolden(testName,"E+",1,2)
    testValues["Natural Occ 2 e+ 2 molden"] = test.getOccupMolden(testName,"E+",2,2)
    testValues["Natural Occ 3 e+ 2 molden"] = test.getOccupMolden(testName,"E+",3,2)
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
