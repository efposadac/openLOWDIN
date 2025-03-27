#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"HF energy" : [-75.895520937848,1E-8],
"HF dipole" : [1.06383820,1E-7],
"HF quadrupole xx" : [-7.21841981,1E-5],
"HF quadrupole yy" : [-4.02024678,1E-5],
"HF quadrupole zz" : [-6.20812552,1E-5],
"CI 1" : [-75.911063510564,1E-7],
"CI dipole" : [1.01277057,1E-7],
"CI quadrupole xx" : [-7.26705299,1E-5],
"CI quadrupole yy" : [-4.21742308,1E-5],
"CI quadrupole zz" : [-6.30621482,1E-5]
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["CI 1"] = test.getCIEnergy(testName,1)
    testValues["HF dipole"] = test.getSCFDipole(testName,"total","A.U.")
    testValues["CI dipole"] = test.getCIDipole(testName,"total","A.U.")
    quad=test.getSCFQuadrupole(testName)
    testValues["HF quadrupole xx"] = quad[0]
    testValues["HF quadrupole yy"] = quad[1]
    testValues["HF quadrupole zz"] = quad[2]
    quad=test.getCIQuadrupole(testName)
    testValues["CI quadrupole xx"] = quad[0]
    testValues["CI quadrupole yy"] = quad[1]
    testValues["CI quadrupole zz"] = quad[2]
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
