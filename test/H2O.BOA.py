#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"HF energy" : [-76.010288769789,1E-8],
"HF dipole" : [1.01124369,1E-7],
"HF quadrupole xx" : [-7.25729057,1E-5],
"HF quadrupole yy" : [-4.07831587,1E-5],
"HF quadrupole zz" : [-6.25445146,1E-5],
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["HF dipole"] = test.getSCFDipole(testName,"total","A.U.")
    quad=test.getSCFQuadrupole(testName)
    testValues["HF quadrupole xx"] = quad[0]
    testValues["HF quadrupole yy"] = quad[1]
    testValues["HF quadrupole zz"] = quad[2]
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
