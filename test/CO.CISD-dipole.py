#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"HF energy" : [-112.737336707933,1E-8],
"CISD energy" : [-112.957030012578,1E-6],
"HF dipole z" : [-0.33130728,1E-3],
"CISD dipole z" :  [0.13640258,1E-3],
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["CISD energy"] = test.getCIEnergy(testName,1)
    testValues["HF dipole"] = test.getSCFDipole(testName,"total","A.U.")
    testValues["HF dipole z"] = test.getSCFDipole(testName,"z","DEBYE")
    testValues["CISD dipole"] = test.getCIDipole(testName,"total","A.U.")
    testValues["CISD dipole z"] = test.getCIDipole(testName,"z","DEBYE")
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
