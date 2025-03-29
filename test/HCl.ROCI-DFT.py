#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"KS energy" : [-460.740603429054,1E-6],
"CI 1" : [-460.763930172830,1E-5],
"CI 2" : [-460.763759214512,1E-5],
"H_1 Kin 1" : [0.004037448408,1E-5],
"H_1 Kin 2" : [0.004123732853,1E-5],
"E-/H_1 Corr 1" : [-0.026900920283,1E-4],
"E-/H_1 Corr 2" : [-0.026900920393,1E-4],
"scaled CI 1" : [-460.753060151213,1E-5],
"scaled CI 2" : [-460.752978293672,1E-5],
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["KS energy"] = test.getSCFTotalEnergy(testName)
    testValues["CI 1"] = test.getCIEnergy(testName,1)
    testValues["CI 2"] = test.getCIEnergy(testName,2)
    testValues["H_1 Kin 1"] = test.getNOCIKineticEnergy(testName,"H_1",1)
    testValues["H_1 Kin 2"] = test.getNOCIKineticEnergy(testName,"H_1",2)
    testValues["E-/H_1 Corr 1"] = test.getNOCIDFTcorrEnergy(testName,"E-","H_1",1)
    testValues["E-/H_1 Corr 2"] = test.getNOCIDFTcorrEnergy(testName,"E-","H_1",2)
    testValues["scaled CI 1"] = test.getNOCIDFTscaledEnergy(testName,1)
    testValues["scaled CI 2"] = test.getNOCIDFTscaledEnergy(testName,2)
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
