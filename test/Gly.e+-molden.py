#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"HF energy" : [-281.379249259318,1E-8],
"e+HOMO molden" : [-3.48415E-03,1E-4],
"e-HOMO molden" : [-5.40431E-01,1E-4],
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["e+HOMO molden"] = test.getHOMOmolden(testName,"E+")
    testValues["e-HOMO molden"] = test.getHOMOmolden(testName,"E-")
    testValues["e-HOMO"] = test.getHFeigenvalues(testName,"E-",1)
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
