#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"HF energy" : [-153.536301213424,1E-6],
"Embedded HF energy" : [-153.118946505217,1E-6],
"CISD energy" : [-153.310526831872,1E-6],
"MP2 energy" : [-153.352174869549,1E-6],
"P2 E+" : [-5.328000,1E-3]
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["Embedded HF energy"] = test.getEmbeddedHFEnergy(testName)
    testValues["CISD energy"] = test.getCIEnergy(testName,1)
    testValues["MP2 energy"] = test.getMP2Energy(testName)
    testValues["P2 E+"] = test.getP2orbEnergy(testName,"E+",1)
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
