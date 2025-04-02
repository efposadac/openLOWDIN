#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"HF energy" : [-534.848423377576,1E-6],
"Embedded HF energy" : [-534.846714877471,1E-6],
"MP2 energy" : [-535.073024924670,1E-6],
"P2 H_1" : [-15.6191,1E-3]
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["Embedded HF energy"] = test.getEmbeddedHFEnergy(testName)
    testValues["MP2 energy"] = test.getMP2Energy(testName)
    testValues["P2 H_1"] = test.getP2orbEnergy(testName,"H_1",1)
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
