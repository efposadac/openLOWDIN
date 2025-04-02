#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"HF energy" : [-75.93125555356,1E-8],
"MP2 energy" : [-76.095442398277,1E-5],
"e-AlphaBetaRepulsion" : [23.2395577028,1E-3]
}
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["MP2 energy"] = test.getMP2Energy(testName)
    testValues["e-AlphaBetaRepulsion"] = test.getSCFCoulombEnergy(testName,"E-ALPHA","E-BETA")
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
