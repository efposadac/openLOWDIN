#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"KS energy" : [-0.454704120729,1E-6],
"U+/E- Corr energy" : [-0.061929601444,1E-4],
"U+/ext pot energy" : [0.010775182673,1E-4],
"E-/ext pot energy" : [0.000679766944,1E-4],
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["KS energy"] = test.getSCFTotalEnergy(testName)
    testValues["U+/E- Corr energy"] = test.getDFTCorrEnergy(testName,"E-ALPHA","U+")
    testValues["U+/ext pot energy"] = test.getSCFExtPotEnergy(testName,"U+")
    testValues["E-/ext pot energy"] = test.getSCFExtPotEnergy(testName,"E-ALPHA")
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
