#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
    "HF energy" : [-0.651377261423,1E-8],
    "X0.5+ Ext Pot" : [0.010606635479,1E-4],
    "X0.5+/X0.5+ Hartree" : [1.217813835967,1E-4],
    "X0.5+ Exchange" : [-1.214419735313,1E-4],
    "X0.5+/Fixed interact." : [-0.035396418562,1E-4]
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["X0.5+ Ext Pot"] = test.getSCFExtPotEnergy(testName,"X0.5+")
    testValues["X0.5+/X0.5+ Hartree"] = test.getSCFCoulombEnergy(testName,"X0.5+","X0.5+")
    testValues["X0.5+ Exchange"] = test.getSCFExchangeEnergy(testName,"X0.5+")
    testValues["X0.5+/Fixed interact."] = test.getSCFPointEnergy(testName,"X0.5+")
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
