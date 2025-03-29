#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"HF energy" : [-0.651377267238,1E-8],
"X0.5+ Ext Pot" : [0.005339817047,1E-4],
"Y0.5+ Ext Pot" : [0.005265789260,1E-4],
"X0.5+/Y0.5+ Hartree" : [0.003393498016,1E-4],
"X0.5+/Fixed interact." : [-0.025239485153,1E-4],
"Y0.5+/Fixed interact." : [-0.010156190638,1E-4]
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["X0.5+ Ext Pot"] = test.getSCFExtPotEnergy(testName,"X0.5+")
    testValues["Y0.5+ Ext Pot"] = test.getSCFExtPotEnergy(testName,"Y0.5+")
    testValues["X0.5+/Y0.5+ Hartree"] = test.getSCFCoulombEnergy(testName,"X0.5+","Y0.5+")
    testValues["X0.5+/Fixed interact."] = test.getSCFPointEnergy(testName,"X0.5+")
    testValues["Y0.5+/Fixed interact."] = test.getSCFPointEnergy(testName,"Y0.5+")
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
