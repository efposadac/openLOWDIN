#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"HF energy" : [29.59021457,1E-8],
"TrCI energy" : [-0.49952730,1E-4],
"E-ALPHA Kinetic" : [0.49964274,1E-4],
"H_1 Kinetic" : [0.00028259,1E-5],
"E-ALPHA/H_1 Hartree" : [-0.99945262,1E-4],
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["TrCI energy"] = test.getTrCIEnergyContribution(testName,"Hamiltonian","","")
    testValues["E-ALPHA Kinetic"] = test.getTrCIEnergyContribution(testName,"Kinetic","E-ALPHA","")
    testValues["H_1 Kinetic"] = test.getTrCIEnergyContribution(testName,"Kinetic","H_1","")
    testValues["E-ALPHA/H_1 Hartree"] = test.getTrCIEnergyContribution(testName,"Hartree","E-ALPHA","H_1")
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
