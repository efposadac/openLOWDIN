#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
    "Total mass" :                  [4.0383E+04, 5E-4],
    "E- Kinetic energy" :           [36.558549208752,1E-6],
    "MUON Kinetic energy" :           [4.823223401186,1E-6],
    "C_12 Kinetic energy" :           [0.020362665915,1E-6],
    "HE_4 Kinetic energy" :           [0.042798597334,1E-6],
    "H_1 Kinetic energy" :            [0.018810346458,1E-6],
    "H_2 Kinetic energy" :            [0.013654378944,1E-6],
    "H_3 Kinetic energy" :            [0.011128664979,1E-6],
    "HF energy" :                     [-74.168958869716,1E-8]
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["Total mass"] = test.getTotalMass(testName)
    testValues["E- Kinetic energy"] = test.getSCFKineticEnergy(testName,"E-")
    testValues["MUON Kinetic energy"] = test.getSCFKineticEnergy(testName,"MUON")
    testValues["C_12 Kinetic energy"] = test.getSCFKineticEnergy(testName,"C_12")
    testValues["HE_4 Kinetic energy"] = test.getSCFKineticEnergy(testName,"HE_4")
    testValues["H_1 Kinetic energy"] = test.getSCFKineticEnergy(testName,"H_1")
    testValues["H_2 Kinetic energy"] = test.getSCFKineticEnergy(testName,"H_2")
    testValues["H_3 Kinetic energy"] = test.getSCFKineticEnergy(testName,"H_3")
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
