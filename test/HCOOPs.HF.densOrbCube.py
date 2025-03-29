#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"HF energy" : [-188.362545831570,1E-7],
"Num E- in density cube" : [24.0,1E-1],
"Num E+ in density cube" : [1.0,1E-2],
"Num E+ in orbital cube" : [1.0,1E-2],
"Num E- in orbital cube" : [1.0,1E-2],
"Num E+ in 2D orbital plot" : [0.00021298,1E-4],
"Num E+ in 3D orbital plot" : [0.03267345,3E-2],
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["Num E- in density cube"] = test.getParticlesInDensCube(testName,"E-")
    testValues["Num E+ in density cube"] = test.getParticlesInDensCube(testName,"E+")
    testValues["Num E- in orbital cube"] = test.getParticlesInOrbCube(testName,"E-",12)
    testValues["Num E+ in orbital cube"] = test.getParticlesInOrbCube(testName,"E+",1)
    testValues["Num E+ in 2D orbital plot"] = test.getParticlesInOrbPlot2D(testName,"E+",1,False)
    testValues["Num E+ in 3D orbital plot"] = test.getParticlesInOrbPlot3D(testName,"E+",1,False)
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
