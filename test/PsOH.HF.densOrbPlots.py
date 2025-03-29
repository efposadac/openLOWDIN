#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"HF energy" : [-75.587683637788,1E-8],
"Num e- in 3D densplot" : [10.0,1E-0],
"Num e- in 3D orbplot" :  [2.0,5E-1],
"Num e+ in 3D densplot" : [1.0,5E-2],
"Num e+ in 3D orbplot" :  [1.0,5E-2],
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["Num e- in 3D densplot"] = test.getParticlesInDensPlot3D(testName,"E-",True)
    testValues["Num e+ in 3D densplot"] = test.getParticlesInDensPlot3D(testName,"E+",True)
    testValues["Num e- in 3D orbplot"] = test.getParticlesInOrbPlot3D(testName,"E-",3,True)
    testValues["Num e+ in 3D orbplot"] = test.getParticlesInOrbPlot3D(testName,"E+",1,True)
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
