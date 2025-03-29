#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"HF energy" : [-99.635031860198,1E-8],
"Num e- in 2D densplot" : [10.0,1E-2],
"Num e- in 2D orbplot" :  [1.0,1E-2],
"Num e+ in 2D densplot" : [1.0,1E-3],
"Num e+ in 2D orbplot" :  [1.0,1E-2],
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["Num e- in 2D densplot"] = test.getParticlesInDensPlot2D(testName,"E-",True)
    testValues["Num e+ in 2D densplot"] = test.getParticlesInDensPlot2D(testName,"E+",True)
    testValues["Num e- in 2D orbplot"] = test.getParticlesInOrbPlot2D(testName,"E-",2,True)
    testValues["Num e+ in 2D orbplot"] = test.getParticlesInOrbPlot2D(testName,"E+",1,True)
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
