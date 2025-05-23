#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues={
    "Orb1Positron_KT" :     [-3.918673,1E-3],
    "Orb1Positron_EP2" :    [-4.579523,1E-3],
    "Orb1Positron_P3" :     [-4.904336,1E-3],
    "Orb1Positron_EP3" :    [-4.830791,1E-3],
    "Orb1Positron_OVGF-A" : [-4.915085,1E-3],
    "Orb1Positron_OVGF-B" : [-4.945115,1E-3],
    "Orb1Positron_OVGF-C" : [-4.900723,1E-3],
    "Orb1Positron_RENP3" :  [-4.946967,1E-3]
    }
    return refValues

def getTestValues(testValues,testName):
    P3values = test.getP3results(testName,"E+",1)
    testValues["Orb1Positron_KT"] = P3values[0]
    testValues["Orb1Positron_EP2"] = P3values[1]
    testValues["Orb1Positron_P3"] = P3values[2]
    testValues["Orb1Positron_EP3"] = P3values[3]
    testValues["Orb1Positron_OVGF-A"] = P3values[4]
    testValues["Orb1Positron_OVGF-B"] = P3values[5]
    testValues["Orb1Positron_OVGF-C"] = P3values[6]
    testValues["Orb1Positron_RENP3"] = P3values[7]
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
