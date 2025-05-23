#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues={
    "MP2 energy" : [-100.301839364887,1E-6]                
}
    return refValues

def getTestValues(testValues,testName):
    testValues["MP2 energy"] = test.getMP2Energy(testName)
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
