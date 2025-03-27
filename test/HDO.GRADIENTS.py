#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues={
    "refGradHx" : [0.000000000000,1E-5],
    "refGradHy" :[ 0.061073535862,1E-5],
    "refGradHz" : [-0.041117203100,1E-5],
    "refGradDx" :[ 0.000000000000,1E-5],
    "refGradDy" : [-0.055359213575,1E-5],
    "refGradDz" : [-0.037187905984,1E-5],
    "refGradOx" :[ 0.000000000000,1E-5],
    "refGradOy" : [-0.005714322288,1E-5],
    "refGradOz" :[ 0.078305109084,1E-5],
}
    return refValues

def getTestValues(testValues,testName):
    grad=test.getGradient(testName)
    testValues["refGradHx"] = grad[0][0]
    testValues["refGradHy"] = grad[0][1]
    testValues["refGradHz"] = grad[0][2]
    testValues["refGradDx"] = grad[1][0]
    testValues["refGradDy"] = grad[1][1]
    testValues["refGradDz"] = grad[1][2]
    testValues["refGradOx"] = grad[2][0]
    testValues["refGradOy"] = grad[2][1]
    testValues["refGradOz"] = grad[2][2]
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
