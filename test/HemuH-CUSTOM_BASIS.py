#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"HF energy" : [-343.383218426575,1E-8],
"U-HOMO" : [-371.889981903188,1E-1],
"H_1-HOMO" : [-1.019346186407,1E-4],
"He_4-HOMO" : [-652.365841581870,1E-1],
"e-HOMO" : [-0.585408097570,1E-4],
}                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["U-HOMO"] = test.getHFeigenvalues(testName,"MUON",1)
    testValues["H_1-HOMO"] = test.getHFeigenvalues(testName,"H_1",1)
    testValues["He_4-HOMO"] = test.getHFeigenvalues(testName,"HE_4",1)
    testValues["e-HOMO"] = test.getHFeigenvalues(testName,"E-",1)
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
