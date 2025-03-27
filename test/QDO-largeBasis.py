#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
        "HF energy" : [0.00,1E-8],
        "QDO zero energy" : [75.00,1E-8],
        "KT 2" : [50.0,1E-7],
        "KT 5" : [100.0,1E-7],
        "RMS R_c" : [0.002139,1E-4]        
    }                       
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
    testValues["KT 2"] = test.getHFeigenvalues(testName,"QP+",2)
    testValues["KT 5"] = test.getHFeigenvalues(testName,"QP+",5)
    testValues["QDO zero energy"] =test.getSCFTotalEnergyPlusQDO(testName)
    testValues["RMS R_c"] = test.getRMSradius(testName,"QP+")
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
