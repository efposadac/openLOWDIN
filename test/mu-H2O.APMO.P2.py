#!/usr/bin/env python
#The corresponding input file is testName.lowdin
#The functions setReferenceValues and getTestValues are specific for this test
#The common procedures are found in lowdinTestFunctions.py
import sys
import lowdinTestFunctions as test
def setReferenceValues():
    refValues = {
"HF energy" : [-6669.258866831546,1E-8],
#"Orb1mu_P2" : [-180478.762974,1E-3],
#"Orb2mu_P2" : [-45540.156168,1E-2],
"Orb1H1a_P2" : [-18.625816,1E-3],
"Orb2H1a_P2" : [-55.835238,1E-2],
"Orb1H2b_P2" : [-13.634227,1E-3],
"Orb2H2b_P2" : [-35.117567,1E-2]
}
    return refValues

def getTestValues(testValues,testName):
    testValues["HF energy"] = test.getSCFTotalEnergy(testName)
#    testValues["Orb1mu_P2"] = test.getP2orbEnergy(testName,"U-",1)
#    testValues["Orb2mu_P2"] = test.getP2orbEnergy(testName,"U-",2)
    testValues["Orb1H1a_P2"] = test.getP2orbEnergy(testName,"H-A_1",1)
    testValues["Orb2H1a_P2"] = test.getP2orbEnergy(testName,"H-A_1",2)
    testValues["Orb1H2b_P2"] = test.getP2orbEnergy(testName,"H-B_2",1)
    testValues["Orb2H2b_P2"] = test.getP2orbEnergy(testName,"H-B_2",2)
    return 

if __name__ == '__main__':
    testName = sys.argv[0][:-3]
    test.performTest(testName,setReferenceValues,getTestValues)
