#!/usr/bin/env python
import lowdinTestFunctions as test
import os
import sys
from colorstring import *

def getDefaultLowdinBin():
    lowdinbin = "openlowdin"
    if(os.path.isfile("../CONFIG")):
        configFile = open("../CONFIG", "r")
        for line in configFile:
            if "EXENAME" in line:
                lowdinbin = line.split()[2]
                break
    return lowdinbin

def runCalculation(lowdinbin,inputName):
    status = os.system(lowdinbin + " -i " + inputName)
    if status:
        print(testName + str_red(" ... NOT OK"))
        sys.exit(1)

def verifyResults(refValues,testValues):
    passTest = True

    for value in refValues:
        diffValue = abs(refValues[value][0] - testValues[value]) 
        if ( diffValue <= refValues[value][1] ):
            passTest = passTest * True
        else :
            passTest = passTest * False
            print("%s %.8f %.8f %.2e" % ( value, refValues[value][0], testValues[value], diffValue))

    return passTest

def exitPrintResults(passTest,testName):
    if passTest:
        print(testName + str_green(" ... OK"))
        sys.exit(0)
    else:
        print(testName + str_red(" ... NOT OK"))
        sys.exit(1)

def getHFTotalEnergy(outputRead):
    energy=1.0E16
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "TOTAL ENERGY =" in line:
            energy = float(line.split()[3])
            break
    return energy

def getHFeigenvalues(outputRead,species,number):
    eigenval=1.0E16
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "Eigenvalues for: "+species in line:
            for j in range(i,len(outputRead)):
                linej = outputRead[j]
                if len(linej.split()) > 0:
                    if str(number) in linej.split()[0]:
                        eigenval = float(linej.split()[1])
                        break
    return eigenval
