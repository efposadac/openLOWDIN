#!/usr/bin/env python
#This file defines the procedures common to the openlowdin software tests
#Felix Moncada, Mar/2025
import os
import sys

def performTest(testName,setRefValuesFunc,getTestValuesFunc):
    # Run calculation
    if len(sys.argv)==2:
        lowdinbin = sys.argv[1]
    else:
        lowdinbin = getDefaultLowdinBin()
    runLowdinCalculation(lowdinbin,testName)

    # Get reference values and tolerance
    refValues=setRefValuesFunc()

    # Get Test results
    testValues = dict(refValues) #copy 
    for value in testValues: 
        testValues[value] = 0 #reset
    getTestValuesFunc(testValues,testName)

    # Check results and finish, return error if test fail
    passTest=verifyResults(refValues,testValues)
    exitPrintResults(passTest,testName)

def str_green(string):
    return chr(27) + "[1;32m" + string + chr(27) + "[0m"

def str_red(string):
    return chr(27) + "[1;31m" + string + chr(27) + "[0m"

def getDefaultLowdinBin():
    lowdinbin = "openlowdin"
    if(os.path.isfile("../CONFIG")):
        configFile = open("../CONFIG", "r")
        for line in configFile:
            if "EXENAME" in line:
                lowdinbin = line.split()[2]
                break
    return lowdinbin

def runLowdinCalculation(lowdinbin,testName):
    status = os.system(lowdinbin + " -i " + testName+".lowdin")
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

def getSCFTotalEnergy(testName):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    energy=1.0E16
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "TOTAL ENERGY =" in line:
            energy = float(line.split()[3])
            break
    output.close()
    return energy

def getDFTTotalExcCorrEnergy(testName):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    energy=1.0E16
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "Total Exchange Correlation energy" in line:
            energy = float(line.split()[5])
            break
    output.close()
    return energy

def getSCFTotalExtPotEnergy(testName):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    energy=1.0E16
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "Total External Potential energy =" in line:
            energy = float(line.split()[5])
            break
    output.close()
    return energy

def getSCFKineticEnergy(testName,species):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    energy=1.0E16
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if species+" Kinetic energy" in line:
            energy = float(line.split()[4])
            break
    output.close()
    return energy

def getSCFPointEnergy(testName,species):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    energy=1.0E16
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if species+"/Fixed interact." in line:
            energy = float(line.split()[4])
            break
    output.close()
    return energy

def getSCFExtPotEnergy(testName,species):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    energy=1.0E16
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if species+" Ext Pot" in line:
            energy = float(line.split()[5])
            break
    output.close()
    return energy

def getSCFExchangeEnergy(testName,species):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    energy=1.0E16
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if species+" Exchange" in line:
            energy = float(line.split()[4])
            break
    output.close()
    return energy

def getSCFCoulombEnergy(testName,species,ospecies):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    energy=1.0E16
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if species+"/"+ospecies+" Hartree" in line:
            energy = float(line.split()[4])
            break
    output.close()
    return energy

def getDFTExcCorrEnergy(testName,species):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    energy=1.0E16
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if species+" Exc.Corr. energy" in line:
            energy = float(line.split()[4])
            break
    output.close()
    return energy

def getDFTCorrEnergy(testName,species,ospecies):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    energy=1.0E16
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if species+"/"+ospecies+" Corr. energy" in line:
            energy = float(line.split()[4])
            break
    output.close()
    return energy

def getEmbeddedHFEnergy(testName):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    energy=1.0E16
    flag=1
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "TOTAL ENERGY =" in line and flag==1:
            flag=2
            continue
        if "TOTAL ENERGY =" in line and flag==2:
            energy = float(line.split()[3])
            break
    output.close()
    return energy

def getCIEnergy(testName,state):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    energy=1.0E16
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "STATE:   "+str(state)+" ENERGY =" in line:
            energy = float(line.split()[4])
            break
    output.close()
    return energy

def getMP2Energy(testName):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    energy=1.0E16
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "E(MP2) =" in line:
            energy = float(line.split()[2])
            break
    output.close()
    return energy

def getEN2Energy(testName):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    energy=1.0E16
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "E(EN2) =" in line:
            energy = float(line.split()[2])
            break
    output.close()
    return energy

def getNSEN2Energy(testName):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    energy=1.0E16
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "E(NS-EN2) =" in line:
            energy = float(line.split()[2])
            break
    output.close()
    return energy

def getSCIPT2Energy(testName):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    energy=1.0E16
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "E_SCI + E_PT2 :" in line:
            energy = float(line.split()[4])
            break
    output.close()
    return energy

def getP2orbEnergy(testName,species,number):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    energy=1.0E16
    speciesFlag=False
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "Results for spin-orbital: "+str(number)+" of species:" in line and species in line :
            speciesFlag=True
        if "Optimized second order pole:" in line and speciesFlag:
            energy = float(line.split()[4])
            break
    output.close()
    return energy

def getP2orbScaledEnergy(testName,species,number,scalingType):
    if scalingType == "SCS":
        query="FactorOS: 1.20000 FactorSS: 0.33333"
    elif scalingType == "SOS":
        query="FactorOS: 1.30000 FactorSS: 0.00000"
    else:
        return "error"
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    energy=1.0E16
    speciesFlag=False
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "Results for spin-orbital: "+str(number)+" of species:" in line and species in line :
            speciesFlag=True
        if query in line :
            energy = float(outputRead[i+1].split()[4])
            break
    output.close()
    return energy

def getHFeigenvalues(testName,species,number):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
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
    output.close()
    return eigenval

def getSCFTotalEnergyPlusQDO(testName):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    energy=1.0E16
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "TOTAL ENERGY plus QDO ZPE =" in line:
            energy = float(line.split()[6])
            break
    output.close()
    return energy

def getRMSradius(testName,species):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    radius=1.0E16
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "EXPECTED RMS RADIUS OF QUANTUM SPECIES" in line:
            for j in range(i,len(outputRead)):
                linej = outputRead[j]
                if "QP+" in linej:
                    radius= float(linej.split()[1])
                    break
    output.close()
    return radius

def getTotalMass(testName):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    mass=1.0E16
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "MASS (m_e)" in line:
            mass = float(line.split()[3])
            break
    output.close()
    return mass

def getSCFDipole(testName,component,units):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    dipole=1.0E16
    dipoleflag=False
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "DIPOLE: ("+units+")" in line:
            dipoleflag=True
        if "Total Dipole:" in line and dipoleflag:
            if component == "x":
                dipole= float(line.split()[2])
            elif component == "y":
                dipole= float(line.split()[3])
            elif component == "z":
                dipole= float(line.split()[4])
            else:
                dipole= float(line.split()[5])                
            break
    output.close()
    return dipole

def getCIDipole(testName,component,units):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    dipole=1.0E16
    dipoleflag=False
    ciflag=False
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "We are calculating properties for E-ALPHA in the CI ground state" in line:
            ciflag=True
        if "DIPOLE: ("+units+")" in line and ciflag:
            dipoleflag=True
        if "Total Dipole:" in line and dipoleflag and ciflag:
            if component == "x":
                dipole= float(line.split()[2])
            elif component == "y":
                dipole= float(line.split()[3])
            elif component == "z":
                dipole= float(line.split()[4])
            else:
                dipole= float(line.split()[5])                
            break
    output.close()
    return dipole

def getSCFQuadrupole(testName):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    quadrupole=[]
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "Total Quadrupole:" in line:
            for j in range(2,len(line.split())):
                quadrupole.append(float(line.split()[j]))
            break
    output.close()
    return quadrupole

def getCIQuadrupole(testName):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    ciflag=False
    quadrupole=[]
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "We are calculating properties for E-ALPHA in the CI ground state" in line:
            ciflag=True
        if "Total Quadrupole:" in line and ciflag:
            for j in range(2,len(line.split())):
                quadrupole.append(float(line.split()[j]))
            break
    output.close()
    return quadrupole

def getParticlesInGrid(testName,species):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    number=1.0E16
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "Number of "+species+" particles in the final grid" in line:
            number = float(line.split()[8])
            break
    output.close()
    return number

def getOccupMolden(testName,species,number):
    molden = open(testName+"."+species+".molden", "r")
    moldenRead = molden.readlines()
    v=0
    eigenv=1.0E16
    for i in range(0,len(moldenRead)):
        line = moldenRead[i]
        if "Ene=" in line:
            v+=1
        if "Occup=" in line and v == number :
            eigenv = float(line.split()[1])
            break
    molden.close()
    return eigenv

def getHOMOmolden(testName,species):
    molden = open(testName+"."+species+".molden", "r")
    moldenRead = molden.readlines()
    v=0
    eigenv=[]
    energy=1.0E16
    for i in range(0,len(moldenRead)):
        line = moldenRead[i]
        if "Ene=" in line:
            eigenv.append(float(line.split()[1]))
            v+=1
        if "Occup=" in line and "0.0" in line :
            energy = eigenv[v-2]
            break
    molden.close()
    return energy

def getOccupiedOrbitalsMolden(testName,species):
    molden = open(testName+"."+species+".molden", "r")
    moldenRead = molden.readlines()
    v=-1
    flag=0
    eigenvec=[]
    for i in range(0,len(moldenRead)):
        line = moldenRead[i]
        if "Ene=" in line:
            eigenvec.append([])
            v+=1
        if flag==1 and "=" not in line:
            eigenvec[v].append(float(line.split()[1]))
        if "[MO]" in line:
            flag=1
        if "Occup=" in line and "0.0" in line :
            break
    molden.close()
    return eigenvec

def getSCFiterations(testName):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    number=1.0E16
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if " Total energy converged after" in line:
            number = float(line.split()[4])
            break
        if " Total density converged after" in line:
            number = float(line.split()[4])
            break
        if " Total energy and density converged after" in line:
            number = float(line.split()[6])
            break        
    output.close()
    return number

def getParticlesInDensPlot2D(testName,species,radialNorm):
    densplot = open(testName+"."+species+".2D.dens", "r")
    densplotRead = densplot.readlines()
    sumPart=0
    for i in range(0,len(densplotRead)):
        line = densplotRead[i]
        if i > 1:
            values = line.split()
            if i == 3: x1=float(values[0])
            if i == 4: x2=float(values[0])
            if len(values) > 1:
                if radialNorm and float(values[0]) > 0.0:
                    sumPart+=float(values[0])**2*float(values[1])
                elif not radialNorm:
                    sumPart+=float(values[1])
    if radialNorm:
        sumPart=4.0*3.14159265359*sumPart*(x2-x1)
    else:
        sumPart=sumPart*(x2-x1)**3
    densplot.close()
    return sumPart

def getParticlesInOrbPlot2D(testName,species,number,radialNorm):
    orbplot = open(testName+"."+species+".2D.orb"+str(number), "r")
    orbplotRead = orbplot.readlines()
    sumPart=0
    for i in range(0,len(orbplotRead)):
        line = orbplotRead[i]
        if i > 1:
            values = line.split()
            if i == 3: x1=float(values[0])
            if i == 4: x2=float(values[0])
            if len(values) > 1:
                if radialNorm and float(values[0]) > 0.0:
                    sumPart+=float(values[0])**2*float(values[1])**2
                elif not radialNorm:
                    sumPart+=float(values[1])**2
    if radialNorm:
        sumPart=4.0*3.14159265359*sumPart*(x2-x1)
    else:
        sumPart=sumPart*(x2-x1)**3
    orbplot.close()
    return sumPart

def getParticlesInDensPlot3D(testName,species,radialNorm):
    densplot = open(testName+"."+species+".3D.dens", "r")
    densplotRead = densplot.readlines()
    sumPart=0
    for i in range(0,len(densplotRead)):
        line = densplotRead[i]
        if i > 1:
            values = line.split()
            if i == 3: x1=float(values[1])
            if i == 4: x2=float(values[1])
            if radialNorm and len(values) > 1:
                sumPart+=(float(values[0])**2+float(values[1])**2)**(1.0/2.0)*float(values[2])
            elif len(values) > 1:
                sumPart+=float(values[2])
    if radialNorm:
        sumPart=2.0*sumPart*(x2-x1)**2
    else:
        sumPart=sumPart*(x2-x1)**3
    densplot.close()
    return sumPart

def getParticlesInOrbPlot3D(testName,species,number,radialNorm):
    orbplot = open(testName+"."+species+".3D.orb"+str(number), "r")
    orbplotRead = orbplot.readlines()
    sumPart=0
    for i in range(0,len(orbplotRead)):
        line = orbplotRead[i]
        if i > 1:
            values = line.split()
            if i == 3: x1=float(values[1])
            if i == 4: x2=float(values[1])
            if radialNorm and len(values) > 1:
                sumPart+=(float(values[0])**2+float(values[1])**2)**(1.0/2.0)*float(values[2])**2
            elif len(values) > 1:
                sumPart+=float(values[2])**2
    if radialNorm:
        sumPart=2.0*sumPart*(x2-x1)**2
    else:
        sumPart=sumPart*(x2-x1)**3
    orbplot.close()
    return sumPart

def getParticlesInDensCube(testName,species):
    cube = open(testName+"."+species+".dens.cub", "r")
    cubeRead = cube.readlines()
    sumPart=0
    for i in range(0,len(cubeRead)):
        line = cubeRead[i]
        if i == 3: step=float(line.split()[1])
        if i > 10:
            values = line.split()
            for j in range(0,len(values)):
                sumPart+=float(values[j])
    sumPart=sumPart*step**3
    cube.close()
    return sumPart

def getParticlesInOrbCube(testName,species,number):
    cube = open(testName+"."+species+".orb"+str(number)+".cub", "r")
    cubeRead = cube.readlines()
    sumPart=0
    for i in range(0,len(cubeRead)):
        line = cubeRead[i]
        if i == 3: step=float(line.split()[1])
        if i > 10:
            values = line.split()
            for j in range(0,len(values)):
                sumPart+=float(values[j])**2
    sumPart=sumPart*step**3
    cube.close()
    return sumPart

def getP3results(testName,species,number):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    eigenval=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]    
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "SUMMARY OF PROPAGATOR RESULTS FOR THE SPIN-ORBITAL: "+str(number)+"  OF SPECIES:"+species in line:
            eigenval[0]=float(outputRead[i + 4].split()[1]) #KT
            eigenval[1]=float(outputRead[i + 5].split()[1]) #EP2
            eigenval[2]=float(outputRead[i + 6].split()[1]) #P3
            eigenval[3]=float(outputRead[i + 7].split()[1]) #EP3
            eigenval[4]=float(outputRead[i + 8].split()[1]) #OVGFA
            eigenval[5]=float(outputRead[i + 9].split()[1]) #OVGFB
            eigenval[6]=float(outputRead[i + 10].split()[1]) #OVGFC
            eigenval[7]=float(outputRead[i + 11].split()[1]) #RENP3
            break
    output.close()
    return eigenval

def getNaturalOrbOcc(testName,species,number):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "  Natural Orbitals in state:            "+str(number)+"  for: "+species in line:
            occupation = float(outputRead[i+2].split()[number-1]) #Works only for the five first orbitals

    output.close()
    return occupation

def getNaturalOrb(testName,species,number):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    eigenvec=[]    
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "  Natural Orbitals in state:            "+str(number)+"  for: "+species in line:
            for j in range(1,len(outputRead)-i-3):
                linej = outputRead[i+3+j]
                if len(linej.split()) == 0:
                    basisSize=j-1
                    break
            for j in range(1,basisSize+1):
                linej = outputRead[i+3+j]
                eigenvec.append(abs(float(linej.split()[2+number]))) #Works only for the five first orbitals
            break
    output.close()
    return eigenvec

def getSpeciesOccupation(testName,species):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    flagC=0
    occupation=1E16
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "CONSTANTS OF COUPLING" in line:
            flagC=1
        if species in line and flagC==1:
            occupation=float(line.split()[4])
            break
    output.close()
    return occupation

def getSpeciesEta(testName,species):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    flagC=0
    eta=1E16
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "CONSTANTS OF COUPLING" in line:
            flagC=1
        if species in line and flagC==1:
            eta=float(line.split()[2])
            break
    output.close()
    return eta

def getHFCoefficient(testName):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    coeff=1.0E16
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "HF COEFFICIENT =" in line:
            coeff = abs(float(line.split()[3]))
            break
    output.close()
    return coeff

def getNOCIKineticEnergy(testName,species,number):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    energy=1.0E16
    stateFlag=False
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "STATE:   "+str(number)+" ENERGY =" in line:
            stateFlag=True
        if species+" Kinetic energy =" in line and stateFlag:
            energy = float(line.split()[4])
            break
    output.close()
    return energy

def getNOCIDFTcorrEnergy(testName,species,ospecies,number):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    energy=1.0E16
    stateFlag=False
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "STATE:   "+str(number)+" ENERGY =" in line:
            stateFlag=True
        if species+"/"+ospecies+" DFTcorrelation energy =" in line and stateFlag:
            energy = float(line.split()[4])
            break
    output.close()
    return energy

def getNOCIDFTscaledEnergy(testName,number):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    energy=1.0E16
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if "STATE:   "+str(number)+" SCALED ENERGY =" in line:
            energy = float(line.split()[5])
            break
    output.close()
    return energy

def getGradient(testName):
    output = open(testName+".out", "r")
    outputRead = output.readlines()
    grad=[]
    gradFlag=False
    natom=-1
    for i in range(0,len(outputRead)):
        line = outputRead[i]
        if len(line.split()) == 3 and gradFlag:
            grad.append([])
            natom+=1
            for j in range(0,3):
                grad[natom].append(float(line.split()[j]))
        if "dE/dx            dE/dy            dE/dz" in line:
            gradFlag=True
        if "END ENERGY GRADIENTS" in line:
            break
    output.close()
    return grad
