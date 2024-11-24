#!/usr/bin/env python
from __future__ import print_function
import os
import sys
from colorstring import *

if len(sys.argv)==2:
    lowdinbin = sys.argv[1]
else:
    lowdinbin = "lowdin2"

testName = sys.argv[0][:-3]
inputName = testName + ".lowdin"
outputName = testName + ".out"
molden1Name = testName + ".E-.molden"                                                                
# Reference values and tolerance

refValues = {
"HF energy" : [-76.056328834832,1E-8],
"e-HOMO" : [-5.08882E-01,1E-4],
"eigvec,1,1": [0.9939786,0.001],
"eigvec,1,2": [0.00321949,0.001],
"eigvec,1,3": [0.02029773,0.001],
"eigvec,1,8": [0.00208321,0.001],
"eigvec,1,11": [0.0017893,0.001],
"eigvec,1,14": [0.0015845,0.001],
"eigvec,1,18": [0.00156669,0.001],
"eigvec,1,20": [0.00113877,0.001],
"eigvec,1,24": [0.0033928,0.001],
"eigvec,1,26": [0.00158062,0.001],
"eigvec,1,48": [0.00110448,0.001],
"eigvec,1,56": [0.00263608,0.001],
"eigvec,1,58": [0.00160313,0.001],
"eigvec,1,65": [0.00106948,0.001],
"eigvec,1,66": [0.00125158,0.001],
"eigvec,2,1": [0.0078336,0.001],
"eigvec,2,2": [0.30870689,0.001],
"eigvec,2,3": [0.19694376,0.001],
"eigvec,2,4": [0.56914251,0.001],
"eigvec,2,5": [0.30272071,0.001],
"eigvec,2,8": [0.02142744,0.001],
"eigvec,2,11": [0.06790512,0.001],
"eigvec,2,14": [0.02891775,0.001],
"eigvec,2,17": [0.02847843,0.001],
"eigvec,2,18": [0.00191158,0.001],
"eigvec,2,19": [0.00182459,0.001],
"eigvec,2,20": [0.00238583,0.001],
"eigvec,2,24": [0.01976506,0.001],
"eigvec,2,25": [0.00425676,0.001],
"eigvec,2,26": [0.0098827,0.001],
"eigvec,2,30": [0.03229255,0.001],
"eigvec,2,31": [0.02079304,0.001],
"eigvec,2,32": [0.01922866,0.001],
"eigvec,2,38": [0.00116329,0.001],
"eigvec,2,44": [0.0022189,0.001],
"eigvec,2,48": [0.04645561,0.001],
"eigvec,2,51": [0.0155471,0.001],
"eigvec,2,54": [0.03489843,0.001],
"eigvec,2,57": [0.0016457,0.001],
"eigvec,2,65": [0.00194347,0.001],
"eigvec,2,66": [0.00124899,0.001],
"eigvec,2,67": [0.00679824,0.001],
"eigvec,2,71": [0.1195507,0.001],
"eigvec,2,72": [0.00111498,0.001],
"eigvec,2,73": [0.11955076,0.001],
"eigvec,2,74": [0.00111497,0.001],
"eigvec,3,7": [0.14103306,0.001],
"eigvec,3,10": [0.26698666,0.001],
"eigvec,3,13": [0.29854243,0.001],
"eigvec,3,16": [0.14680475,0.001],
"eigvec,3,23": [0.00612365,0.001],
"eigvec,3,29": [0.01201828,0.001],
"eigvec,3,35": [0.03167148,0.001],
"eigvec,3,37": [0.0020827,0.001],
"eigvec,3,40": [0.00100915,0.001],
"eigvec,3,43": [0.00271485,0.001],
"eigvec,3,47": [0.03543525,0.001],
"eigvec,3,50": [0.00106762,0.001],
"eigvec,3,53": [0.01304568,0.001],
"eigvec,3,62": [0.01280384,0.001],
"eigvec,3,64": [0.00520432,0.001],
"eigvec,3,71": [0.21875249,0.001],
"eigvec,3,72": [0.11175294,0.001],
"eigvec,3,73": [0.21875254,0.001],
"eigvec,3,74": [0.11175293,0.001],
"eigvec,4,1": [0.0029427,0.001],
"eigvec,4,2": [0.11064105,0.001],
"eigvec,4,3": [0.09231006,0.001],
"eigvec,4,4": [0.17254931,0.001],
"eigvec,4,5": [0.23227883,0.001],
"eigvec,4,8": [0.15648245,0.001],
"eigvec,4,11": [0.29035693,0.001],
"eigvec,4,14": [0.32398334,0.001],
"eigvec,4,17": [0.24327021,0.001],
"eigvec,4,18": [0.00390379,0.001],
"eigvec,4,20": [0.00378954,0.001],
"eigvec,4,25": [0.00339092,0.001],
"eigvec,4,30": [0.00946407,0.001],
"eigvec,4,31": [0.01003612,0.001],
"eigvec,4,32": [0.05371892,0.001],
"eigvec,4,38": [0.00166599,0.001],
"eigvec,4,44": [0.00343828,0.001],
"eigvec,4,48": [0.02940923,0.001],
"eigvec,4,51": [0.00174679,0.001],
"eigvec,4,54": [0.01368111,0.001],
"eigvec,4,56": [0.00429996,0.001],
"eigvec,4,57": [0.00137211,0.001],
"eigvec,4,58": [0.00364046,0.001],
"eigvec,4,65": [0.00181726,0.001],
"eigvec,4,66": [0.00196488,0.001],
"eigvec,4,67": [0.00609422,0.001],
"eigvec,4,71": [0.13676351,0.001],
"eigvec,4,72": [0.07750444,0.001],
"eigvec,4,73": [0.13676345,0.001],
"eigvec,4,74": [0.07750437,0.001],
"eigvec,5,6": [0.17966322,0.001],
"eigvec,5,9": [0.33034402,0.001],
"eigvec,5,12": [0.38223722,0.001],
"eigvec,5,15": [0.32672634,0.001],
"eigvec,5,22": [0.00385062,0.001],
"eigvec,5,28": [0.00134941,0.001],
"eigvec,5,34": [0.04247752,0.001],
"eigvec,5,36": [0.00350524,0.001],
"eigvec,5,39": [0.00103526,0.001],
"eigvec,5,42": [0.00117129,0.001],
"eigvec,5,46": [0.00723327,0.001],
"eigvec,5,49": [0.02445057,0.001],
"eigvec,5,52": [0.01703417,0.001],
"eigvec,5,63": [0.00285564,0.001],
"eigvec,5,69": [0.00492437,0.001]
}                       

testValues = dict(refValues) #copy 
for value in testValues: #reset
    testValues[value] = 0 #reset
    
# Run calculation

status = os.system(lowdinbin + " -i " + inputName)

if status:
    print(testName + str_red(" ... NOT OK"))
    sys.exit(1)

output = open(outputName, "r")
outputRead = output.readlines()

# Values
for i in range(0,len(outputRead)):
    line = outputRead[i]
    if "TOTAL ENERGY =" in line:
        testValues["HF energy"] = float(line.split()[3])

output.close()

molden1 = open(molden1Name, "r")
molden1Read = molden1.readlines()
v=0
eigenv=[]
flag=0
for i in range(0,len(molden1Read)):
    line = molden1Read[i]
    if "Ene=" in line:
        eigenv.append(float(line.split()[1]))
        v+=1
    if flag==1 and "=" not in line:
        if abs(float(line.split()[1])) >= 0.001:
            string="eigvec,"+str(v)+","+str(line.split()[0])
            testValues[string] = abs(float(line.split()[1]))
    if "[MO]" in line:
        flag=1
    if "Occup=" in line and "0.0" in line :
        testValues["e-HOMO"] = eigenv[v-2]
        break
molden1.close()

passTest = True

for value in refValues:
    diffValue = abs(refValues[value][0] - testValues[value]) 
    if ( diffValue <= refValues[value][1] ):
        passTest = passTest * True
    else :
        passTest = passTest * False
        print("%s %.8f %.8f %.2e" % ( value, refValues[value][0], testValues[value], diffValue))

if passTest :
    print(testName + str_green(" ... OK"))
else:
    print(testName + str_red(" ... NOT OK"))
    sys.exit(1)

