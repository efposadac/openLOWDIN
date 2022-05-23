#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import subprocess
import numpy as np
import time
import signal

#:::::::::::::::::::::::::::::::::::
# Cargar en una lista el archivo
#:::::::::::::::::::::::::::::::::::

basisName = sys.argv[1]

aS = float(sys.argv[2])
aP = float(sys.argv[3])
aD = float(sys.argv[4])
aF = float(sys.argv[5])
aG = float(sys.argv[6])
aH = float(sys.argv[7])
nS = int(sys.argv[8])
nP = int(sys.argv[9])
nD = int(sys.argv[10])
nF = int(sys.argv[11])
nG = int(sys.argv[12])
nH = int(sys.argv[13])

out = open (basisName,"w")

##:::::::::::::::::::::::::::::::::::
# Diccionario de letras de momento angular
#:::::::::::::::::::::::::::::::::::

angularMomentDic = {
'S':'0',
'P':'1',
'D':'2',
'F':'3',
'G':'4',
'H':'5',
's':'0',
'p':'1',
'd':'2',
'f':'3',
'g':'4',
'h':'5'
}


#:::::::::::::::::::::::::::::::::::
# Saludo :D
#:::::::::::::::::::::::::::::::::::

#print (" Script para crear una base even-tempered en formato de Lowdin1")
#print (" alpha_{i+1} = alpha_{i}*beta")
#print (" Nombre de la base: " + basisName )
#print (" particleName: " + particleName )
#print (" particlSymbol: " + particleSymbol )
#print (" basisType: " + basisType )
#print (" alpha: " + str(alpha) )
#print (" beta: " + str(beta) )
#print (" numberOfFunctions: " + numberOfFunctions )
#print (" angularMoment: " + angularMoment )

#:::::::::::::::::::::::::::::::::::
# MethodName in TASK
#:::::::::::::::::::::::::::::::::::

headerSpecies = list()
headerSpecies.append('O-HYDROGEN H ('+basisName+') BASIS TYPE: 1\n')
headerSpecies.append('O-EALONE EA- ('+basisName+') BASIS TYPE: 1\n')
headerSpecies.append('O-POSITRON E+ ('+basisName+') BASIS TYPE: 1\n')	

coreBasisSet = ("""1 0 3
33.87000000 0.00606800
5.09500000 0.04530800
1.15900000 0.20282200
2 0 1
0.32580000 1.00000000
3 0 1
0.10270000 1.00000000
4 1 1
1.40700000 1.00000000
5 1 1
0.38800000 1.00000000
6 2 1
1.05700000 1.00000000
""")

coreBasisSet = ""

coeficiente = 1.0000000
nCore = 0

for species in range(0,len(headerSpecies)) :

    out.write (headerSpecies[species])	
    out.write ('#\n')
    out.write (str(nS+nP+nD+nF+nG+nH+nCore)+"\n")
    out.write (coreBasisSet)
    
    counter = nCore
    interval = int(nS)
    for i in range(0,interval+0) :
    
    	counter = counter +1
    	alpha = aS*(3.16227766**(-i))
    	#alpha = aS*(3.16227766**(i-1))
    	out.write (str(counter)+" 0 1\n")
    	out.write ("%.8f %.8f\n" % (alpha, coeficiente ) )
    
    interval = int(nP)
    for i in range(0,interval+0) :
    
    	counter = counter +1
    	alpha = aP*(3.16227766**(-i))
    	#alpha = aP*(3.16227766**(i-1))
    	out.write (str(counter)+" 1 1\n")
    	out.write ("%.8f %.8f\n" % (alpha, coeficiente ) )
    
    interval = int(nD)
    for i in range(0,interval+0) :
    
    	counter = counter +1
    	alpha = aD*(3.16227766**(-i))
    	#alpha = aD*(3.16227766**(i-1))
    	out.write (str(counter)+" 2 1\n")
    	out.write ("%.8f %.8f\n" % (alpha, coeficiente ) )
    
    interval = int(nF)
    for i in range(0,interval+0) :
    
    	counter = counter +1
    	alpha = aF*(3.16227766**(-i))
    	#alpha = aD*(3.16227766**(i-1))
    	out.write (str(counter)+" 3 1\n")
    	out.write ("%.8f %.8f\n" % (alpha, coeficiente ) )
    
    interval = int(nG)
    for i in range(0,interval+0) :
    
    	counter = counter +1
    	alpha = aG*(3.16227766**(-i))
    	#alpha = aD*(3.16227766**(i-1))
    	out.write (str(counter)+" 4 1\n")
    	out.write ("%.8f %.8f\n" % (alpha, coeficiente ) )
    
    interval = int(nH)
    for i in range(0,interval+0) :
    
    	counter = counter +1
    	alpha = aH*(3.16227766**(-i))
    	#alpha = aD*(3.16227766**(i-1))
    	out.write (str(counter)+" 5 1\n")
    	out.write ("%.8f %.8f\n" % (alpha, coeficiente ) )

    out.write ('\n')





#:::::::::::::::::::::::::::::::::::
# Cerrando los archivos	
#:::::::::::::::::::::::::::::::::::

out.close()
