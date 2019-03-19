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

alphaS = float(sys.argv[2])
alphaP = float(sys.argv[3])
alphaD = float(sys.argv[4])
numberOfFunctions = int(sys.argv[5])

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
's':'0',
'p':'1',
'd':'2',
'f':'3',
'g':'4'
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

out.write ('O-POSITRON E+ ('+basisName+') BASIS TYPE: 2\n')	
out.write ('#\n')
out.write (str(numberOfFunctions*3)+"\n")
coeficiente = 1.0000000

counter = 0
interval = int(numberOfFunctions)
for i in range(1,interval+1) :

	counter = counter +1
	alpha = alphaS*(3**(i-1))
	out.write (str(counter)+" 0 1\n")
	out.write ("%.8f %.8f\n" % (alpha, coeficiente ) )

for i in range(1,interval+1) :

	counter = counter +1
	alpha = alphaP*(3**(i-1))
	out.write (str(counter)+" 1 1\n")
	out.write ("%.8f %.8f\n" % (alpha, coeficiente ) )

for i in range(1,interval+1) :

	counter = counter +1
	alpha = alphaD*(3**(i-1))
	out.write (str(counter)+" 2 1\n")
	out.write ("%.8f %.8f\n" % (alpha, coeficiente ) )

#:::::::::::::::::::::::::::::::::::
# Cerrando los archivos	
#:::::::::::::::::::::::::::::::::::

out.close()
