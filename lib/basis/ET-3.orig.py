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
nS = int(sys.argv[7])
nP = int(sys.argv[8])
nD = int(sys.argv[9])
nF = int(sys.argv[10])
nG = int(sys.argv[11])

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

#out.write ('O-EALONE EA- ('+basisName+') BASIS TYPE: 1\n')	
out.write ('O-HYDROGEN H ('+basisName+') BASIS TYPE: 1\n')	
out.write ('#\n')
out.write (str(nS+nP+nD+nF+nG)+"\n")
coeficiente = 1.0000000

counter = 0
interval = int(nS)
for i in range(0,interval+0) :

	counter = counter +1
	alpha = aS*(3.162**(-i))
	#alpha = aS*(3.162**(i-1))
	out.write (str(counter)+" 0 1\n")
	out.write ("%.8f %.8f\n" % (alpha, coeficiente ) )

interval = int(nP)
for i in range(0,interval+0) :

	counter = counter +1
	alpha = aP*(3.162**(-i))
	#alpha = aP*(3.162**(i-1))
	out.write (str(counter)+" 1 1\n")
	out.write ("%.8f %.8f\n" % (alpha, coeficiente ) )

interval = int(nD)
for i in range(0,interval+0) :

	counter = counter +1
	alpha = aD*(3.162**(-i))
	#alpha = aD*(3.162**(i-1))
	out.write (str(counter)+" 2 1\n")
	out.write ("%.8f %.8f\n" % (alpha, coeficiente ) )

interval = int(nF)
for i in range(0,interval+0) :

	counter = counter +1
	alpha = aF*(3.162**(-i))
	#alpha = aD*(3.162**(i-1))
	out.write (str(counter)+" 3 1\n")
	out.write ("%.8f %.8f\n" % (alpha, coeficiente ) )

interval = int(nG)
for i in range(0,interval+0) :

	counter = counter +1
	alpha = aG*(3.162**(-i))
	#alpha = aD*(3.162**(i-1))
	out.write (str(counter)+" 4 1\n")
	out.write ("%.8f %.8f\n" % (alpha, coeficiente ) )

out.write ('\n')
out.write ('O-POSITRON E+ ('+basisName+') BASIS TYPE: 1\n')	
out.write ('#\n')
out.write (str(nS+nP+nD+nF+nG)+"\n")

coeficiente = 1.0000000

counter = 0
interval = int(nS)
for i in range(0,interval+0) :

	counter = counter +1
	alpha = aS*(3.162**(-i))
	#alpha = aS*(3.162**(i-1))
	out.write (str(counter)+" 0 1\n")
	out.write ("%.8f %.8f\n" % (alpha, coeficiente ) )

interval = int(nP)
for i in range(0,interval+0) :

	counter = counter +1
	alpha = aP*(3.162**(-i))
	#alpha = aP*(3.162**(i-1))
	out.write (str(counter)+" 1 1\n")
	out.write ("%.8f %.8f\n" % (alpha, coeficiente ) )

interval = int(nD)
for i in range(0,interval+0) :

	counter = counter +1
	alpha = aD*(3.162**(-i))
	#alpha = aD*(3.162**(i-1))
	out.write (str(counter)+" 2 1\n")
	out.write ("%.8f %.8f\n" % (alpha, coeficiente ) )

interval = int(nF)
for i in range(0,interval+0) :

	counter = counter +1
	alpha = aF*(3.162**(-i))
	#alpha = aD*(3.162**(i-1))
	out.write (str(counter)+" 3 1\n")
	out.write ("%.8f %.8f\n" % (alpha, coeficiente ) )

interval = int(nG)
for i in range(0,interval+0) :

	counter = counter +1
	alpha = aG*(3.162**(-i))
	#alpha = aD*(3.162**(i-1))
	out.write (str(counter)+" 4 1\n")
	out.write ("%.8f %.8f\n" % (alpha, coeficiente ) )



#:::::::::::::::::::::::::::::::::::
# Cerrando los archivos	
#:::::::::::::::::::::::::::::::::::

out.close()
