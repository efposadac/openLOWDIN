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

nombreInputxml = sys.argv[1]

inputxml = open (nombreInputxml,"r")
inputxmlRead = inputxml.readlines()

out = open (nombreInputxml[:-4],"w")

#:::::::::::::::::::::::::::::::::::
# Funciones para manejo de mayusculas y minusculas
#:::::::::::::::::::::::::::::::::::

def Minuscula(String):
	String = String.lower()
	return String

def Mayuscula(String):
	String = String.upper()
	return String

##:::::::::::::::::::::::::::::::::::
# Diccionario de letras de momento angular
#:::::::::::::::::::::::::::::::::::

angularMoment = {
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

print (" Script para convertir una base de xml a texto plano")
print (" Nombre de la entrada: "+nombreInputxml)
print (" Nombre de la salida: "+nombreInputxml[:-4])

#:::::::::::::::::::::::::::::::::::
# MethodName in TASK
#:::::::::::::::::::::::::::::::::::
numberContractions = 0
size = 0
numPrim = 0 

numLine = 0

for line in inputxmlRead :
	numLine = numLine + 1                ## Un contador para que mas adelante no sea necesario leer toda la lista de nuevo
	line = line.replace(' =','=')        ## Para evitar problemas...
	line = line.replace('= ','=')

	auxline = line.split()               ## Separar la linea del archivo por espacios en blanco
	auxline = auxline[1:]                ## Quitar el primer campo. No es necesario

	if "<!--!" in line and "-->" :
		line = line.replace("<!--!","")
		line = line.replace("-->","")
		out.write(line)


	if "<Basis" in line:
		for data in auxline :			## Leer la nueva lista que surge al dividir una linea por espacios en blanco
			data = data.replace('>','')		## Quitando caracteres no deseados 
			if "name" in data :
				basisName = data.split('"')[1]	## Buscar cada campo que se desea sacar
			if "type" in data :
				basisType = data.split('"')[1]

	if "<Atom" in line:

		for data in auxline :	
			data = data.replace('>','')
			if "name" in data :
				name = data.split('"')[1]
				numFunction = 0 		## Iterador del numero de contracciones por átomo. Se reinicia en cada átomo
			if "symbol" in data :
				symbol = data.split('"')[1]
			if "numberContractions" in data :
				numberContractions = int(data.split('"')[1])
			if "contraction" in data :
				contraction = data.split('"')[1]

		out.write ("O-"+Mayuscula(name)+" "+Mayuscula(symbol)+" ("+Mayuscula(basisName)+") BASIS TYPE: "+basisType+"\n")  # Escribir
		out.write ("# "+contraction+"\n")
		out.write (str(numberContractions)+"\n")

		for line2 in inputxmlRead[numLine:] :		## Leer lineas del archivo desde que encontramos atom hasta que cierra atom
			auxline = line2.split()			
			auxline = auxline[1:]
				
			if "<ContractedGaussian" in line2:
				numFunction = numFunction + 1 	## Arranca el contador para indicar en que función se está leyendo 

				for data in auxline :	
					data = data.replace('>','')

					if "size" in data and numberContractions > 0:  	## Debe haber mas de una contracción para buscarlas...
						size = int(data.split('"')[1])
						numPrim = 0				## Reinicia el contador para el # de primitivas
					if "shellCode" in data :
						shellCode = data.split('"')[1]
				exponent = np.zeros((size),dtype=float)			## Creando arreglos
				coefficient = np.zeros((size),dtype=float)		
				out.write(str(numFunction)+" "+angularMoment[shellCode]+" "+str(size)+"\n")

			if "PrimitiveGaussian" in line2 :
				for data in auxline :	
					data = data.replace('>','')

					if "exponent" in data and size > 0:
						numPrim = numPrim + 1			## Contando el de primitiva de la contracción
						exponent[numPrim-1] = float(data.split('"')[1])
					if "coefficient" in data and size > 0:
						coefficient[numPrim-1] = float(data.split('"')[1])
				if numPrim == size :					## Acabamos de leer las primitivas
					for i in range(size) :
						out.write("%.8f"%exponent[i] +" "+"%.8f"%coefficient[i]+"\n")	## Escribir exp y coeff

			if "</Atom" in line2:
				out.write("\n")
				break		## Llegamos al fin del átomo, Se puede parar la segunda lectura del archivo

#:::::::::::::::::::::::::::::::::::
# Cerrando los archivos	
#:::::::::::::::::::::::::::::::::::

inputxml.close()
out.close()
