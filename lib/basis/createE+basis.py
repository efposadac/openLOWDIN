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

nombreInputEl = sys.argv[1]

inputEl = open (nombreInputEl,"r")
inputElRead = inputEl.readlines()

#:::::::::::::::::::::::::::::::::::
# Saludo :D
#:::::::::::::::::::::::::::::::::::

print (" Script para generar una base positronica a partir de una electronica")
print (" Nombre de la entrada: "+nombreInputEl)
print (" Nombre de la salida: E+-Z-"+nombreInputEl)


#:::::::::::::::::::::::::::::::::::
# El trabajo sucio
#:::::::::::::::::::::::::::::::::::
atom="trololo"
out=open(atom,"w")

for line in inputElRead :
	auxline = line.split()               ## Separar la linea del archivo por espacios en blanco
        if len(auxline) > 0:
                if not("#" in auxline[0]) :
                        if "O-" in auxline[0] :
                                atom=auxline[1]
                                out.close()
                                out=open("E+-"+atom+"-"+nombreInputEl,"w")
                                out.write("O-POSITRON E+ ("+atom+"-"+nombreInputEl+") BASIS TYPE:2 \n")
                        
                if atom != "trololo" and not("O-" in auxline[0]):
                        for i in range( 0, len(auxline)):
                                out.write(str(auxline[i])+" ")
                                if(i == len(auxline)-1): out.write("\n")


out.close()
                        
