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
numberOfFunctionsS = int(sys.argv[5])
numberOfFunctionsP = int(sys.argv[6])
numberOfFunctionsD = int(sys.argv[7])

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


for n in range(0,3):
    out.write (str(numberOfFunctionsS+numberOfFunctionsP+numberOfFunctionsD)+"\n")
    coeficiente = 1.0000000

    interval = int(numberOfFunctionsS)
    for i in range(1,interval+1) :

            alpha = alphaS*(3**(i-1))
            out.write (  "s    1\n")
            out.write ("1  %.8f %.8f\n" % (alpha, coeficiente ) )

    interval = int(numberOfFunctionsP)
    for i in range(1,interval+1) :

            alpha = alphaP*(3**(i-1))
            out.write ( "p    1\n")
            out.write ("1  %.8f %.8f\n" % (alpha, coeficiente ) )

    interval = int(numberOfFunctionsD)
    for i in range(1,interval+1) :

            alpha = alphaD*(3**(i-1))
            out.write ( "d    1\n")
            out.write ("1  %.8f %.8f\n" % (alpha, coeficiente ) )

    out.write ("\n")

#:::::::::::::::::::::::::::::::::::
# Cerrando los archivos	
#:::::::::::::::::::::::::::::::::::

out.close()
