#!/usr/bin/python

import numpy as np

global string
string = np.array([[1,2,3],[10,20,30],[100,200,300]])
stringtype = np.array([[0,1,2],[0,1,2],[0,1,2]])
auxarray = np.array([0,0,0])


auxarray2 = np.array([0,1,2])
print auxarray2[0:]
print auxarray2[1:]
print auxarray2[2:]


def recursive1 ( i, n, size, inistring,count, auxarray) :
    i = i + 1
    if ( i < n ) :
        for j in range(0,size):
            newstring = str(string[i-1][j]) 
            auxstring = inistring+" "+newstring 
            auxarray[i-1] = stringtype[i-1][j]
            count, auxstring = recursive1(i,n,size,auxstring, count, auxarray)
    else :
        for j in range(0,size):
            newstring = str(string[i-1][j]) 
            auxstring = inistring+" "+newstring 
            auxarray[i-1] = stringtype[i-1][j]
            count = count + 1 
            print "tot<=2",auxstring, "->", auxarray, " = ", sum(auxarray)
            #auxarray[:] = 0

    return count, auxstring

def recursive2 ( i, n, size, inistring,count, auxarray) :
    i = i + 1
    print auxarray
    if ( i < n ) :
        if  ( sum(auxarray) <= 2 ) :
            for j in range(0,size-sum(auxarray)):
                newstring = str(string[i-1][j]) 
                auxstring = inistring+" "+newstring 
                auxarray[i-1] = stringtype[i-1][j]
                auxarray[i:] = 0
                count, auxstring = recursive2(i,n,size,auxstring, count, auxarray)
        else :
            auxstring = ""
    else :
        if  ( sum(auxarray) <= 2 ) :
            for j in range(0,size-sum(auxarray)):
                newstring = str(string[i-1][j]) 
                auxstring = inistring+" "+newstring 
                auxarray[i-1] = stringtype[i-1][j]
                auxarray[i:] = 0
                count = count + 1 
                print "tot<=2",auxstring, "->", auxarray, " = ", sum(auxarray)
        else :
            auxstring = ""

    return count, auxstring
    

nspecies = 3
size = 3
count = 0
count, finalstring = recursive1(0, nspecies, size, "", count, auxarray)

print count
    
auxarray[:] = 0
nspecies = 3
size = 3
count = 0
count, finalstring = recursive2(0, nspecies, size, "", count, auxarray)
print count

#count, finalstring = recursive (0, nspecies, size, "", count)
#print count
#count = 0

#c = 0
#
#def recursive ( i, n, size, inistring,count) :
#    i = i + 1
#    for j in range(0,size):
#        newstring = str(string[i-1][j]) 
#        auxstring = inistring+" "+newstring 
#        if ( i < n ) :
#            count, auxstring = recursive(i,n,size,auxstring, count)
#        else :
#            count = count + 1 
#            print "tot",auxstring
#
#    return count, auxstring
#
#
#print "="*20
#for j in range(0,3):
#    auxstring1 = str(string[0][j])
#    for k in range(0,3):
#        auxstring2 = auxstring1 + " "+ str(string[1][k])
#        for l in range(0,3):
#            auxstring3 = auxstring2 +" "+str(string[2][l])
#            print auxstring3
#
