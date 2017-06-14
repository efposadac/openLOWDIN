#!/usr/bin/python

import numpy as np

global sizearray 
sizearray = np.array([4,4,4])
array = np.array([1,1,1])

print sizearray
print "="*20
c = 0
for i in range(1,sizearray[0]):
    for j in range(1,sizearray[1]):
        for k in range(1,sizearray[2]):
            c = c + 1
            auxi = (float(c) / float((sizearray[2]-1)*(sizearray[1]-1)))

            auxii = auxi - int(auxi)
            auxjj = auxii* ((sizearray[2]-1)*(sizearray[1]-1)) #new c
            ii = int(auxi) + 1

            auxi2 = (float(auxjj) /float ((sizearray[2]-1)))
            auxii2 = auxi2 - int(auxi2) 
            auxjj2 = auxii2* ((sizearray[2]-1))
            ii2 = int(auxi2) + 1

            print i,j,k,c, k+(j-1)*(sizearray[2]-1)+(i-1)*(sizearray[2]-1)*(sizearray[1]-1),ii,ii2,int(auxjj2)

def getoneindex ( n, array ) :

    index = 0    
    for i in range(0,n) :
        ssize = 1
        for j in range(i+1,n):
            ssize = ssize * (sizearray[j]-1)
        index = index + (array[i] -1) * ssize 
    index = index + 1
    return index

n = 3
c = 0

for i in range(1,sizearray[0]):
    for j in range(1,sizearray[1]):
        for k in range(1,sizearray[2]):
            array[0] = i
            array[1] = j
            array[2] = k
            

            c = c + 1
            cc = getoneindex(n,array)
            print  array, c, "->", cc

def expandindex ( i, n, index, array ) :
    i = i + 1
    if ( i < n ) :
        ssize = 1
        for j in range(i,n) :
            ssize = ssize * (sizearray[j]-1)
        print ssize
        auxindex = float(index)/float(ssize)
        newindex = (auxindex - int(auxindex))*ssize
        array[i-1] = int(auxindex) + 1
        array = expandindex ( i, n, newindex, array )
    else: 
        print "a",array
        print "index", index
        ssize = (sizearray[i-1]-1)
        auxindex = float(index)/float(ssize)
        array[i-1] = int(  (auxindex - int(auxindex)) * (sizearray[i-1]-1)  )

    return array

#            auxi2 = (float(auxjj) /float ((sizearray[2]-1)))
# auxii2 = auxi2 - int(auxi2) 
 #           auxjj2 = auxii2* ((sizearray[2]-1))



n = 3
index = 16
i = 0
print expandindex ( i, n, index, array )
