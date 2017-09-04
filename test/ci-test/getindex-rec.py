#!/usr/bin/python

import numpy as np

global sizearray 
sizearray = np.array([13,13,13])
array = np.array([1,1,1])

print sizearray
print "="*20

def getoneindex ( n, array ) :

    index = 0    
    for i in range(0,n) :
        ssize = 1
        for j in range(i+1,n):
            ssize = ssize * (sizearray[j]-1)
        index = index + (array[i] -1) * ssize 
    index = index + 1
    return index

def getoneindex2 ( n, array ) :

    index = 0    
    for i in range(0,array[0]) :
        for j in range(0,array[1]):
            for k in range(0,array[2]):
                print i,j,k
                index = index + 1

    return index


n = 3
c = 0

#for i in range(1,sizearray[0]):
#    for j in range(1,sizearray[1]):
#        for k in range(1,sizearray[2]):
#            if ( i == 1 or j == 1 or k == 1) :
#                array[0] = i
#                array[1] = j
#                array[2] = k
#
#                c = c + 1
#                cc = getoneindex(n,array)
#                ccc = getoneindex2(n,array)
#                print  array[0],"\t",array[1],"\t",array[2],"\t", c, "->", cc, ccc



ccc = getoneindex2(n,[1,3,1])
print ccc
