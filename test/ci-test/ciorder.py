#!/usr/bin/python

import numpy as np

sizearray = np.array([13,13,13])
ciorder = np.array([2,2,1])

maxci = sum(ciorder)
print maxci
for ci in range(1,maxci+1) :
    c = 0
    for i in range(0,ciorder[0]+1):
        for j in range(0,ciorder[1]+1):
            for k in range(0,ciorder[2]+1):
                totalorder = i + j + k
                if ( totalorder <= ci ) : 
                    c = c + 1
                    print i,j,k, "\t", totalorder, "|", (i)*3*2+(j)*2+k+1

    print "total ",c 
