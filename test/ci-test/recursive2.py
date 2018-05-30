#!/usr/bin/python

import numpy as np

global string
string = np.array([[1,2,3],[10,20,30],[100,200,300]])
print string

n = 0
for i in xrange(0,3):
    for j in xrange(0,3):
        for k in xrange(0,3):
            n = n + 1
            print n, string[0][i], string[1][j], string[2][k]

indexConf = [0,0,-1]
ni = 3
nj = 3
nk = 3

ii = 0
jj = 0 
kk = 0

for x in xrange(0,ni*nj*nk):
    s = 2
    indexConf[s] = indexConf[s] + 1
    if kk == nk :
        indexConf[s] = 0
        indexConf[s-1] = indexConf[s-1] + 1
        kk = 0
    if jj == nj*nk :
        indexConf[s] = 0
        indexConf[s-1] = 0
        indexConf[s-2] = indexConf[s-2] + 1
        jj = 0
    kk = kk + 1
    jj = jj + 1
    #print indexConf
    print string[0][indexConf[0]], string[1][indexConf[1]], string[2][indexConf[2]]
        
