#!/usr/bin/python

import numpy as np

global string
string = np.array([[1,2],[10,20,30],[100,200,300]])
print string

n = 0
for i in xrange(0,2):
    for j in xrange(0,3):
        for k in xrange(0,3):
            n = n + 1
            print n, string[0][i], string[1][j], string[2][k]

nspecies = 3
indexConf = [0,0,-1]
nsize = [2,3,3]

auxindexconf = [0,0,0]

totalsize = 1
for i in xrange(0,nspecies):
    totalsize = totalsize*nsize[i]
print "totalsize", totalsize

for x in xrange(0,totalsize):
    s = nspecies-1
    indexConf[s] = indexConf[s] + 1
    for i in xrange(nspecies-1,0,-1):
        auxsize = 1
        # get size 
        for j in xrange(i,nspecies):
            auxsize = auxsize*nsize[j]
        print i
        if auxindexconf[i] == auxsize:

            for j in xrange(i,nspecies):
                indexConf[j] = 0

            auxindexconf[i] = 0
            indexConf[i-1] = indexConf[i-1] + 1

        auxindexconf[i] = auxindexconf[i] + 1
    #print indexConf
    conf = ""
    for i in xrange(0,nspecies):
        conf = conf + " " + str(string[i][indexConf[i]])
    print conf

        
