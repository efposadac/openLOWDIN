#!/usr/bin/python

import numpy as np

global sizearray 
sizearray = np.array([13,13,13])
array = np.array([1,1,1])

print sizearray
print "="*20

c = 0
for i in range(1,sizearray[0]):
    for j in range(1,sizearray[1]):
        for k in range(1,sizearray[2]):
            c = c + 1

cc = 1
for s in range(0,2+1):
    c = 0
    for i in range(1,sizearray[s]):
        c = c + cc
        print c
    cc = c
    print "cc",cc
