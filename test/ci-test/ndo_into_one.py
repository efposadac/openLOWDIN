#!/usr/bin/python

import numpy as np
import time

ni = 5
nj = 5

start1 = time.time()
n = 0
a = 0 
for i in xrange(1,ni+1):
    for j in xrange(1,nj+1):
        a = a + 1
#        print a,i,j
end1 = time.time()

start2 = time.time()
i = 1
j = 0
for x in xrange(1,ni*nj+1):
    j = j + 1
#    print i,j
    if j == nj :
        j = 0
        i = i + 1
end2 = time.time()

start3 = time.time()
for x in xrange(1,ni*nj+1):
    i = (x-1)/nj+1
    j =(x-1)%nj+1
    print x, i, j, i*j
end3 = time.time()

print end1-start1
print end2-start2
print end3-start3
