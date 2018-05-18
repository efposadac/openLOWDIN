#!/usr/bin/env python

import numpy
import sys

inputName = sys.argv[1]
n = int(sys.argv[2]) + 1
inputFile = open(inputName, "r" )

inputLines = inputFile.readlines()

omptime = numpy.zeros(n)
omptime_n = 0
omptime_average = numpy.zeros(n)

for i in xrange(0,len(inputLines)) :
    if "omptime" in inputLines[i] :
        omptime_n = omptime_n + 1
        for j in xrange(1,n ) :
            omptime[j] = omptime[j] + float(inputLines[i+j].split()[1])
            #print j, omptime_n, float(inputLines[i+j].split()[1]), omptime[j]

print "-------"
print " Total / Average", omptime_n
for i in xrange(1,n):
    omptime_average[i] = omptime[i]/omptime_n
    print omptime[i], omptime_average[i] 

         




