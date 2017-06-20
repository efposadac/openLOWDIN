#!/usr/bin/python

import numpy as np

global occ
global basis
occ = 5
basis = 10
    

for i in range (1,occ+1) :
    for j in range (i+1,occ+1) :
        for a in range (occ+1,basis+1) :
            for b in range (a +1,basis+1) :
                print i, j, "|", a,b
print "-"*20
for i in range (1,occ+1) :
    for a in range (occ+1,basis+1) :
        for j in range (i+1,occ+1) :
            for b in range (a +1,basis+1) :
                print i, j, "|", a,b

def configuration ( occupied, virtual, ci, level ):

    ci = ci + 1
    if ( ci == 1 and ci < level ) : #first
        for i in range (occupied[ci-1]+1,occ+1) :
            for a in range (virtual[ci-1]+1,basis+1) :
                occupied[ci-1] = i 
                virtual[ci-1] = a 
                configuration ( occupied, virtual, ci, level)
            virtual[:] = occ

    if ( ci > 1 and ci < level ) : #mid
        for i in range (occupied[ci-2]+1,occ+1) :
            for a in range (virtual[ci-2]+1,basis+1) :
                occupied[ci-1] = i 
                virtual[ci-1] = a 
                configuration ( occupied, virtual, ci, level)
            #virtual[:] = occ

    else :  #final
        for i in range (occupied[ci-2]+1,occ+1) :
            for a in range (virtual[ci-2]+1,basis+1) :
                occupied[ci-1] = i 
                virtual[ci-1] = a 
                print "\t"*int(ci+2),occupied, virtual

            if ( ci == 1 ) : # only singles
                virtual[:] = occ

    return

print "-"*20

global occupiedindex 
global virtualindex 
level = 5

for i in range(1,level+1):
    print i
    occupiedindex = np.zeros((i),dtype=np.int)
    virtualindex = np.zeros((i),dtype=np.int)
    virtualindex[:] = occ 
    print virtualindex
    if ( i <= occ ):
        configuration ( occupiedindex, virtualindex, 0, i )

