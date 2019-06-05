import sys
import os
import numpy as np
import subprocess
sys.path.append( os.environ['AFQMCLAB_DIR']+"/scripts/supercubic" )
from setHoping import *

#Model Parameter
latt_n   = [4,4]
ktwist   = [0.0,0.0]
t1       = 1.0
U        = -4.0
hp       = 0.0
Nup      = 8

print Nup

#Set lattice information
latt = Latt_class( latt_n )
latt.write("latt_param")
up_i, up_j, up_K = HubbardNearestNeighborHopping(latt, ktwist, t1)
#up_i, up_j, up_K = HubbardNearestNeighborHoppingOpenOneDimension(latt, ktwist, t1, 1)

tmatrix = np.zeros( (latt.L, latt.L), dtype='complex', order='F' )
for i in range( len(up_K) ):
    tmatrix[ up_i[i], up_j[i] ] += up_K[i]

Dmatrix = np.zeros( (latt.L, latt.L), dtype='complex', order='F' )
for i in range( len(up_K) ):
    Dmatrix[ up_i[i], up_j[i] ] += -hp / ( 2.0 * np.sqrt(2.0) )

f = open("model_param", 'w')
f.write( '{:16d} \n'.format(latt.L) )
f.write( '{:16d} \n'.format(Nup) )
for i in range( latt.L ):
    for j in range( latt.L ):
        f.write( '{:26.18e} {:26.18e} \n'.format( tmatrix[j,i].real,  tmatrix[j,i].imag ) )
for i in range( latt.L ):
    for j in range( latt.L ):
        f.write( '{:26.18e} {:26.18e} \n'.format( Dmatrix[j,i].real,  Dmatrix[j,i].imag ) )
for i in range( latt.L ):
    f.write( '{:26.18e} \n'.format( U ) )
f.close()

#Method Parameter
initialType               =  "readOrderParameter"  # "setFromModel", "readWaveFunction", "readOrderParameter"
convergeTolerance         =  1e-5
maxIterateStep            =  100
annealMagnitude           =  0.3
annealStep                =  10
relaxMagnitude            =  0.6                 # 1.0 fully relax to new order paramter, 0.0 not update
initMu                    =  0.3
deltaMu                   =  1.0
Neta                      =  1e-6
seed                      =  985456376             # 1. read file, 0. random, else is seeds
initOrderParamType        =  "random"              # zero random

#write method_param
f = open('bcs_param', 'w')
f.write(   '{:>26} \n'.format(initialType       ) )
f.write('{:26.18e} \n'.format(convergeTolerance ) )
f.write(   '{:26d} \n'.format(maxIterateStep    ) )
f.write('{:26.18e} \n'.format(annealMagnitude   ) )
f.write(   '{:26d} \n'.format(annealStep        ) )
f.write('{:26.18e} \n'.format(relaxMagnitude    ) )
f.write('{:26.18e} \n'.format(initMu            ) )
f.write('{:26.18e} \n'.format(deltaMu           ) )
f.write('{:26.18e} \n'.format(Neta              ) )
f.write(   '{:26d} \n'.format(seed              ) )
f.close()

#Run script to generate orderParameter file.
pythonFileDir =  os.path.dirname(os.path.realpath(__file__))
if initOrderParamType is not None:
    subprocess.call( "python {0}/initialOrderParameter.py {1}".format(pythonFileDir, initOrderParamType),shell=True)
