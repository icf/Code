import sys
import os
import subprocess
sys.path.append( os.environ['AFQMCLAB_DIR']+"/scripts/supercubic" )
from setHoping import *

#Model Parameter
latt_n   = [4,4]
ktwist   = [0.0,0.0]
t1       = 1.0
t2       = 0.1
pin      = 0.25
U        = 12.0
mu       = 0.0
Ntot     = 14
UpDnFlag = 0    # 0 up=dn, 1 up=conj(dn) ==> different twist

#Set lattice information
latt = Latt_class( latt_n )
latt.write("latt_param")

up_i, up_j, up_K = Hubbard2DNextNearestNeighborHoppingOpenY(latt, ktwist, t2)
up_i2, up_j2, up_K2 = HubbardNearestNeighborHopping(latt, ktwist, t1)

#up_i, up_j, up_K = HubbardNearestNeighborHoppingOpenOneDimension(latt, ktwist, t1, openD=0)
if UpDnFlag == 0:
    dn_i = up_i; dn_j = up_j; dn_K = up_K
    dn_i2 = up_i2; dn_j2 = up_j2; dn_K2 = up_K2
elif UpDnFlag == 1:
    dn_i, dn_j, dn_K = HubbardNearestNeighborHopping(latt, -np.array(ktwist), t1)
#    up_i, up_j, up_K = HubbardNearestNeighborHoppingOpenOneDimension(latt, ktwist, t1, openD=0)
else:
    print( "WRONG!!! Do not know UpDnFlag!!!" )
    sys.exit(1)

up_i3, up_j3, up_K3 = Hubbard2DRowAFMPinning(latt, 0, pin)
up_i4, up_j4, up_K4 = Hubbard2DRowAFMPinning(latt, latt_n[1]-1, -1*pin)

dn_i3, dn_j3, dn_K3 = Hubbard2DRowAFMPinning(latt, 0, -1*pin)
dn_i4, dn_j4, dn_K4 = Hubbard2DRowAFMPinning(latt, latt_n[1]-1, pin)
#---------------------------------------------------------------------------------------

Kmatrix = np.zeros( (2*latt.L, 2*latt.L), dtype='complex', order='F' )
for i in range( len(up_K) ):
    Kmatrix[ up_i[i], up_j[i] ] += up_K[i]
for i in range( len(up_K2) ):
    Kmatrix[ up_i2[i], up_j2[i] ] += up_K2[i]
for i in range( len(up_K3) ):
    Kmatrix[ up_i3[i], up_j3[i] ] += up_K3[i]
for i in range( len(up_K4) ):
    Kmatrix[ up_i4[i], up_j4[i] ] += up_K4[i]
for i in range( len(dn_K) ):
    Kmatrix[ dn_i[i]+latt.L, dn_j[i]+latt.L ] += dn_K[i]
for i in range( len(dn_K2) ):
    Kmatrix[ dn_i2[i]+latt.L, dn_j2[i]+latt.L ] += dn_K2[i]
for i in range( len(dn_K3) ):
    Kmatrix[ dn_i3[i]+latt.L, dn_j3[i]+latt.L ] += dn_K3[i]
for i in range( len(dn_K4) ):
    Kmatrix[ dn_i4[i]+latt.L, dn_j4[i]+latt.L ] += dn_K4[i]

f = open("model_param", 'w')
f.write( '{:16d} \n'.format(latt.L) )
f.write( '{:16d} \n'.format(Ntot) )
for i in range( 2*latt.L ):
    for j in range( 2*latt.L ):
       f.write( '{:26.18e} {:26.18e} \n'.format( Kmatrix[j,i].real,  Kmatrix[j,i].imag ) )
for i in range( latt.L ):
    f.write( '{:26.18e} \n'.format( mu ) )
for i in range( latt.L ):
    f.write( '{:26.18e} \n'.format( 0.0 ) )  #hx
for i in range( latt.L ):
    f.write( '{:26.18e} \n'.format( 0.0 ) )  #hy
for i in range( latt.L ):
    f.write( '{:26.18e} \n'.format( 0.0 ) )  #hz
for i in range( latt.L ):
    f.write( '{:26.18e} \n'.format( U ) )
f.close()


#Method Parameter
initialType               =  "readOrderParameter"  # "setFromModel", "readWaveFunction", "readOrderParameter"
convergeType              =  "energy"        # "energy", "orderParameter"
convergeTolerance         =  1e-10
maxIterateStep            =  100
annealMagnitude           =  0.3
annealStep                =  0
relaxMagnitude            =  1.0          # 1.0 fully relax to new order paramter, 0.0 not update
seed                      =  985456376    # 1. read file, 0. random, else is seeds
initOrderParamType        =  "xafm"       # None, zafm, xafm, yafm, random

#write method_param
f = open('ghf_param', 'w')
f.write(   '{:>26} \n'.format(initialType       ) )
f.write(   '{:>26} \n'.format(convergeType      ) )
f.write('{:26.18e} \n'.format(convergeTolerance ) )
f.write(   '{:26d} \n'.format(maxIterateStep    ) )
f.write('{:26.18e} \n'.format(annealMagnitude   ) )
f.write(   '{:26d} \n'.format(annealStep        ) )
f.write('{:26.18e} \n'.format(relaxMagnitude    ) )
f.write(   '{:26d} \n'.format(seed              ) )
f.close()

#Run script to generate orderParameter file.
pythonFileDir =  os.path.dirname(os.path.realpath(__file__))
if initOrderParamType is not None:
    subprocess.call( "python {0}/initialOrderParameter.py {1}".format(pythonFileDir, initOrderParamType),shell=True)
