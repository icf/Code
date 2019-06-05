import sys
import os
sys.path.append( os.environ['AFQMCLAB_DIR']+"/scripts/supercubic" )
from setHoping import *

#Model Parameter
latt_n   = [4,4]
ktwist   = [0.0,0.0]
t1       = 1.0
t2       = 0.0
pin      = 0.0
U        = -4.0
mu       = 0.0
Ntot     = 16
UpDnFlag = 0    # 0 up=dn, 1 up=conj(dn) ==> different twist

#---------------------------------------------------------------------------------------
#Set lattice information
latt = Latt_class( latt_n )
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
dt                        = 0.01
thermalSize               = 100
writeNumber               = 10
measureNumberPerWrite     = 2
measureSkipStep           = 5
walkerSizePerThread       = 50
backPropagationStep       = 10
decompType                = "densitySpin"  # "densityCharge", "densitySpin", "hopCharge", "hopSpin"
forceType                 = "dynamicForce"   # "dynamicForce", "constForce"
forceCap                  = 1.5
initialPhiTFlag           = "setFromModel"   #"setFromModel", "setRandomly", "readFromFile": phiT.dat
initialSCPhiTFlag         = "setFromGHF_readOccupancy"   #"setFromDensity_Analytical", "setFromDensity_VMGpBCS","setFromDensity_VMGpBCS_withGHF_orbital","setFromGHF","setFromGHF_readOccupancy","setFromDensity_VMGpBCS_withInput_orbital": only work in pBCS->GpBCS form 
scSteps                   = 2  
initialWalkerFlag         = "setFromModel"   #"setFromModel", "setRandomly", "sampleFromPhiT \\only work for hao's style: still may has some problems","readFromFile","readAllWalkers"
mgsStep                   = 10
popControlStep            = 10
ET                        = -10
ETAdjustStep              = 10
ETAdjustMaxSize           = 100
seed                      = 985456376        # -1. read file, 0. random, else is seeds

#write method_param
f = open('afqmc_param', 'w')
f.write(" {:<36} {:<26.18e} \n".format("dt", dt ) )
f.write(" {:<36} {:<26} \n".format("thermalSize", thermalSize) )
f.write(" {:<36} {:<26} \n".format("writeNumber", writeNumber) )
f.write(" {:<36} {:<26} \n".format("measureNumberPerWrite", measureNumberPerWrite) )
f.write(" {:<36} {:<26} \n".format("measureSkipStep", measureSkipStep) )
f.write(" {:<36} {:<26} \n".format("walkerSizePerThread",walkerSizePerThread) )
f.write(" {:<36} {:<26} \n".format("backPropagationStep",backPropagationStep) )
f.write(" {:<36} {:<26} \n".format("decompType",decompType) )
f.write(" {:<36} {:<26} \n".format("forceType", forceType) )
f.write(" {:<36} {:<26.18e} \n".format("forceCap", forceCap) )
f.write(" {:<36} {:<26} \n".format("initialPhiTFlag",initialPhiTFlag) )
f.write(" {:<36} {:<26} \n".format("initialSCPhiTFlag",initialSCPhiTFlag) )
f.write(" {:<36} {:<26} \n".format("scSteps",scSteps) )
f.write(" {:<36} {:<26} \n".format("initialWalkerFlag",initialWalkerFlag) )
f.write(" {:<36} {:<26} \n".format("mgsStep", mgsStep) )
f.write(" {:<36} {:<26} \n".format("popControlStep", popControlStep) )
f.write(" {:<36} {:<26.18e} \n".format("ET", ET) )
f.write(" {:<36} {:<26} \n".format("ETAdjustStep", ETAdjustStep) )
f.write(" {:<36} {:<26} \n".format("ETAdjustMaxSize", ETAdjustMaxSize) )
f.write(" {:<36} {:<26} \n".format("seed", seed) )
f.close()

#write constForce_param
f = open('constForce_param', 'w')
for i in range( latt.L ):
    f.write( '{:26.18e} \n'.format( 0.0 ) )
f.close()
