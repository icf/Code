import sys
import os
sys.path.append( os.environ['AFQMCLAB_DIR']+"/scripts/supercubic" )
from lattClass import *

if len(sys.argv) < 2:
    sys.exit("Missing arguments!!! Example: python initialOrderParameter.py random")
typ = sys.argv[1]

latt  = Latt_class.read("latt_param")

if typ=="zero":
    pairing = np.zeros(latt.L, dtype=complex)

elif typ=="random":
    pairing = np.zeros(latt.L, dtype=complex)
    for i in range(latt.L):
        pairing[i] = np.random.random()

np.savetxt("pairing.dat", np.column_stack(( np.real(pairing), np.imag(pairing) )),fmt=('%26.16e', '%26.16e'), header="{:26d} \n{:26d}".format(1, latt.L), comments='')