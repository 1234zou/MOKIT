from pyscf import gto
from py2fch_direct import mol2fch

mol = gto.M(atom='O 0.0 0.0 0.1; H 0.0 0.0 1.0', basis='cc-pvtz', charge=-1).build()

mol2fch(mol)
