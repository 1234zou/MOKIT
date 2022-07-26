from pyscf import gto
from py2fch_direct import mol2fch
from py2fch import py2fch

mol = gto.M(atom='O 0.0 0.0 0.1; H 0.0 0.0 1.0', basis='cc-pvtz', charge=-1).build()
mf = mol.RHF().run()

#S = mol.intor('int1e_ovlp')
#print(S)
mol2fch(mol, 'test.fch')
py2fch('test.fch', mol.nao, mol.nao, mf.mo_coeff, 'a', mf.mo_energy, 0)
