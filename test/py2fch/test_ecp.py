from pyscf import gto
from mokit.lib.py2fch_direct import fchk

mol = gto.M(atom='Cu 0.0 0.0 0.1; Au 0.0 0.0 3.1', 
        basis= 'lanl2dz', ecp = 'lanl2dz',
        spin=2).build()
mf = mol.UHF()
mf.conv_tol *= 0.1
mf.run()

fchk(mf, 'testcu.fch')
#mol2fch(mol, 'testcu.fch', uhf=True)
#py2fch('testcu.fch', mol.nao, mol.nao, mf.mo_coeff[0], 'a', mf.mo_energy[0], False, False)
#py2fch('testcu.fch', mol.nao, mol.nao, mf.mo_coeff[1], 'b', mf.mo_energy[1], False, False)
