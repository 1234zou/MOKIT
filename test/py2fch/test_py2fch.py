from pyscf import gto

mol = gto.M(atom='O 0.0 0.0 0.1; H 0.0 0.0 1.0', 
        basis='cc-pvtz', 
        charge=-1,
        cart=True
        ).build()
mf = mol.RHF().run()

from mokit.lib.py2fch_direct import fchk
fchk(mf, 'test.fch')

# another equivalent way, for debugging
#from mokit.lib.py2fch_direct import mol2fch
#from mokit.lib.py2fch import py2fch
#mol2fch(mol, 'test.fch')
#py2fch('test.fch', mol.nao, mol.nao, mf.mo_coeff, 'a', mf.mo_energy, 0)

