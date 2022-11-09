from pyscf import gto
from py2fch_direct import fchk

mol = gto.M(atom='O 0.0 0.0 0.1; H 0.0 0.0 1.0', 
        basis='cc-pvtz', 
        charge=-1,
        verbose=4
        ).build()

mf = mol.RHF().run()
fchk(mf, 'test.fch')

mf2 = mol.UHF().run()
fchk(mf2, 'test_uhf.fch')

from pyscf import mcscf, cc
mf3 = mcscf.CASSCF(mf, 5, 8)
mf3.natorb = True
mf3.kernel()
fchk(mf3, 'test_cas.fch')

mf4 = cc.CCSD(mf).set_frozen().run()
fchk(mf4, 'test_cc.fch') # will fail because not implemented
