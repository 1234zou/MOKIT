from pyscf import gto

mol = gto.M(atom='O 0.0 0.0 0.1; H 0.0 0.0 1.0', basis='cc-pvtz', charge=-1).build()
mf = mol.RHF()
mf.conv_tol *= 0.1 # tighten converge tolerence
mf.run()

from mokit.lib.py2fch_direct import fchk
fchk(mf, '07-OH.fch')

# To run this example, at least py2fch need to be compiled by `make py2fch`
