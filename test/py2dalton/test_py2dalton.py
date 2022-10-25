from pyscf import gto, scf
from py2dalton import py2dalton

mol = gto.M(atom='''
O  -0.49390246   0.93902438   0.0
H   0.46609754   0.93902438   0.0
H  -0.81435705   1.84396021   0.0
''',
basis='cc-pVDZ')

mf = scf.RHF(mol).run()
py2dalton(mf, 'h2o.dal')
