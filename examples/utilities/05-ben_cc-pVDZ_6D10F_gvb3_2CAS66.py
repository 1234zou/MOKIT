from pyscf import gto, scf, mcscf
from mokit.lib.fch2py import fch2py
from mokit.lib.py2fch import py2fch
from shutil import copyfile

mol = gto.M()
#  12 atom(s)
mol.atom = '''
C           0.00000000        1.38776300        0.00000000
C           1.20183800        0.69388200        0.00000000
C           1.20183800       -0.69388200        0.00000000
C           0.00000000       -1.38776300        0.00000000
C          -1.20183800       -0.69388200        0.00000000
C          -1.20183800        0.69388200        0.00000000
H           0.00000000        2.46867400        0.00000000
H           2.13793400        1.23433700        0.00000000
H           2.13793400       -1.23433700        0.00000000
H           0.00000000       -2.46867400        0.00000000
H          -2.13793400       -1.23433700        0.00000000
H          -2.13793400        1.23433700        0.00000000
'''

mol.basis = {
'C': gto.basis.parse('''
C     S 
    0.666500000E+04    0.691583963E-03
    0.100000000E+04    0.532579615E-02
    0.228000000E+03    0.270607210E-01
    0.647100000E+02    0.101656846E+00
    0.210600000E+02    0.274574824E+00
    0.749500000E+01    0.448294319E+00
    0.279700000E+01    0.284902611E+00
    0.521500000E+00    0.151948592E-01
C     S 
    0.666500000E+04   -0.293269653E-03
    0.100000000E+04   -0.231803547E-02
    0.228000000E+03   -0.114997860E-01
    0.647100000E+02   -0.468267270E-01
    0.210600000E+02   -0.128466169E+00
    0.749500000E+01   -0.301266272E+00
    0.279700000E+01   -0.255630702E+00
    0.521500000E+00    0.109379336E+01
C     S 
    0.159600000E+00    0.100000000E+01
C     P 
    0.943900000E+01    0.569792516E-01
    0.200200000E+01    0.313207212E+00
    0.545600000E+00    0.760376742E+00
C     P 
    0.151700000E+00    0.100000000E+01
C     D 
    0.550000000E+00    0.100000000E+01
'''),
'C': gto.basis.parse('''
C     S 
    0.666500000E+04    0.691583963E-03
    0.100000000E+04    0.532579615E-02
    0.228000000E+03    0.270607210E-01
    0.647100000E+02    0.101656846E+00
    0.210600000E+02    0.274574824E+00
    0.749500000E+01    0.448294319E+00
    0.279700000E+01    0.284902611E+00
    0.521500000E+00    0.151948592E-01
C     S 
    0.666500000E+04   -0.293269653E-03
    0.100000000E+04   -0.231803547E-02
    0.228000000E+03   -0.114997860E-01
    0.647100000E+02   -0.468267270E-01
    0.210600000E+02   -0.128466169E+00
    0.749500000E+01   -0.301266272E+00
    0.279700000E+01   -0.255630702E+00
    0.521500000E+00    0.109379336E+01
C     S 
    0.159600000E+00    0.100000000E+01
C     P 
    0.943900000E+01    0.569792516E-01
    0.200200000E+01    0.313207212E+00
    0.545600000E+00    0.760376742E+00
C     P 
    0.151700000E+00    0.100000000E+01
C     D 
    0.550000000E+00    0.100000000E+01
'''),
'C': gto.basis.parse('''
C     S 
    0.666500000E+04    0.691583963E-03
    0.100000000E+04    0.532579615E-02
    0.228000000E+03    0.270607210E-01
    0.647100000E+02    0.101656846E+00
    0.210600000E+02    0.274574824E+00
    0.749500000E+01    0.448294319E+00
    0.279700000E+01    0.284902611E+00
    0.521500000E+00    0.151948592E-01
C     S 
    0.666500000E+04   -0.293269653E-03
    0.100000000E+04   -0.231803547E-02
    0.228000000E+03   -0.114997860E-01
    0.647100000E+02   -0.468267270E-01
    0.210600000E+02   -0.128466169E+00
    0.749500000E+01   -0.301266272E+00
    0.279700000E+01   -0.255630702E+00
    0.521500000E+00    0.109379336E+01
C     S 
    0.159600000E+00    0.100000000E+01
C     P 
    0.943900000E+01    0.569792516E-01
    0.200200000E+01    0.313207212E+00
    0.545600000E+00    0.760376742E+00
C     P 
    0.151700000E+00    0.100000000E+01
C     D 
    0.550000000E+00    0.100000000E+01
'''),
'C': gto.basis.parse('''
C     S 
    0.666500000E+04    0.691583963E-03
    0.100000000E+04    0.532579615E-02
    0.228000000E+03    0.270607210E-01
    0.647100000E+02    0.101656846E+00
    0.210600000E+02    0.274574824E+00
    0.749500000E+01    0.448294319E+00
    0.279700000E+01    0.284902611E+00
    0.521500000E+00    0.151948592E-01
C     S 
    0.666500000E+04   -0.293269653E-03
    0.100000000E+04   -0.231803547E-02
    0.228000000E+03   -0.114997860E-01
    0.647100000E+02   -0.468267270E-01
    0.210600000E+02   -0.128466169E+00
    0.749500000E+01   -0.301266272E+00
    0.279700000E+01   -0.255630702E+00
    0.521500000E+00    0.109379336E+01
C     S 
    0.159600000E+00    0.100000000E+01
C     P 
    0.943900000E+01    0.569792516E-01
    0.200200000E+01    0.313207212E+00
    0.545600000E+00    0.760376742E+00
C     P 
    0.151700000E+00    0.100000000E+01
C     D 
    0.550000000E+00    0.100000000E+01
'''),
'C': gto.basis.parse('''
C     S 
    0.666500000E+04    0.691583963E-03
    0.100000000E+04    0.532579615E-02
    0.228000000E+03    0.270607210E-01
    0.647100000E+02    0.101656846E+00
    0.210600000E+02    0.274574824E+00
    0.749500000E+01    0.448294319E+00
    0.279700000E+01    0.284902611E+00
    0.521500000E+00    0.151948592E-01
C     S 
    0.666500000E+04   -0.293269653E-03
    0.100000000E+04   -0.231803547E-02
    0.228000000E+03   -0.114997860E-01
    0.647100000E+02   -0.468267270E-01
    0.210600000E+02   -0.128466169E+00
    0.749500000E+01   -0.301266272E+00
    0.279700000E+01   -0.255630702E+00
    0.521500000E+00    0.109379336E+01
C     S 
    0.159600000E+00    0.100000000E+01
C     P 
    0.943900000E+01    0.569792516E-01
    0.200200000E+01    0.313207212E+00
    0.545600000E+00    0.760376742E+00
C     P 
    0.151700000E+00    0.100000000E+01
C     D 
    0.550000000E+00    0.100000000E+01
'''),
'C': gto.basis.parse('''
C     S 
    0.666500000E+04    0.691583963E-03
    0.100000000E+04    0.532579615E-02
    0.228000000E+03    0.270607210E-01
    0.647100000E+02    0.101656846E+00
    0.210600000E+02    0.274574824E+00
    0.749500000E+01    0.448294319E+00
    0.279700000E+01    0.284902611E+00
    0.521500000E+00    0.151948592E-01
C     S 
    0.666500000E+04   -0.293269653E-03
    0.100000000E+04   -0.231803547E-02
    0.228000000E+03   -0.114997860E-01
    0.647100000E+02   -0.468267270E-01
    0.210600000E+02   -0.128466169E+00
    0.749500000E+01   -0.301266272E+00
    0.279700000E+01   -0.255630702E+00
    0.521500000E+00    0.109379336E+01
C     S 
    0.159600000E+00    0.100000000E+01
C     P 
    0.943900000E+01    0.569792516E-01
    0.200200000E+01    0.313207212E+00
    0.545600000E+00    0.760376742E+00
C     P 
    0.151700000E+00    0.100000000E+01
C     D 
    0.550000000E+00    0.100000000E+01
'''),
'H': gto.basis.parse('''
H     S 
    0.130100000E+02    0.334987264E-01
    0.196200000E+01    0.234800801E+00
    0.444600000E+00    0.813682958E+00
H     S 
    0.122000000E+00    0.100000000E+01
H     P 
    0.727000000E+00    0.100000000E+01
'''),
'H': gto.basis.parse('''
H     S 
    0.130100000E+02    0.334987264E-01
    0.196200000E+01    0.234800801E+00
    0.444600000E+00    0.813682958E+00
H     S 
    0.122000000E+00    0.100000000E+01
H     P 
    0.727000000E+00    0.100000000E+01
'''),
'H': gto.basis.parse('''
H     S 
    0.130100000E+02    0.334987264E-01
    0.196200000E+01    0.234800801E+00
    0.444600000E+00    0.813682958E+00
H     S 
    0.122000000E+00    0.100000000E+01
H     P 
    0.727000000E+00    0.100000000E+01
'''),
'H': gto.basis.parse('''
H     S 
    0.130100000E+02    0.334987264E-01
    0.196200000E+01    0.234800801E+00
    0.444600000E+00    0.813682958E+00
H     S 
    0.122000000E+00    0.100000000E+01
H     P 
    0.727000000E+00    0.100000000E+01
'''),
'H': gto.basis.parse('''
H     S 
    0.130100000E+02    0.334987264E-01
    0.196200000E+01    0.234800801E+00
    0.444600000E+00    0.813682958E+00
H     S 
    0.122000000E+00    0.100000000E+01
H     P 
    0.727000000E+00    0.100000000E+01
'''),
'H': gto.basis.parse('''
H     S 
    0.130100000E+02    0.334987264E-01
    0.196200000E+01    0.234800801E+00
    0.444600000E+00    0.813682958E+00
H     S 
    0.122000000E+00    0.100000000E+01
H     P 
    0.727000000E+00    0.100000000E+01
''')}

# Remember to check the charge and spin
mol.charge = 0
mol.spin = 0
mol.verbose = 4
mol.cart = True
mol.build()

mf = scf.RHF(mol)
mf.max_cycle = 1
mf.kernel()

# read MOs from .fch(k) file
nbf = mf.mo_coeff.shape[0]
nif = mf.mo_coeff.shape[1]
mf.mo_coeff = fch2py('03-ben_cc-pVDZ_6D10F_gvb3.fch', nbf, nif, 'a')
# read done

mc = mcscf.CASSCF(mf, 6, 6)
mc.verbose = 5
mc.natorb = True
mc.kernel()

# export MOs into a new .fch(k) file
copyfile('03-ben_cc-pVDZ_6D10F_gvb3.fch','03-ben_cc-pVDZ_6D10F_CASSCF66_NO.fch')
py2fch('03-ben_cc-pVDZ_6D10F_CASSCF66_NO.fch',nbf,nif,mc.mo_coeff,'a',mc.mo_occ,True,True)
