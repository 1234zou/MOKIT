from pyscf import gto, scf
from mokit.lib.fch2py import fch2py
from mokit.lib.ortho import check_cghf_orthonormal
from mokit.lib.py2fch import py2fch_cghf
from shutil import copyfile

mol = gto.M()
# 1 atom(s)
mol.atom = '''
Sg1             0.00000000        0.00000000        0.00000000
'''

mol.basis = {
'Sg1': gto.basis.parse('''
Sg    P 
    0.670770000E-01    0.100000000E+01
Sg    P 
    0.175343000E-01    0.100000000E+01
Sg    D 
    0.254217990E+02    0.100000000E+01
Sg    D 
    0.259876260E+01    0.100000000E+01
Sg    D 
    0.205607850E+01    0.100000000E+01
Sg    D 
    0.777535300E+00    0.100000000E+01
Sg    D 
    0.290524700E+00    0.100000000E+01
Sg    D 
    0.975429000E-01    0.100000000E+01
Sg    F 
    0.818033850E+01    0.100000000E+01
Sg    F 
    0.370892090E+01    0.100000000E+01
Sg    F 
    0.166504880E+01    0.100000000E+01
Sg    F 
    0.770885700E+00    0.100000000E+01
Sg    F 
    0.335943800E+00    0.100000000E+01
''')}

mol.ecp = {
'Sg1': gto.basis.parse_ecp('''
Sg nelec  78
Sg ul
2      2.6431999200       -2.8641700700
2      6.0714001700      -24.8322163000
2     15.7412005000      -65.2087936000
2     42.4580002000     -230.8303070000
2    142.6430970000     -581.9624020000
1    474.1123960000      -72.6012726000
Sg S
2      0.8629999800       12.0596018000
2      1.0542000500      -40.3617897000
2      1.4940999700       83.7126770000
2      2.3296001000     -160.9820400000
2      4.0693998300      368.8686220000
2      7.5616998700     -378.9625550000
1     12.4813995000      154.1953740000
0     92.5594025000       13.1549702000
Sg P
2      1.4038000100      -26.9109764000       44.7906820000
2      1.6789000000      104.8424300000     -146.1085120000
2      2.2133998900     -225.3716580000      218.2392400000
2      3.2416999300      430.5595700000     -165.1309100000
2      5.0552001000     -450.6157230000       35.4875920000
2      8.1288003900      441.5639340000       48.3004620000
1     23.1009998000       36.7362557000      -11.1249820000
0     21.4099007000       15.1648626000        0.9273400000
Sg D
2      1.4240000200       46.3465157000        4.7736880000
2      1.7034000200     -151.6163330000      -17.0418770000
2      2.3399000200      313.8243100000       40.6165360000
2      3.5046999500     -357.8924870000      -62.3186650000
2      5.5816998500      368.8264770000       53.2819980000
2      8.8429002800     -140.9747620000      -27.6744170000
1     14.3563004000       74.6530151000        0.9829150000
0     46.2206993000        7.2974090600       -0.0982670000
Sg F
2      0.2009000000        0.0007810000       -0.0009326670
2      2.5058000100      -24.0355587000        1.0319386700
2      3.1310000400       67.8944473000       -5.5016973300
2      4.5275001500     -121.3988270000       14.1711120000
2      7.3171000500      232.5162050000      -24.2446193000
2     12.3481998000     -253.1587070000       26.3350140000
1     19.3710995000      110.2540970000       -0.6296486670
0    119.0846020000        8.2313900000        0.0608646670
''')}

# Remember to check the charge and spin
mol.charge = 0
mol.spin = 4
mol.verbose = 4
mol.build()

mf = scf.GHF(mol)
mf.with_soc = True
dm = mf.get_init_guess(key='1e') + 0j
mf.max_cycle = 1
mf.kernel(dm0=dm)

# read MOs from .fch(k) file
nbf = mf.mo_coeff.shape[0]
nif = mf.mo_coeff.shape[1]
mf.mo_coeff.real = fch2py('08-test1198.fch', nbf, nif, 'r')
mf.mo_coeff.imag = fch2py('08-test1198.fch', nbf, nif, 'i')
# read done

# check if input MOs are orthonormal
S = mol.intor_symmetric('int1e_ovlp')
check_cghf_orthonormal(nbf, nif, mf.mo_coeff, S)

dm = mf.make_rdm1()
mf.max_cycle = 10
mf.kernel(dm)  # will converge in 2 cycles
mf.stability() # internal instability detected

dm = mf.make_rdm1() + 0.01j # small perturbation to the density matrix
mf.max_cycle = 100
mf.kernel(dm0=dm)
mf.stability() # should be stable now

# save MOs of the stable complex GHF wave function into a new .fch file
copyfile('08-test1198.fch', '09-test1198.fch')
py2fch_cghf('09-test1198.fch', nbf, nif, mf.mo_coeff, mf.mo_energy, False)

