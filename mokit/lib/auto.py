
'''
automr utilities and APIs for automr workflow
'''

import shutil
from mokit.lib.gaussian import get_bfirst_from_mol, load_mol_from_fch
from mokit.lib.rwwfn import (
    read_nbf_and_nif_from_fch,
    read_na_and_nb_from_fch,
    read_eigenvalues_from_fch,
    calc_diag_gross_pop,
    get_mo_center_from_pop
)
from mokit.lib.fch2py import fch2py
from mokit.lib.py2fch import py2fch
    
def find_antibonding_orb(mol, mo, i1=0, i2=0, i3=0, start_from_one=False,
                         ao_ovlp=None, popm='lowdin'):
    '''
    Construct antibonding orbitals for a set of bonding orbitals.
    mol: PySCF molecule object
    mo: molecular orbital coefficients
    k1, k2, k3: orbital indices
    start_from_one: whether the given orbital indices starts from 1 (Fortran
     convention). start_from_one=False means starting from 0 (Python convention).
    ao_ovlp: AO overlap integral matrix. If `None` is given, it will be calculated
     using mol.intor_symmetric('int1e_ovlp') below, otherwise it will be directly
     used.
    popm: population method, 'mulliken' or 'lowdin'
    Note:
    1) mo(:,i3:nif) will be updated, where mo(:,i3:i3+i2-i1) are generated anti-
     bonding orbitals, and mo(:,i3+i2-i1+1:nif) are remaining virtual orbitals.
    2) all MOs are still orthonormalized after calling this function.
    '''
    from mokit.lib.gaussian import get_ao_dip
    from mokit.lib.ortho import check_orthonormal
    from mokit.lib.wfn_analysis import find_antibonding_orbitals

    if start_from_one is True:
        k1 = i1
        k2 = i2
        k3 = i3
    else:
        k1 = i1+1
        k2 = i2+1
        k3 = i3+1
    npair = k2 - k1 + 1
    nbf = mo.shape[0]
    nif = mo.shape[1]

    if k2<k1 or k3<=k2 or nif<k3+npair-1:
        print('i1= %d, i2= %d, i3= %d, nif= %d' %(i1, i2, i3, nif))
        raise ValueError('Wrong orbital indices.')

    natom = mol.natm
    bfirst = get_bfirst_from_mol(mol)

    if ao_ovlp is None:
        S = mol.intor_symmetric('int1e_ovlp')
    else:
        S = ao_ovlp
    pop = calc_diag_gross_pop(natom, nbf, npair, bfirst, S, mo[:,k1-1:k2], popm)
    mo_center = get_mo_center_from_pop(natom, npair, pop)
    center, ao_dip = get_ao_dip(mol, fix_center=True)
    new_mo = find_antibonding_orbitals(k1, k2, k3, natom, nbf, nif, bfirst,
                                       mo_center, S, ao_dip, mo)
    check_orthonormal(nbf, nif, new_mo, S)
    return new_mo


def find_antibonding_orb_in_fch(fchname, i1, i2, i3, start_from_one=False, popm='lowdin'):
    '''
    Construct antibonding orbitals for a set of bonding orbitals. The original MOs
    are stored in fchname.
    i1, i2, i3: orbital indices
    start_from_one: whether the given orbital indices starts from 1 (Fortran
     convention). start_from_one=False means starting from 0 (Python convention).
    popm: population method, 'mulliken' or 'lowdin'
    Note:
    1) mo(:,i3:nif) will be updated, where mo(:,i3:i3+i2-i1) are generated anti-
     bonding orbitals, and mo(:,i3+i2-i1+1:nif) are remaining virtual orbitals.
    2) all MOs are still orthonormalized after calling this function.
    '''
    mol = load_mol_from_fch(fchname)
    nbf, nif = read_nbf_and_nif_from_fch(fchname)
    mo = fch2py(fchname, nbf, nif, 'a')
    new_mo = find_antibonding_orb(mol, mo, i1, i2, i3, start_from_one, None, popm)
    new_fch = fchname[0:fchname.rindex('.fch')]+'_a.fch'
    shutil.copyfile(fchname, new_fch)
    ev = read_eigenvalues_from_fch(fchname, nif, 'a')
    py2fch(new_fch, nbf, nif, new_mo, 'a', ev, False, False)

