
'''
automr utilities
'''

from mokit.lib.gaussian import loc, load_mol_from_fch
from mokit.lib.rwwfn import read_nbf_and_nif_from_fch, read_na_and_nb_from_fch, orb_resemble_ref1, get_occ_from_na_nb
from pyscf import gto
from mokit.lib.fch2py import fch2py
from mokit.lib.lo import gen_loc_ini_guess, pm, boys
from mokit.lib.gaussian import BOHR2ANG
import numpy as np

def proj_vir2virlmo(hf_fch, minbas_fch, minbas_loc_fch, mo, mol, localm='pm'):
    na, nb = read_na_and_nb_from_fch(hf_fch)
    na2, nb2 = read_na_and_nb_from_fch(minbas_fch)
    nbf2, nif2 = read_nbf_and_nif_from_fch(minbas_fch)
    mol2 = load_mol_from_fch(minbas_fch)
    nbf, nif = mo.shape
    nvir = nif - na
    npair = nif2 - na2
    vir_idx = range(na,nif)
    vir_idx2 = range(na2,nif2)
    loc(fchname=minbas_fch,idx=vir_idx2,method=localm)
    cross_s = gto.intor_cross('int1e_ovlp', mol, mol2)
    mo2 = fch2py(minbas_loc_fch, nbf2, nif2, 'a')
    new_vir = orb_resemble_ref1(nbf, nvir, mo[:,vir_idx], nbf2, npair, mo2[:,vir_idx2], cross_s)
    mo[:,vir_idx] = new_vir.copy()
    mo_occ = get_occ_from_na_nb(nif, na, nb)
    return mo, mo_occ, npair

def loc_ini_guess(mol, occ_mo, nval):
    natom = mol.natm
    nbf = occ_mo.shape[0]
    S = mol.intor_symmetric('int1e_ovlp')
    chosen = np.ones(natom, dtype=bool)
    bfirst = np.ones(natom+1, dtype=np.int32)
    bfirst[1:] = mol.aoslice_by_atom()[:,3] + 1
    nmo1, lmo_ini = gen_loc_ini_guess(natom,nbf,nval,chosen,bfirst,S,occ_mo)
    return nmo1, lmo_ini

def loc_driver(mol, lmo_ini, nval, method='pm', ao_dip=None, dis_tol=17.0, conv_tol=1e-5):
    '''
    TODO: can lmo_ini and nval be optional?
    '''
    natom = mol.natm
    nbf = mol.nao
    bfirst = np.ones(natom+1, dtype=np.int32)
    bfirst[1:] = mol.aoslice_by_atom()[:,3] + 1
    dis = BOHR2ANG*gto.inter_distance(mol)
    if method == 'pm':
        S = mol.intor_symmetric('int1e_ovlp')
        occ_lmo = pm(natom,nbf,nval,bfirst,dis,lmo_ini,S,'mulliken',dis_tol,conv_tol)
    elif method == 'boys':
        if ao_dip is None:
            # TODO: probably we can calculate ao_dip here
            raise ValueError('ao_dip should be provided for boys localization')
        occ_lmo = boys(natom,nbf,nval,bfirst,dis,lmo_ini,ao_dip,dis_tol,conv_tol)
    else:
        raise ValueError(f'Localization method {method} cannot be recognized.')
    return occ_lmo

def find_root_by_ss(mc, nroots, target_root, target_ss, iroot_init=-1):
    '''
    target_root: the target root index of certain spin, 1 for S1, 1 for T1
    iroot_init: set -1 for same spin excitation, 0 for different spin
    '''
    iroot = iroot_init
    for j in range(nroots):
        ss = mc.fcisolver.spin_square(mc.ci[j], mc.ncas, mc.nelecas)
        if abs(ss[0] - target_ss) < 1e-4:
            iroot = iroot + 1
        if iroot == target_root:
            return j
    raise ValueError(f'Cannot find root with target ss {target_ss} in {nroots} roots')

