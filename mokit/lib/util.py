#!/usr/bin/env python

'''
utility functions
1) needed by mokit.lib.gaussian, mokit.lib.auto
2) needed by automr workflows
3) misc
''' 

import shutil
import numpy as np
from mokit.lib.gaussian import get_bfirst_from_mol, load_mol_from_fch
from mokit.lib.rwwfn import (
    read_nbf_and_nif_from_fch,
    read_na_and_nb_from_fch,
    read_eigenvalues_from_fch
)
from mokit.lib.fch2py import fch2py
from mokit.lib.py2fch import py2fch

ON_criteria = 1e-5

def get_idx_from_noon(suno_out, noon, nbf, nif, nopen, thres=ON_criteria):
    ndb = np.count_nonzero(noon > 2.0-thres)
    nocc = np.count_nonzero(noon > thres)
    idx = [ndb+1, nocc+1, nopen]
    with open(suno_out, 'w') as f:
        f.write("nbf=    %5d\n" % nbf)
        f.write("nif=    %5d\n" % nif)
        f.write("ON_criteria= %11.7f\n" % ON_criteria)
        f.write("uno_thres= %11.7f\n" % thres)
        f.write("ndb=      %5d\n" % ndb)
        f.write("nocc=      %5d\n" % nocc)
        f.write("nopen=      %5d\n" % nopen)
        f.write("idx=      %5d%5d%5d\n" % (idx[0],idx[1],idx[2]))
    return idx


def get_npair_and_ovidx(idx, nskip_uno=0):
    npair = np.int32((idx[1]-idx[0]-idx[2])/2)
    if npair < 0:
        raise ValueError('npair<0. Impossible.')
    else:
        idx2 = idx[0] + npair - 1
        idx3 = idx2 + idx[2]
        idx1 = idx2 - npair
        idx4 = idx3 + npair
        i = nskip_uno # pair(s) of UNO to be skipped
        occ_idx = range(idx1,idx2-i)
        vir_idx = range(idx3+i,idx4)
    return npair, occ_idx, vir_idx


def get_nvir_from_fch(fchname):
    '''
    Get the number of virtual orbitals from a specified .fch(k) file. For R(O)HF,
    the result is nif-na. For UHF, the result is an integer array which cotains the
    number of alpha virtual orbitals and beta virtual orbitals.
    '''
    from mokit.lib.rwwfn import check_uhf_in_fch

    na, nb = read_na_and_nb_from_fch(fchname)
    nbf, nif = read_nbf_and_nif_from_fch(fchname)
    uhf = check_uhf_in_fch(fchname)
    if uhf == 0: # not UHF type
        nvir = nif - na
    else:        # UHF type
        nvir = np.array([nif-na, nif-nb])
    return nvir

# There exist multiple local minima of vir_idx orbital localization at target
# basis set. Since 20250113 this orbital localization is replaced by an
# orbital localization at the minimal basis set, followed by a virtual space
# projection using orb_resemble_ref1().
def proj_vir2virlmo(na, nb, mo, mol, mb_fch, localm='pm'):
    from pyscf import gto
    from mokit.lib.gaussian import loc
    from mokit.lib.rwwfn import orb_resemble_ref1, get_occ_from_na_nb

    mb_loc_fch = mb_fch[0:mb_fch.rindex('.fch')]+'_LMO.fch'
    na2, nb2 = read_na_and_nb_from_fch(mb_fch)
    nbf2, nif2 = read_nbf_and_nif_from_fch(mb_fch)
    mol2 = load_mol_from_fch(mb_fch)
    cross_S = gto.intor_cross('int1e_ovlp', mol, mol2)
    nbf, nif = mo.shape
    nvir = nif - na
    nvir2 = nif2 - na2
    vir_idx = range(na, nif)
    vir_idx2 = range(na2, nif2)
    loc(fchname=mb_fch, idx=vir_idx2, method=localm)
    mo2 = fch2py(mb_loc_fch, nbf2, nif2, 'a')
    new_vir = orb_resemble_ref1(nbf, nvir, mo[:,vir_idx], nbf2, nvir2,
                                mo2[:,vir_idx2], cross_S)
    mo[:,vir_idx] = new_vir.copy()
    mo_occ = get_occ_from_na_nb(nif, na, nb)
    return mo, mo_occ


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


def get_Fii(cmo, lmo, ao_ovlp, ev):
    # lmo = cmo*u, i.e. u = (cmo^T)*S*lmo
    u = np.einsum('ui,uv,vj->ij', cmo, ao_ovlp, lmo, optimize=True)
    # (U^T)eU
    new_ev = np.einsum('ji,j,ji->i', u, ev, u)
    return new_ev


def get_Fii_native(mf, mo, mo_occ):
    mf_r = mf.to_rhf()
    dm1 = mf_r.make_rdm1(mo, mo_occ)
    fock_ao = mf_r.get_fock(dm=dm1)
    Fii = np.einsum('mi,mn,ni->i', mo, fock_ao, mo, optimize=True)
    return Fii


def get_db_idx_from_molden(molden, on_thres=1e-6):
    '''
    Find indices of doubly occupied orbitals in a specified MOLDEN format file.
    Simple usage::
    >>> from mokit.lib.util import get_db_idx_from_molden
    >>> mo_idx = get_db_idx_from_molden(molden)
    '''
    from mokit.lib.rwwfn import read_nmo_from_molden, read_on_from_molden

    nmo, nmo_b = read_nmo_from_molden(molden)
    occ = read_on_from_molden(molden, nmo)
    return np.where(occ > 2.0-on_thres)[0]


def get_nbond_from_mol_and_lmo(mol, lmo, ao_ovlp=None, popm='lowdin'):
    '''
    Get the number of chemical bonds from the mol object and a set of occupied LMOs.
    Any occupied LMO which is located on more than 1 atom will be counted as one
    chemical bond. Note: a diradical MO is also viewed as a chemical bond.
    '''
    from mokit.lib.rwwfn import calc_diag_gross_pop, get_mo_center_from_pop

    natom = mol.natm
    nbf, nmo = lmo.shape
    bfirst = get_bfirst_from_mol(mol)
    if ao_ovlp is None:
        S = mol.intor_symmetric('int1e_ovlp')
    else:
        S = ao_ovlp
    pop = calc_diag_gross_pop(natom, nbf, nmo, bfirst, S, lmo, popm)
    mo_center = get_mo_center_from_pop(natom, nmo, pop)
    condition = mo_center[0,:] > 1
    indices_gt1 = np.where(condition)[0]
    indices_other = np.where(~condition)[0]
    new_order = np.concatenate([indices_other, indices_gt1])
    lmo_sorted = lmo[:,new_order]
    return np.sum(condition), lmo_sorted


def find_antibonding_orb(mol, mo, i1=0, i2=0, i3=0, start_from_one=False,
                         ao_ovlp=None, popm='lowdin'):
    raise OSError('This function is moved to auto.py. Please use `from mokit.lib.auto import`')


def find_antibonding_orb_in_fch(fchname, i1, i2, i3, start_from_one=False, popm='lowdin'):
    raise OSError('This function is moved to auto.py. Please use `from mokit.lib.auto import`')


def get_t1_t2_from_stc_ri_ccsd(mf, ncore=0, ccd=False, target_error=1e-4, niters=14,
                               diis_space=6, diis_shift=2):
    '''
    Get t1 and t2 amplitudes (based on canonical MOs) from an STC-RI-CCSD calculation.
    The t1 and t2 amplitudes can be used as the initial guess of a (RI-)CCSD calculation.
    Note: the input mf can either be a density fitting (DF) or non-DF RHF object.
    This function can only be applied to the RHF-(RI-)CCSD case currently.
    '''
    nocc = np.count_nonzero(mf.mo_occ > 0)
    nvir = mf.mo_coeff.shape[1] - nocc
    nocc = nocc - ncore

    try:
        import stc_cc
    except ImportError:
        from_stc = False
        print('\nSTC-CC library not detected/installed. Use MP2 amplitudes...')
    else:
        # frozen core CCSD/cc-pVTZ of benzene has nocc=15 and nvir=243
        if nocc > 14 and nvir > 242:
            from_stc = True
            print('\nSTC-CC library detected. Try to use STC-RI-CCSD t1 and t2...')
        else:
            from_stc = False
            print('\nSTC-CC library detected. But both nocc and nvir are small.')
            print('Still use MP2 amplitudes...')

    if from_stc:
        if not hasattr(mf, 'with_df') or mf.with_df is None:
            mf2 = mf.density_fit() # use even-tempered gaussian functions as fitting basis
            mf2.max_cycle = 0
            mf2.kernel()
            t1, t2 = get_t1_t2_from_stc_ri_ccsd2(mf2, ncore, ccd, target_error,
                                                 niters, diis_space, diis_shift)
        else:
            t1, t2 = get_t1_t2_from_stc_ri_ccsd2(mf, ncore, ccd, target_error,
                                                 niters, diis_space, diis_shift)
        return t1, t2, True
    else:
        from pyscf import mp
        mc = mp.MP2(mf)
        mc.verbose = 4
        mc.frozen = ncore
        mc.kernel()
        # TODO: implement t1 using 2e integrals and MP2 t2
        return np.zeros((nocc,nvir)), mc.t2, False


def get_t1_t2_from_stc_ri_ccsd2(mf, ncore=0, ccd=False, target_error=1e-4, niters=14,
                                diis_space=6, diis_shift=2):
    '''
    Get t1 and t2 amplitudes (based on canonical MOs) from an STC-RI-CCSD calculation.
    The t1 and t2 amplitudes can be used as the initial guess of a (RI-)CCSD calculation.
    Note: the input mf must be a DF-RHF object.
    '''
    # See the STC-CC repo (https://github.com/SUSYUSTC/stc_qc_paper).
    # This module is originally modified from stc_qc_paper/run/stc_cc_general.py
    # Thanks to Jiace Sun's help for detailed explanations on the usage of STC-CC.
    # Any user who uses this module is supposed to cite arXiv:2602.17158.
    # We will update the arXiv number by formal DOI and link once it is accepted.

    # TODO: check memory requirement
    # os.environ['PYSCF_MAX_MEMORY'] = '100000'
    from stc_cc import cc, utils, utils_pyscf, sample, stc
    from run import stc_orb
    import torch

    # some minor setups, does not affect result, shoult be good for all systems and machines
    stc.quota_sampling = True
    cc.use_alias_for_X2vvoo = True
    cc.jit_estimation = True
    sample.pack_sign = True
    sample.align_alias = False

    # use S^{-1/2} for orthonormal auxiliary basis, which is more local, and is
    # a better initial guess for localization later.
    utils_pyscf.make_df_eig()

    singles = (not ccd)

    # TODO: check basis set linear dependent case
    nao = mf.mol.nao
    nocc = mf.mol.nelectron // 2
    nvir = nao - nocc
    o = slice(ncore, nocc)
    v = slice(nocc, None)
    nocc = nocc - ncore

    C = mf.mo_coeff
    Cocc = C[:, o]
    Cvir = C[:, v]
    C = np.concatenate((Cocc, Cvir), axis=1)

    # TODO: change to orbital localization of MOKIT, so that we can remove the
    # dependency of stc_orb
    PM_guess = True
    Cocc_guess, Cvir_guess = stc_orb.get_guess(mf.mol, Cocc, Cvir, PM=PM_guess)
    thres = 5e-3  # some tight value for safety
    Cocc_local, Cvir_local, Uaux = stc_orb.get_Rov(mf, Cocc_guess, Cvir_guess, thres)
    key_args = dict(canonical_denominator=True, minimal_stc=True)
    Clocal = np.concatenate((Cocc_local, Cvir_local), axis=1)

    block_size = 20
    padding = (False, True, True)

    # weights_out_power: gamma in the SI, Sec 4.1, Eq. 35
    contractor_args = dict(linear=False, singles=singles, **key_args, weights_out_power=0.5,
                           block_size=block_size, padding=padding, diis_space=diis_space,
                           diis_shift=diis_shift, jit=True, dynamic=True)
    print('options')
    for key, value in contractor_args.items():
        print(key, value)
    print()

    mycc = cc.StochasticCC(mf, Clocal, target_error, Uaux=Uaux, verbose=1, **contractor_args)
    mycc.initialize()
    mycc.set_freezing_mask(max_Eerr=None)

    # some trick to save memory, to allow deleting big input tensor in function calls
    # without leaving a reference outside the function
    wrap_ref = True
    X_tuple = mycc.get_X_tuple_guess(wrap_ref=wrap_ref)
    mycc.initialize_diis()

    vec_prev = mycc.flatten_X(X_tuple)
    Eij_init = mycc.get_Eij(X_tuple)
    E_init = Eij_init.sum()
    print('init energy', E_init.item())
    Es = [E_init]
    Eij_prev = Eij_init
    for it in range(niters):
        X_tuple, vec_raw = mycc.stochastic_CCSD_update(X_tuple, return_raw=False)
        Eij = mycc.get_Eij(X_tuple)
        Eij_diff = Eij - Eij_prev
        E = Eij.sum()
        Es.append(E)
        print(f'it{it} energy', E.item(), 'Eij update', torch.linalg.norm(Eij_diff).item())
        vec = mycc.flatten_X(X_tuple)
        print('vec diff', (torch.linalg.norm(vec - vec_prev) / torch.linalg.norm(vec_prev)).item())
        print()
        Eij_prev = Eij
        vec_prev = vec
        del vec_raw
    del vec, vec_prev

    print('factor', mycc.nsamples_history[-1] / mycc.nsamples_history[0])
    # average the last a few iterations
    print('averaged CCSD energy', torch.stack(Es[-4:]).mean().item())

    # extract t1 and t2 amplitudes of canonical MOs
    X_tuple_unwrapped = utils.get_reference_data(X_tuple)
    X1A, X1B, X2 = X_tuple_unwrapped
    X1A_canonical_full = mycc.X1_local_to_canonical(X1A)
    X1B_canonical_full = mycc.X1_local_to_canonical(X1B)
    X2_canonical_full = mycc.X2_local_to_canonical(X2)
    X1A_canonical = X1A_canonical_full[:nocc, :nvir].numpy() 
    X1B_canonical = X1B_canonical_full[:nocc, :nvir].numpy()
    X1_canonical = 0.5 * (X1A_canonical + X1B_canonical)
    X2_canonical = X2_canonical_full[:nocc, :nvir, :nocc, :nvir].numpy().transpose(0, 2, 1, 3)
    return X1_canonical, X2_canonical


def export_mo_e2txt(fchname):
    '''
    export the data of Alpha Orbital Energies in a .fch file into a plain text file
    '''
    from mokit.lib.rwwfn import export_rarray2txt

    txtname = fchname[0:fchname.rindex('.fch')]+'.txt'
    nbf, nif = read_nbf_and_nif_from_fch(fchname)
    ev = read_eigenvalues_from_fch(fchname, nif, 'a')
    export_rarray2txt(txtname, 'MO Eigenvalues', nif, ev)


def mo_svd_in_fch(fchname1, fchname2, idx1=None, idx2=None):
    '''
    Perform SVD on two sets of MOs in two .fch(k) files.
    idx1/idx2: the 1st/last index of the MO, starts from 0
    '''
    from mokit.lib.wfn_analysis import mo_svd_in2fch

    if idx1 is None:
        idx1 = 0
    if idx2 is None:
        nbf, idx2 = read_nbf_and_nif_from_fch(fchname1)
    print('idx1= %d, idx2= %d' %(idx1, idx2) )
    mo_svd_in2fch(fchname1, fchname2, idx1+1, idx2)


def average_nmr_shielding(nmr_out, atom_list, program='gaussian', paramagnetic=False,
                          start_from_one=False):
    '''
    find (paramagnetic) NMR isotropic shieldings of target atoms in the given
    output file and calculate the average value
    '''
    from mokit.lib.rwwfn import average_nmr_shield_in_gau_log, \
     average_nmr_shield_in_orca_out, average_pnmr_shield_in_orca_pnmr_out

    natom = len(atom_list)
    if start_from_one is False:
        new_list = [i + 1 for i in atom_list]
    else:
        new_list = atom_list

    if program == 'gaussian':
        if paramagnetic:
            raise ValueError('pNMR calculation is not supported in Gaussian.')
        else:
            ave_val = average_nmr_shield_in_gau_log(nmr_out, natom, new_list)
    elif program == 'orca':
        if paramagnetic:
            ave_val = average_pnmr_shield_in_orca_pnmr_out(nmr_out, natom, new_list)
        else:
            ave_val = average_nmr_shield_in_orca_out(nmr_out, natom, new_list)
    else:
        raise ValueError('program not supported.')
    return ave_val


def mo_g_int(fnames, x, na=None, nb=None, trace_PS=False):
    '''
    Generate occupied MOs of a new geometry using Grassmann interpolation. Currently
     only available for R(O)HF and UHF.
    fnames  : a series of .fch(k) files and a .gjf file
    x       : a series of the changed variable (bond distance, angle, dihedral, even
              composite coordinate)
    na      : the number of alpha occupied orbitals
    nb      : the number of beta occupied orbitals
    trace_PS: whether to calculate trace(PS) to check the number of electrons
    Example:
      mo_g_int(['h2o_105.fch', 'h2o_115.fch', 'h2o_120.fch', 'h2o_109_5.gjf'],
               [105.0, 115.0, 120.0, 109.5])
    '''
    from mokit.lib.qchem import read_hf_type_from_fch, construct_vir
    from mokit.lib.mirror_wfn import mo_grassmann_intrplt
    from mokit.lib.rwgeom import replace_coor_in_fch_by_gjf

    nfile = len(fnames)
    if nfile < 2:
        raise ValueError('At least two files must be provided.')
    if len(x) != nfile:
        raise ValueError('Size of arrays fnames and x are not equal.')

    # copy a .fch file and replace coordinates therein by coordinates from .gjf
    gjfname = fnames[nfile-1]
    new_fch = gjfname[0:gjfname.rindex('.gjf')]+'.fch'
    shutil.copyfile(fnames[0], new_fch)
    replace_coor_in_fch_by_gjf(gjfname, new_fch)

    nbf, nif = read_nbf_and_nif_from_fch(fnames[0])
    na0, nb0 = read_na_and_nb_from_fch(fnames[0])
    if na is None:
        na = na0
    if nb is None:
        nb = nb0

    S = np.zeros([nbf,nbf,nfile])
    mo = np.zeros([nbf,na,nfile-1])

    for i in range(nfile-1):
        mol = load_mol_from_fch(fnames[i])
        S[:,:,i] = mol.intor_symmetric('int1e_ovlp')
        coeff = fch2py(fnames[i], nbf, nif, 'a')
        mo[:,:na,i] = coeff[:,:na]

    mol = load_mol_from_fch(new_fch)
    S[:,:,nfile-1] = mol.intor_symmetric('int1e_ovlp')

    # generate alpha occupied MOs of the new geometry
    new_mo = mo_grassmann_intrplt(nbf, na, nfile, x, S, mo)
    coeff0 = np.zeros([nbf,nif])
    coeff0[:,:na] = new_mo[:,:na]

    # construct alpha virtual MOs of the new geometry using PAO
    coeff = construct_vir(nbf, nif, na+1, coeff0, S[:,:,nfile-1])

    # export alpha MOs to .fch file
    mo_e = np.zeros(nif)
    py2fch(new_fch, nbf, nif, coeff, 'a', mo_e, False, False)

    # check the number of alpha electrons
    if trace_PS is True:
        dm = np.dot(coeff[:,:na], coeff[:,:na].transpose())
        ne_a = np.trace(np.dot(dm,S[:,:,nfile-1]))
        print('No. alpha electrons: %.4f' %ne_a)

    # Note: DO NOT merge beta occupied MOs into the array mo for UHF, otherwise MOs
    #       in the array mo are non-orthogonal. In fact, alpha/beta should be dealt
    #       with separately.
    ihf = read_hf_type_from_fch(fnames[0])
    if ihf == 2:   # real UHF
        mo = np.zeros([nbf,nb,nfile-1])
        for i in range(nfile-1):
            coeff = fch2py(fnames[i], nbf, nif, 'b')
            mo[:,:nb,i] = coeff[:,:nb]
        # generate beta occupied MOs of the new geometry
        new_mo = mo_grassmann_intrplt(nbf, nb, nfile, x, S, mo)
        coeff0[:,:nb] = new_mo[:,:nb]
        # construct beta virtual MOs of the new geometry using PAO
        coeff = construct_vir(nbf, nif, nb+1, coeff0, S[:,:,nfile-1])
        # export beta MOs to .fch file
        py2fch(new_fch, nbf, nif, coeff, 'b', mo_e, False, False)
        # check the number of beta electrons
        if trace_PS is True:
            dm = np.dot(coeff[:,:nb], coeff[:,:nb].transpose())
            ne_b = np.trace(np.dot(dm,S[:,:,nfile-1]))
            print('No. beta electrons: %.4f' %ne_b)

