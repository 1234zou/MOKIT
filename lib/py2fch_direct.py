from py2fch import molinfo2fch, molecp2fch
import numpy as np
from pyscf.data import elements, nist
from pyscf.gto.mole import ANG_OF, NPRIM_OF, NCTR_OF, PTR_EXP, PTR_COEFF, \
        gto_norm


def mol2fch(mol, fchname='test.fch', uhf=False, mo=None, trim_zeros=True):
    nbf = mol.nao
    if mo is not None:
        if uhf:
            nif = mo[0].shape[1]
        else:
            nif = mo.shape[1]
    else:
        nif = nbf 
    na, nb = mol.nelec
    #ncontr = mol.nbas
    #print(nbf,ncontr)
    charge = mol.charge
    mult = mol.spin+1
    natom = mol.natm
    ielem = [elements.charge(mol.atom_symbol(a)) for a in range(natom)]
    ncontr = 0
    nprimitive = 0
    shell_type = []
    prim_per_shell = []
    shell2atom_map = []
    exps = []
    ccoeffs = []
    #print(mol._bas)
    for a, sh in enumerate(mol._bas):
        contr_in_sh = sh[NCTR_OF]
        ptr_exp = sh[PTR_EXP]
        ptr_c = sh[PTR_COEFF]
        for c in range(contr_in_sh):
            sh_type = sh[ANG_OF]
            if not mol.cart and sh_type > 1:
                sh_type *= -1
            shell_type.append(sh_type)
            nprim = sh[NPRIM_OF]
            norm = np.array( gto_norm(sh[ANG_OF], mol._env[ptr_exp:ptr_exp+nprim]) )
            sh_exps = np.array( mol._env[ptr_exp:ptr_exp+nprim] )
            sh_coeffs = np.array( mol._env[ptr_c+c*nprim : ptr_c+nprim+c*nprim] )
            if trim_zeros:
                nprim = np.count_nonzero(sh_coeffs)
                non0_idx = np.flatnonzero(sh_coeffs)
                norm = norm[non0_idx]
                sh_exps = sh_exps[non0_idx]
                sh_coeffs = sh_coeffs[non0_idx]
            ncontr += 1
            nprimitive += nprim
            prim_per_shell.append(nprim)
            shell2atom_map.append(mol.bas_atom(a)+1)
            exps.append( sh_exps )
            ccoeffs.append( sh_coeffs / norm )
    #print(mol._basis)
    prim_exp = np.concatenate(exps)
    contr_coeff = np.concatenate(ccoeffs)
    contr_coeff_sp = np.zeros(1)
    virial = 0.0
    tot_e = 0.0
    coor = np.array([mol._atom[a][1] for a in range(natom)]).T
    coor = coor * nist.BOHR
    #print(coor)
    #print(shell_type, prim_per_shell, shell2atom_map)
    
    LenNCZ = 0
    KFirst = np.zeros((natom,10), dtype=int)
    KLast = np.zeros((natom,10), dtype=int)
    Lmax = np.zeros(natom, dtype=int)
    LPSkip = np.ones(natom, dtype=int)
    RNFroz = np.zeros(natom)

    if mol._ecpbas.size > 0:
        #print(mol._ecpbas)
        NLP = []
        CLP = []
        ZLP = []
    
        knum = np.zeros((natom,10), dtype=int)
        for a, item in enumerate(mol._ecpbas):
            iatom, lb, nexp, rorder, _, ptr_exp, ptr_coeff, _ = item
            Lmax[iatom] = max(lb+1, Lmax[iatom])
            LPSkip[iatom] = 0
            for c in range(nexp):
                NLP.append(rorder)
                e = mol._env[ptr_exp+c]
                c = mol._env[ptr_coeff+c]
                ZLP.append(e)
                CLP.append(c)
                knum[iatom][lb+1] += 1
        #print(NLP, ZLP, CLP)
        #print(knum)
        LenNCZ = len(NLP)
        #print(LenNCZ, Lmax)
        kl = 0
        for i in range(natom):
            for k,n in enumerate(knum[i]):
                if n == 0: continue
                kf = kl
                kl = kf + n
                KFirst[i,k] = kf + 1
                KLast[i,k] = kl
        #print(KFirst, KLast)
        RNFroz += np.array(ielem) 
        RNFroz -= mol.atom_charges()
        #print(RNFroz, RNFroz.dtype)
        NLP = np.array(NLP)
        CLP = np.array(CLP)
        ZLP = np.array(ZLP)
    else:
        NLP = np.zeros(0)
        CLP = np.zeros(0)
        ZLP = np.zeros(0)
    
    if LenNCZ > 0:
        molecp2fch(fchname, uhf, nbf, nif, na, nb, ncontr, nprimitive, charge, mult, natom, LenNCZ, 
             ielem, shell_type, prim_per_shell, shell2atom_map, 
             virial, tot_e, coor, prim_exp, contr_coeff, #contr_coeff_sp,
             KFirst, KLast, Lmax, LPSkip, NLP, RNFroz, CLP, ZLP
             )
    else:
        molinfo2fch(fchname, uhf, nbf, nif, na, nb, ncontr, nprimitive, charge, mult, natom, LenNCZ, 
             ielem, shell_type, prim_per_shell, shell2atom_map, 
             virial, tot_e, coor, prim_exp, contr_coeff, #contr_coeff_sp,
             #KFirst, KLast, Lmax, LPSkip, NLP, RNFroz, CLP, ZLP
             )

def fchk(mf, fchname, density=False, overwrite_mol=False):
    from pyscf import scf, mcscf
    import os
    from py2fch import py2fch
    is_uhf = isinstance(mf, scf.uhf.UHF)
    mo = mf.mo_coeff
    if (not os.path.isfile(fchname)) or overwrite_mol:
        mol2fch(mf.mol, fchname, is_uhf, mo)
    if isinstance(mf, scf.hf.SCF):
        if isinstance(mf, (scf.ghf.GHF, scf.dhf.DHF)):
            raise NotImplementedError('GHF/DHF not supported in py2fch')
        if not is_uhf: # ROHF is also RHF here
            py2fch(fchname, mo.shape[0], mo.shape[1], mo, 'a', mf.mo_energy, density)
        else:
            py2fch(fchname, mo[0].shape[0], mo[0].shape[1], mo[0], 'a', mf.mo_energy[0], density)
            py2fch(fchname, mo[1].shape[0], mo[1].shape[1], mo[1], 'b', mf.mo_energy[1], density)
    elif isinstance(mf, mcscf.casci.CASCI):
        py2fch(fchname, mo.shape[0], mo.shape[1], mo, 'a', mf.mo_occ, density)
        
