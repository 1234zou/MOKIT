from py2fch import molinfo2fch
import numpy as np
from pyscf.data import elements, nist
from pyscf.gto.mole import ANG_OF, NPRIM_OF, NCTR_OF, PTR_EXP, PTR_COEFF, \
        gto_norm


def mol2fch(mol, fchname='test.fch'):
    nbf = mol.nao
    nif = nbf 
    na, nb = mol.nelec
    #ncontr = mol.nbas
    #print(nbf,ncontr)
    #nprim = 0
    charge = mol.charge
    mult = mol.spin+1
    natom = mol.natm
    LenNCZ = 0
    ielem = [elements.charge(mol.atom_symbol(a)) for a in range(natom)]
    print('_bas', mol._bas)
    ncontr = 0
    nprimitive = 0
    shell_type = []
    prim_per_shell = []
    shell2atom_map = []
    exps = []
    ccoeffs = []
    for a, sh in enumerate(mol._bas):
        contr_in_sh = sh[NCTR_OF]
        for c in range(contr_in_sh):
            sh_type = sh[ANG_OF]
            if not mol.cart and sh_type > 1:
                sh_type *= -1
            shell_type.append(sh_type)
            nprim = sh[NPRIM_OF]
            prim_per_shell.append(nprim)
            shell2atom_map.append(mol.bas_atom(a)+1)
            ncontr += 1
            nprimitive += nprim
        ptr_exp = sh[PTR_EXP]
        ptr_c = sh[PTR_COEFF]
        for c in range(contr_in_sh):
            norm = gto_norm(sh[ANG_OF], mol._env[ptr_exp:ptr_exp+nprim])
            sh_exps = mol._env[ptr_exp:ptr_exp+nprim]
            sh_coeffs = mol._env[ptr_c+c*nprim : ptr_c+nprim+c*nprim]
            exps.append( np.array(sh_exps)  )
            ccoeffs.append(np.array(sh_coeffs) / np.array(norm) )
    print(shell_type, prim_per_shell, shell2atom_map)
    virial = 0.0
    tot_e = 0.0
    coor = np.array([mol._atom[a][1] for a in range(natom)]).T
    coor = coor * nist.BOHR
    print(coor)

    #print(mol._basis)
    prim_exp = np.concatenate(exps)
    contr_coeff = np.concatenate(ccoeffs)
    contr_coeff_sp = np.zeros(1)

    molinfo2fch(fchname, nbf, nif, na, nb, ncontr, nprimitive, charge, mult, natom, LenNCZ, 
             ielem, shell_type, prim_per_shell, shell2atom_map, 
            # KFirst, KLast, Lmax, LPSkip, NLP, RNFroz, CLP, ZLP,&
             virial, tot_e, coor, prim_exp, contr_coeff, contr_coeff_sp)

