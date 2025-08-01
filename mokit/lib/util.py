import numpy as np
from mokit.lib.rwwfn import get_u, calc_expect_value

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
    npair = np.int64((idx[1]-idx[0]-idx[2])/2)
    if (npair < 0):
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

def get_Fii(nbf, nif, mo_a, mo, ev_a):
    u = get_u(nbf, nif, mo_a, mo)
    new_ev = calc_expect_value(nif, u, ev_a)
    return new_ev

def get_Fii_native(mf, mo, mo_occ):
    mf_r = mf.to_rhf()
    dm1 = mf_r.make_rdm1(mo, mo_occ)
    fock_ao = mf_r.get_fock(dm=dm1)
    Fii = np.einsum('mi,mn,ni->i', mo, fock_ao, mo, optimize=True)
    return Fii

def average_nmr_shielding(nmr_out, atom_list, program='gaussian', paramagnetic=False,
                          start_from_one=False):
    '''
    find (p)NMR isotropic shieldings of target atoms in the given output file and
    calculate the average value
    '''
    from mokit.lib.rwwfn import average_nmr_shield_in_gau_log, \
     average_nmr_shield_in_orca_out, average_pnmr_shield_in_orca_pnmr_out

    natom = len(atom_list)
    if start_from_one is False:
        new_list = [i + 1 for i in atom_list]
    else:
        new_list = atom_list

    if program == 'gaussian':
        if paramagnetic is True:
            raise ValueError('pNMR calculation is not supported in Gaussian.')
        else:
            ave_val = average_nmr_shield_in_gau_log(nmr_out, natom, new_list)
    elif program == 'orca':
        if paramagnetic is True:
            ave_val = average_pnmr_shield_in_orca_pnmr_out(nmr_out, natom, new_list)
        else:
            ave_val = average_nmr_shield_in_orca_out(nmr_out, natom, new_list)
    else:
        raise ValueError('program not supported.')
    return ave_val

