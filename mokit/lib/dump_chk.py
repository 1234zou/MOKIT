import h5py

def dump_scf_no_mol(chkfile, e_tot, mo_energy, mo_coeff, mo_occ):
    fh5 = h5py.File(chkfile, 'w')
    scf = fh5.create_group('scf')
    write_dic(scf, {'e_tot': e_tot,
                    'mo_energy': mo_energy,
                    'mo_occ': mo_occ,
                    'mo_coeff': mo_coeff})
    fh5.close()

def write_dic(root, dic):
    for k, v in dic.items():
        if isinstance(v, dict):
            grp = root.create_group(k)
            write_dic(grp, v)
        else:
            root[k] = v
