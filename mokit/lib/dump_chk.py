import h5py
import numpy as np

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

def dump_scf_for_rest(chkfile, e_tot, mo_energy, mo_coeff, mo_occ,
                      basis4elem_json=None, geom_json=None, cinttype='spheric',
                      num_basis=None, num_states=None, spin_channel=None,
                      spin=None, charge=None, num_elec=None):
    # Write a REST chkfile/guessfile. In addition to the 'scf' group, the new
    # REST chkfile format (>= 2026.1.1) requires the 'molecule' group carrying
    # 'basis4elem', 'geom', 'cinttype' and 'num_elec' as string scalars; see
    # fileop/chkfile.rs of REST.
    fh5 = h5py.File(chkfile, 'w')
    scf = fh5.create_group('scf')
    dic = {'e_tot': e_tot,
           'mo_energy': mo_energy,
           'mo_occ': mo_occ,
           'mo_coeff': mo_coeff}
    if num_basis is not None:
        dic['num_basis'] = np.array([num_basis], dtype=np.int64)
    if num_states is not None:
        dic['num_states'] = np.array([num_states], dtype=np.int64)
    if spin_channel is not None:
        dic['spin_channel'] = np.array([spin_channel], dtype=np.int64)
    if spin is not None:
        dic['spin'] = np.array([spin], dtype=np.float64)
    if charge is not None:
        dic['charge'] = np.array([charge], dtype=np.float64)
    write_dic(scf, dic)

    mol = fh5.create_group('molecule')
    sdt = h5py.string_dtype(encoding='utf-8')
    if basis4elem_json is not None:
        mol.create_dataset('basis4elem', data=basis4elem_json, dtype=sdt)
    if geom_json is not None:
        mol.create_dataset('geom', data=geom_json, dtype=sdt)
    if cinttype is not None:
        mol.create_dataset('cinttype', data=cinttype, dtype=sdt)
    if num_elec is not None:
        mol.create_dataset('num_elec', data=str(float(num_elec)), dtype=sdt)
    fh5.close()
