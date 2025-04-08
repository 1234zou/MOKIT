#!/usr/bin/env python
# written by jxzou at 20210129: subroutines involving Gaussian files

import random, os, shutil
import numpy as np
from mokit.lib.fch2py import fch2py
from mokit.lib.py2fch import py2fch
from mokit.lib.rwwfn import read_nbf_and_nif_from_fch, read_na_and_nb_from_fch

BOHR2ANG = 0.52917721092e0


def ao_dipole_int(mol):
  # mol can be either molecule or cell object
  numerator = np.einsum('z,zx->x', mol.atom_charges(), mol.atom_coords())
  charge_center = (numerator / mol.atom_charges().sum())
  with mol.with_common_origin(charge_center):
    if getattr(mol, 'pbc_intor', None):
      ao_dip = mol.pbc_intor('int1e_r', comp=3, hermi=1)
    else:
      ao_dip = mol.intor_symmetric('int1e_r', comp=3)
  return charge_center, ao_dip


def load_mol_from_fch(fchname):
  '''
  Load the PySCF mol object from a specified Gaussian .fch(k) file

  Simple usage::
  >>> from pyscf import scf
  >>> from mokit.lib.gaussian import load_mol_from_fch
  >>> mol = load_mol_from_fch(fchname='benzene.fch')
  >>> mf = scf.RHF(mol)
  >>> mf.kernel()
  '''
  import sys, importlib
  from mokit.lib.qchem import find_and_del_pyc

  proname = 'gau'+str(random.randint(1,10000))
  tmp_fch = proname+'.fch'
  tmp_py  = proname+'.py'
  shutil.copyfile(fchname, tmp_fch)
  with os.popen('bas_fch2py '+tmp_fch+' -obj') as run:
    null = run.read()
  os.remove(tmp_fch)

  importlib.invalidate_caches() # important
  molpy = importlib.import_module(proname)
  os.remove(tmp_py)
  find_and_del_pyc(proname, sys.version)
  # It is not appropriate to delete the whole `__pycache__` directory since
  # other process(es) might be using this directory. Only the .pyc file (e.g.
  # gau5141.cpython-39.pyc) would be deleted.
  return molpy.mol


def load_mol_from_molden(molden, program):
  '''
  Load the PySCF mol object from a specified .molden file. Be careful that
  .molden file does not have any ECP/PP data.

  Simple usage::
  >>> from pyscf import scf
  >>> from mokit.lib.gaussian import load_mol_from_molden
  >>> mol = load_mol_from_molden(molden='benzene.molden',program='orca')
  >>> mf = scf.RHF(mol).run()
  '''
  with os.popen('molden2fch '+molden+' -'+program.lower()) as run:
    null = run.read()
  fchname = molden[0:molden.rindex('.molden')]+'.fch'
  mol = load_mol_from_fch(fchname)
  return mol


def load_cell_from_fch(fchname):
  '''
  Load the PySCF cell object from a specified Gaussian .fch(k) file. This file
  is supposed to include the wave function of an isolated molecule. The true PBC
  Gaussian .fch file is not supported currently. Since the lattice vectors are
  unknown, they will be set to 50.0 A temporarily. One should set cell.a
  appropriately and call cell.build after using this function. See pbc_loc()
  below for an example.

  Simple usage::
  >>> from pyscf import scf
  >>> from mokit.lib.gaussian import load_cell_from_fch
  >>> cell = load_cell_from_fch(fchname='water64.fch')
  >>> cell.a = np.eye(3)*12.42
  >>> cell.build(parse_arg=False)
  >>> mf = scf.RHF(cell)
  >>> mf.kernel()
  '''
  import sys, importlib
  from mokit.lib.qchem import find_and_del_pyc

  proname = 'gau'+str(random.randint(1,10000))
  tmp_fch = proname+'.fch'
  tmp_py  = proname+'.py'
  shutil.copyfile(fchname, tmp_fch)
  with os.popen('bas_fch2py '+tmp_fch+' -pbc -obj') as run:
    null = run.read()
  os.remove(tmp_fch)

  importlib.invalidate_caches() # important
  cellpy = importlib.import_module(proname)
  os.remove(tmp_py)
  find_and_del_pyc(proname, sys.version)
  return cellpy.cell


def mo_fch2py(fchname):
  '''
  Read MOs from a given Gaussian .fch(k) file, and convert MOs for usage in PySCF

  Simple usage::
  >>> from mokit.lib.gaussian import mo_fch2py
  >>> mf.mo_coeff = mo_fch2py('h2o.fch')
  '''
  from mokit.lib.qchem import read_hf_type_from_fch
  from mokit.lib.fch2py import fch2py, fch2py_cghf

  nbf, nif = read_nbf_and_nif_from_fch(fchname)
  ihf = read_hf_type_from_fch(fchname)
  # 1/2/7/101 for real RHF, real UHF, complex GHF, real ROHF, respectively

  if ihf==1 or ihf==101: # real RHF, ROHF
    mo = fch2py(fchname, nbf, nif, 'a')
  elif ihf == 2:         # real UHF
    mo_a = fch2py(fchname, nbf, nif, 'a')
    mo_b = fch2py(fchname, nbf, nif, 'b')
    mo = (mo_a, mo_b)
  elif ihf == 7:         # complex GHF
    mo = fch2py_cghf(fchname, 2*nbf, 2*nif)
  else:
    raise ValueError('Confused HF_type.')
  return mo


def loc(fchname, idx, method='pm', alpha=True, center_xyz=None, ions_centers=False,
        dis_tol=17.0, conv_tol=1e-5):
  '''
  Perform orbital localization for a specified set of orbitals in a given
  Gaussian .fch(k) file.
  (The following 1e AO-basis integrals are computed using PySCF:
   1) overlap integrals for Pipek-Mezey localization;
   2) dipole integrals for Boys localization.)
  The method can be either 'pm' or 'boys'. This function/module is designed for
  isolated molecules. For PBC systems, please use pbc_loc() instead.
  dis_tol=17.0 A for Boys, water64 cluster
  dis_tol=26.5 A for Boys, water128 cluster
  dis_tol=24.5 A for PM, water128 cluster

  Simple usage::
  >>> # perform Pipek-Mezey localization for occupied PI orbitals of benzene
  >>> # a file named benzene_rhf_LMO.fch will be created
  >>> from mokit.lib.gaussian import loc
  >>> loc(fchname='benzene_rhf.fch',idx=range(6,21))
  '''
  import time
  from pyscf import gto
  from mokit.lib.lo import boys, pm
  from mokit.lib.rwgeom import periodic_table as pt
  from mokit.lib.rwgeom import read_elem_and_coor_from_fch

  t0 = time.perf_counter()
  if dis_tol < 0.1:
    raise ValueError('dis_tol must be a reasonable and positive float number')

  if conv_tol > 0.1:
    raise ValueError('conv_tol must be <= 0.1')
  elif conv_tol < 1e-8:
    raise ValueError('conv_tol must be >= 1e-8')

  if alpha is True:
    spin = 'a'
  else:
    spin = 'b'

  fchname1 = fchname[0:fchname.rindex('.fch')]+'_LMO.fch'
  mol = load_mol_from_fch(fchname)
  nbf, nif = read_nbf_and_nif_from_fch(fchname)
  mo_coeff = fch2py(fchname, nbf, nif, spin)
  nmo = len(idx)

  print('\nOrbital range:', idx)
  print(f"distance threshold (A): {dis_tol:.3f}")
  print(f"convergence threshold: {conv_tol:.7f}")
  print('nmo= %d' % nmo)

  natom = mol.natm
  dis = BOHR2ANG*gto.inter_distance(mol)
  bfirst = np.zeros(natom+1, dtype=np.int32)
  bfirst[0] = 1
  bfirst[1:] = mol.aoslice_by_atom()[:,3] + 1
  S = mol.intor_symmetric('int1e_ovlp')

  if method == 'pm':
    lmo = pm(natom, nbf, nmo, bfirst, dis, mo_coeff[:,idx], S, 'mulliken',
             dis_tol, conv_tol)
    if center_xyz is not None:
      center, ao_dip = ao_dipole_int(mol)
      center *= BOHR2ANG
      ao_dip *= BOHR2ANG
  elif method == 'boys':
    center, ao_dip = ao_dipole_int(mol)
    center *= BOHR2ANG
    ao_dip *= BOHR2ANG
    lmo = boys(natom, nbf, nmo, bfirst, dis, mo_coeff[:,idx], S, ao_dip,
               dis_tol, conv_tol)
  else:
    raise ValueError('Localization method cannot be recognized.')

  if center_xyz is not None:   # print LMO centers into xyz
    mo_center = np.einsum('ui,xuv,vi->xi', lmo, ao_dip, lmo, optimize=True)
    # np.einsum with optimize=True seems slightly faster than calc_ctdc_diag
    for i in range(nmo):
      mo_center[:,i] = mo_center[:,i] + center
    if ions_centers is True:
      k = natom + nmo
      elem = np.full(k, 'X ', dtype='U2')
      elem0, nuc, coor0, ch, mu = read_elem_and_coor_from_fch(fchname, natom)
      elem[:natom] = elem0
      coor = np.zeros((3,k))
      coor[:,:natom] = coor0
      coor[:,natom:] = mo_center
      pt.write_xyz(k, elem, coor, center_xyz, np.zeros((3,3)))
    else:
      elem = np.full(nmo, 'X ', dtype='U2')
      pt.write_xyz(nmo, elem, mo_center, center_xyz, np.zeros((3,3)))

  mo_coeff[:,idx] = lmo.copy()
  noon = np.zeros(nif)
  shutil.copyfile(fchname, fchname1)
  py2fch(fchname1, nbf, nif, mo_coeff, spin, noon, False, False)
  print('Localized orbitals exported to file '+fchname1)

  t1 = time.perf_counter()
  elapsed_time = t1 - t0
  print(f"Localization time(sec): {elapsed_time:.2f}")


def pbc_mo_center(cell, loc_orb):
  '''
  Calculate orbital centers/centroid of gamma-point localized orbitals
  '''
  from pyscf import lib
  from pyscf.pbc.dft import gen_grid

  nmo = loc_orb.shape[1]
  ngrids = np.prod(cell.mesh)
  blksize = 2400*56
  grids = gen_grid.UniformGrids(cell)
  # Note: grids are generated around the origin (0,0,0) by default
  # grids.coords.shape is (ngrids,3)
  # G*r
  Gr = np.dot(cell.reciprocal_vectors(), grids.coords.T)
  # e^(-iG*r)
  r = np.exp(-Gr*1j)
  # MO values on grids
  mo_on_grid = np.empty((ngrids,nmo))
  for ip0, ip1 in lib.prange(0, ngrids, blksize):
    ao = cell.pbc_eval_gto('GTOval', grids.coords[ip0:ip1])
    mo_on_grid[ip0:ip1] = np.dot(ao, loc_orb)

  # M_nn = < w_n | e^(-iG*r) | w_n >, size (nmo,3)
  dip_complex = np.einsum('gi,gi,g,xg->ix', mo_on_grid, mo_on_grid, grids.weights, r)

  # Any complex number z can be written as z = C*e^(i*t) = C*(cos(t) + i*sin(t)),
  # where C is a real number and t is the angle theta, so
  #  ln(z) = ln(C*e^(i*t)) = ln(C) + i*t
  #  Im(ln(z)) = t
  # Interestingly, the result does not depend on C and is just theta.
  # Im(In(M_nn)) is np.angle(dip_complex)
  # -Im(In(M_nn))/(2*PI) is
  dip = -np.angle(dip_complex)/(2*np.pi)

  # Now dip[i,:] are fractional coordinates of the i-th MO center (which can be
  # seen from subsequent np.dot(a, dip.T)). These fractional coordinates can be
  # positive/negative. Since the box is located in the subspace (x>0,y>0,z>0),
  # we may hope that all MO centers are also located in that subspace for a better
  # visual inspection. This can be done via translation by adding lattice vectors.
  # Multiple integer times of lattice vectors can be applied, but usually one a/b/c
  # is enough.
  dip[dip < 0] += 1.0
  # a = cell.lattice_vectors()
  mo_center = BOHR2ANG*np.dot(cell.lattice_vectors(), dip.T)
  return mo_center, dip_complex


def pbc_loc(molden, box, method='berry', wannier_xyz=None, ions_centers=False,
            projected=None, dis_tol=27.0, conv_tol=1e-5, save_lmo=False,
            old_fch=None):
  '''
  Perform orbital localization for a specified set of orbitals in a given
  CP2K .molden file. The method can be either 'berry' or 'pm'.
  dis_tol=17.0/27.0 are sufficient for water / SnO2-H2O, respectively.
  Current limitations:
   1) only gamma-point; 2) only CP2K molden; 3) only Alpha spin.
  Only cubic cells have been tested so far.
  dis_tol=10.0 A for Berry, water128 box
  dis_tol= 9.5 A for Berry, water512 box
  dis_tol=27.0 A for Berry, SnO2-H2O box

  Simple usage::
  >>> # perform Boys orbital localization for water64 box
  >>> from mokit.lib.gaussian import pbc_loc
  >>> pbc_loc('water64-MOS-1_0.molden',box=np.eye(3)*12.42)
  '''
  import time
  from pyscf.pbc.df.ft_ao import ft_aopair
  from mokit.lib.rwwfn import read_lat_vec_from_file
  from mokit.lib.rwgeom import periodic_table as pt
  from mokit.lib.rwgeom import read_elem_and_coor_from_fch
  from mokit.lib.lo import berry, boys, pm, calc_dis_mat_from_coor_pbc
  from mokit.lib.ortho import check_orthonormal

  t0 = time.perf_counter()
  if dis_tol < 0.1:
    raise ValueError('dis_tol must be a reasonable and positive float number')

  #if conv_tol > 0.1:
  #  raise ValueError('conv_tol must be <= 0.1')
  #elif conv_tol < 1e-8:
  #  raise ValueError('conv_tol must be >= 1e-8')

  # determine the lattice vectors
  if isinstance(box, np.ndarray):
    lat_vec = box
  elif isinstance(box, str):
    lat_vec = read_lat_vec_from_file(box)
  else:
    raise ValueError('datatype of box cannot be identified.')

  proname = molden[0:molden.rindex('.molden')]
  fchname = proname+'.fch'
  lmo_fch = proname+'_LMO.fch'
  if wannier_xyz is None:
    wannier_xyz = proname+'_wannier.xyz'

  if old_fch is None:
    with os.popen('molden2fch '+molden+' -cp2k') as run:
      null = run.read()
  else:
    fchname = old_fch

  cell = load_cell_from_fch(fchname)
  cell.a = lat_vec
  cell.build(parse_arg=False)
  nbf, nif = read_nbf_and_nif_from_fch(fchname)
  na, nmo = read_na_and_nb_from_fch(fchname)

  print('Lattice vectors\n', lat_vec)
  print(f"distance threshold (A): {dis_tol:.3f}")
  print(f"convergence threshold: {conv_tol:.7f}")
  print('nmo= %d' % nmo)

  lines = cell.atom.strip().split('\n')
  coor_s = [line.split()[1:4] for line in lines]
  coor = np.transpose(np.array(coor_s, dtype=float))
  natom = cell.natm
  # Currently `gto.inter_distance` cannot be used here, since it calculates the
  # inter-atomic distances of an isolated molecule.
  dis = calc_dis_mat_from_coor_pbc(natom, cell.a, coor)
  bfirst = np.zeros(natom+1, dtype=np.int32)
  bfirst[0] = 1
  bfirst[1:] = cell.aoslice_by_atom()[:,3] + 1

  # determine which atoms are chosen (default: all atoms)
  # LMOs centered on chosen atoms would be constructed.
  if projected is None:
    chosen = np.ones(natom, dtype=bool)
  else:
    chosen = np.zeros(natom, dtype=bool)
    chosen[projected] = True

  mo = fch2py(fchname, nbf, nif, 'a')
  S = cell.pbc_intor('int1e_ovlp', hermi=1)
  # do not use cell.intor_symmetric('int1e_ovlp') here, since it returns the AO
  # overlap of an isolated molecule

  if method == 'berry':
    # Note: ao_zdip is a (double) complex array with size (3,nmo,nmo)
    ao_zdip = ft_aopair(cell, Gv=cell.reciprocal_vectors())
    ao_zdip *= BOHR2ANG
    lmo = berry(natom, nbf, nmo, bfirst, dis, mo[:,:nmo], S, ao_zdip, dis_tol,
                conv_tol)
  elif method == 'boys':
    raise ValueError('Boys orbital localization not implemented yet.')
    center, ao_dip = ao_dipole_int(cell)
    ao_dip *= BOHR2ANG
    lmo = boys(natom, nbf, nmo, bfirst, dis, mo[:,:nmo], S, ao_dip, dis_tol,
               conv_tol)
  elif method == 'pm':
    lmo = pm(natom, nbf, nmo, bfirst, dis, mo[:,:nmo], S, 'mulliken', dis_tol,
             conv_tol)
  else:
    raise ValueError('Localization method cannot be recognized.')

  # print LMO centers into xyz
  if method != 'berry':
    ao_zdip = ft_aopair(cell, Gv=cell.reciprocal_vectors())
    # Theoretically, the array ao_zdip is in unit Bohr currently, and is supposed
    # to be multiplied by the constant BOHR2ANG. But actually there is no need
    # to do that since np.angle will lead to a dimensionless array, where the
    # constant BOHR2ANG will be cancelled.
  mo_zdip = np.einsum('ui,xuv,vi->xi', lmo, ao_zdip, lmo, optimize=True)
  mo_dip = -np.angle(mo_zdip)/(2*np.pi)
  mo_dip[mo_dip < 0] += 1.0
  mo_center = BOHR2ANG*np.dot(cell.lattice_vectors(), mo_dip)
  if ions_centers is True:
    k = natom + nmo
    elem = np.full(k, 'X ', dtype='U2')
    elem0, nuc, coor0, ch, mu = read_elem_and_coor_from_fch(fchname, natom)
    elem[:natom] = elem0
    coor = np.zeros((3,k))
    coor[:,:natom] = coor0
    coor[:,natom:] = mo_center
    pt.write_xyz(k, elem, coor, wannier_xyz, cell.a)
  else:
    elem = np.full(nmo, 'X ', dtype='U2')
    pt.write_xyz(nmo, elem, mo_center, wannier_xyz, cell.a)

  # update MOs and print them into .fch
  check_orthonormal(nbf, nmo, lmo, S)
  mo[:,:nmo] = lmo.copy()
  if save_lmo is True:
    noon = np.zeros(nif)
    shutil.copyfile(fchname, lmo_fch)
    py2fch(lmo_fch, nbf, nif, mo, 'a', noon, False, False)
    print('Localized orbitals exported to file '+lmo_fch)

  if old_fch is None:
    os.remove(fchname)
  t1 = time.perf_counter()
  elapsed_time = t1 - t0
  print(f"Localization time(sec): {elapsed_time:.1f}")


def uno(fchname):
  '''
  Generate UHF natural orbitals(UNOs) from a given Gaussian .fch(k) file
  (AO-basis overlap integrals are computed using PySCF)
  
  Simple usage::
  >>> # generate UNOs for a UHF wave function of benzene
  >>> # a file named benzene_uhf_UNO.fch will be created
  >>> from mokit.lib.gaussian import uno
  >>> uno(fchname='benzene_uhf.fch')
  '''
  from mokit.lib.uno import uno as uhf_no
  from mokit.lib.rwwfn import construct_vir

  os.system('fch_u2r '+fchname)
  fchname0 = fchname[0:fchname.rindex('.fch')]+'_r.fch'
  fchname1 = fchname[0:fchname.rindex('.fch')]+'_UNO.fch'
  outname = str(random.randint(1,10000))+'.out'
  os.rename(fchname0, fchname1)
  nbf, nif = read_nbf_and_nif_from_fch(fchname)
  na, nb = read_na_and_nb_from_fch(fchname)
  alpha_mo = fch2py(fchname, nbf, nif, 'a')
  beta_mo  = fch2py(fchname, nbf, nif, 'b')
  mol = load_mol_from_fch(fchname)
  S = mol.intor_symmetric('int1e_ovlp')
  idx, noon, alpha_coeff = uhf_no(outname,nbf,nif,na,nb,alpha_mo,beta_mo,S,1e-5)
  alpha_coeff = construct_vir(nbf, nif, idx[1], alpha_coeff, S)
  os.remove(outname)
  py2fch(fchname1, nbf, nif, alpha_coeff, 'a', noon, True, True)
  print('UNOs exported to file '+fchname1)


def permute_orb(fchname, orb1, orb2):
  '''
  Permute two orbitals in a given Gaussian .fch(k) file.
  Note: orb1/orb2 are in Fortran convention (starts from 1)
  '''
  from mokit.lib.rwwfn import read_mo_from_fch, write_mo_into_fch, \
    read_eigenvalues_from_fch, write_eigenvalues_to_fch

  nbf, nif = read_nbf_and_nif_from_fch(fchname)
  mo = read_mo_from_fch(fchname, nbf, nif, 'a')
  mo1 = mo[:,orb1-1].copy()
  mo2 = mo[:,orb2-1].copy()
  mo[:,orb1-1] = mo2.copy()
  mo[:,orb2-1] = mo1.copy()

  ev = read_eigenvalues_from_fch(fchname, nif, 'a')
  r = ev[orb1-1]
  ev[orb1-1] = ev[orb2-1]
  ev[orb2-1] = r

  write_mo_into_fch(fchname, nbf, nif, 'a', mo)
  write_eigenvalues_to_fch(fchname, nif, 'a', ev, True)


def lin_comb_two_mo(fchname, orb1, orb2):
  '''
  Perform root2/2 (mo1+mo2) and root2/2 (mo1-mo2) unitary transformation for two
   specified MOs in a Gaussian .fch file.
  When the sigma and pi orbitals of a multiple bond are mixed (banana bond), this
   can be used to make them separated.
  Note: orb1/orb2 are in Python convention (starts from 0)
  '''
  import math
  from mokit.lib.rwwfn import read_mo_from_fch, write_mo_into_fch

  nbf, nif = read_nbf_and_nif_from_fch(fchname)
  mo = read_mo_from_fch(fchname, nbf, nif, 'a')
  mo1 = mo[:,orb1].copy()
  mo2 = mo[:,orb2].copy()
  cons = 0.5*math.sqrt(2.0)
  mo[:,orb1] = cons*(mo1 + mo2)
  mo[:,orb2] = cons*(mo1 - mo2)
  write_mo_into_fch(fchname, nbf, nif, 'a', mo)


def get_dipole(fchname, itype=1):
  '''
  Calculate the dipole moment using density in .fch(k) file
  itype=1/3/5/7 for Total SCF/CI/MP2/CC Density. Default: itype=1
  '''
  from mokit.lib.lo import get_e_dipole_using_density_in_fch
  from mokit.lib.rwgeom import read_natom_from_fch, read_elem_and_coor_from_fch, \
                               get_nuc_dipole
  # calculate nuclear dipole
  natom = read_natom_from_fch(fchname)
  elem, nuc, coor, charge, mult = read_elem_and_coor_from_fch(fchname, natom)
  n_dip = get_nuc_dipole(natom, nuc, coor)
  print('\n Dipole moment from nuclear charges (a.u.):', n_dip)

  # call Gaussian to calculate dipole integrals and the electronic dipole
  e_dip = get_e_dipole_using_density_in_fch(fchname, itype)
  print(' Dipole moment from electrons (a.u.):', e_dip)

  # total electric dipole moment
  dipole = e_dip + n_dip
  print(' Dipole moment (a.u.):', dipole)
  print(' Dipole moment (Debye):', dipole*2.541746231)
  return dipole


def gen_fcidump(fchname, nacto, nacte, mem=4000):
  '''
  generate a FCIDUMP file using the provided .fch(k) file
  nacto: the number of active orbitals
  nacte: the number of active electrons
  mem: total memory, in MB
  '''
  from pyscf import scf, mcscf, ao2mo
  from pyscf.tools.fcidump import from_integrals

  # load the mol object from a given .fch(k) file
  mol = load_mol_from_fch(fchname)

  if mol.spin == 0:
    mf = scf.RHF(mol)
  else:
    mf = scf.ROHF(mol)

  # do 1-cycle R(O)HF to make necessary arrays allocated
  mf.max_cycle = 1
  mf.max_memory = mem
  mf.kernel()

  # read MOs from a given .fch(k) file
  nbf, nif = read_nbf_and_nif_from_fch(fchname)
  mf.mo_coeff = fch2py(fchname, nbf, nif, 'a')

  # generate integrals and create FCIDUMP
  mc = mcscf.CASCI(mf, nacto, nacte)
  eri_cas = mc.get_h2eff()
  eri_cas = ao2mo.restore(8, eri_cas, nacto)
  h1eff, ecore = mc.get_h1eff()
  int_file = fchname[0:fchname.rindex('.fch')]+'.FCIDUMP'
  from_integrals(int_file, h1eff, eri_cas, nacto, nacte, ecore, ms=mol.spin)


def make_orb_resemble(target_fch, ref_fch, nmo=None, align=False):
  '''
  make a set of target MOs resembles the reference MOs
  (Different basis set in two .fch files are allowed, but their geometries
   should be identical or very similar. If two geometries are in different
   orientations, remember to set align=True)
  target_fch: the .fch file which holds MOs to be updated
  ref_fch: the .fch file which holds reference MOs
  nmo: indices 1~nmo MOs in ref_fch will be labeled as reference MOs
  align: whether to align two molecules
  If nmo is not given, it will be set as na (the number of alpha electrons)
  '''
  from pyscf import gto
  from mokit.lib.rwgeom import read_natom_from_fch, read_coor_from_fch
  from mokit.lib.mirror_wfn import rotate_atoms_wfn2
  from mokit.lib.mo_svd import orb_resemble
  from mokit.lib.qchem import read_hf_type_from_fch

  nmo_given = True
  if nmo is None:
    nmo_given = False
    nmo, nb = read_na_and_nb_from_fch(ref_fch)
    # nmo default: the number of alpha occupied orbitals

  if align is True:
    natom = read_natom_from_fch(target_fch)
    coor = read_coor_from_fch(target_fch, natom)
    ref_fch1 = ref_fch[0:ref_fch.rindex('.fch')]+'_rot.fch'
    rotate_atoms_wfn2(ref_fch, natom, coor, ref_fch1)
  else:
    ref_fch1 = ref_fch

  mol1 = load_mol_from_fch(target_fch)
  mol2 = load_mol_from_fch(ref_fch1)
  ao_S1 = mol1.intor_symmetric('int1e_ovlp')
  cross_S = gto.intor_cross('int1e_ovlp', mol1, mol2)

  nbf1, nif1 = read_nbf_and_nif_from_fch(target_fch)
  nbf2, nif2 = read_nbf_and_nif_from_fch(ref_fch1)
  mo2 = fch2py(ref_fch1, nbf2, nif2, 'a')

  # rotate alpha MOs of target molecule at target basis to resemble known orbitals
  mo1 = orb_resemble(nbf1, nif1, nbf2, nmo, mo2[:,:nmo], ao_S1, cross_S)
  noon = np.zeros(nif1)
  py2fch(target_fch, nbf1, nif1, mo1, 'a', noon, False, False)

  # If UHF, rotate beta MOs, too
  ihf = read_hf_type_from_fch(ref_fch1)
  if ihf==1 or ihf==101: # real R(O)HF
    if align is True:
      os.remove(ref_fch1)
  elif ihf == 2: # UHF
    if nmo_given is False:
      na, nmo = read_na_and_nb_from_fch(ref_fch1)
      # reset nmo as the number of beta occupied orbitals
    mo2 = fch2py(ref_fch1, nbf2, nif2, 'b')
    if align is True:
      os.remove(ref_fch1)
    mo1 = orb_resemble(nbf1, nif1, nbf2, nmo, mo2[:,:nmo], ao_S1, cross_S)
    py2fch(target_fch, nbf1, nif1, mo1, 'b', noon, False, False)
  else:
    raise NotImplementedError('Unsupported HF type in make_orb_resemble')


def proj2target_basis(fchname, target_basis='cc-pVTZ', nmo=None, cart=False):
  '''
  Project MOs of the original basis set onto the target basis set.
  cart: True/False for Cartesian-type or spherical harmonic type functions
  '''
  from pyscf import scf
  from mokit.lib.qchem import read_hf_type_from_fch
  from mokit.lib.rwwfn import gen_no_from_density_and_ao_ovlp
  from mokit.lib.lo import get_nmo_from_ao_ovlp
  from mokit.lib.py2fch_direct import fchk

  mol = load_mol_from_fch(fchname)
  mol.basis = target_basis
  mol.cart = cart
  mol.build(parse_arg=False)

  ihf = read_hf_type_from_fch(fchname)
  if ihf == 1:     # real RHF
    mf = scf.RHF(mol)
  elif ihf == 2:   # UHF
    mf = scf.UHF(mol)
  elif ihf == 101: # real ROHF
    mf = scf.ROHF(mol)
  else:
    raise NotImplementedError('proj2target_basis supports only R(O)HF/UHF currently.')

  S = mol.intor_symmetric('int1e_ovlp')
  nbf = S.shape[0]
  nif = get_nmo_from_ao_ovlp(nbf, S)
  if nif < nbf:
    mf2 = mf.copy()
    mf = scf.remove_linear_dep_(mf2, threshold=1.1e-6, lindep=1.1e-6)

  dm0 = mf.get_init_guess(mol, '1e')
  if ihf == 1:   # real RHF
    mf.mo_energy, mf.mo_coeff = gen_no_from_density_and_ao_ovlp(nbf, nif, dm0, S)
  elif ihf == 2: # UHF
    dm0 = dm0[0] + dm0[1]
    mo_e_a, alpha_mo = gen_no_from_density_and_ao_ovlp(nbf, nif, dm0, S)
    mf.mo_energy = (mo_e_a, mo_e_a)
    mf.mo_coeff = (alpha_mo, alpha_mo)
  if ihf == 101: # real ROHF
    dm0 = dm0[0] + dm0[1]
    mf.mo_energy, mf.mo_coeff = gen_no_from_density_and_ao_ovlp(nbf, nif, dm0, S)
  target_fch = fchname[0:fchname.rindex('.fch')]+'_proj.fch'
  fchk(mf, target_fch)
  make_orb_resemble(target_fch, fchname, nmo=nmo, align=False)


def mo_svd_in_fch(fchname1, fchname2, idx1=None, idx2=None):
  '''
  Perform SVD on two sets of MOs in two .fch(k) files.
  idx1/idx2: the 1st/last index of the MO, starts from 0
  '''
  from mokit.lib.wfn_analysis import mo_svd_in2fch
  from mokit.lib.rwwfn import read_nif_from_fch

  if idx1 is None:
    idx1 = 0
  if idx2 is None:
    idx2 = read_nif_from_fch(fchname1)
  print('idx1= %d, idx2= %d' %(idx1, idx2) )
  mo_svd_in2fch(fchname1, fchname2, idx1+1, idx2)


def export_mo_e2txt(fchname):
  '''
  export the data of Alpha Orbital Energies in a .fch file into a plain text file
  '''
  from mokit.lib.rwwfn import read_eigenvalues_from_fch, export_rarray2txt

  txtname = fchname[0:fchname.rindex('.fch')]+'.txt'
  nbf, nif = read_nbf_and_nif_from_fch(fchname)
  ev = read_eigenvalues_from_fch(fchname, nif, 'a')
  export_rarray2txt(txtname, 'MO Eigenvalues', nif, ev)


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

  # check the numeber of alpha electrons
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
    # check the numeber of beta electrons
    if trace_PS is True:
      dm = np.dot(coeff[:,:nb], coeff[:,:nb].transpose())
      ne_b = np.trace(np.dot(dm,S[:,:,nfile-1]))
      print('No. beta electrons: %.4f' %ne_b)

