#!/usr/bin/env python
# written by jxzou at 20210129: subroutines involving Gaussian files

import random, os, shutil
import numpy as np
from mokit.lib.fch2py import fch2py
from mokit.lib.py2fch import py2fch
from mokit.lib.rwwfn import read_nbf_and_nif_from_fch
from mokit.lib.lo import boys, pm

def load_mol_from_fch(fchname):
  '''
  Load the PySCF mol object from a given Gaussian .fch(k) file

  Simple usage::
  >>> from pyscf import scf
  >>> from mokit.lib.gaussian import load_mol_from_fch
  >>> mol = load_mol_from_fch(fchname='benzene.fch')
  >>> mf = scf.RHF(mol).run()
  '''
  import importlib

  proname = 'gau'+str(random.randint(1,10000))
  tmp_fch = proname+'.fch'
  tmp_py  = proname+'.py'
  shutil.copyfile(fchname, tmp_fch)
  os.system('bas_fch2py '+tmp_fch)
  os.remove(tmp_fch)

  with open(tmp_py, 'r+') as fp:
    lines = fp.readlines()
    fp.seek(0)
    fp.truncate()
    fp.writelines('from pyscf import gto\n')
    for i, line in enumerate(lines):
      if ('mol.build' in line):
        j = i + 1
        break
    fp.writelines(lines[3:j])

  molpy = importlib.import_module(proname)
  importlib.invalidate_caches()
  # invalidate_caches() is needed, otherwise a second call of load_mol_from_fch
  # will probably lead to the error "ModuleNotFoundError: No module named 'gauxxx'"
  os.remove(tmp_py)
  shutil.rmtree('__pycache__')
  return molpy.mol


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


def loc(fchname, idx, method='pm', alpha=True):
  '''
  Perform orbital localization for a specified set of orbitals in a given
  Gaussian .fch(k) file.
  (The following 1e AO-basis integrals are computed using PySCF:
   1) overlap integrals for Pipek-Mezey localization;
   2) dipole integrals for Boys localization.)
  The method can be either 'pm' or 'boys'.

  Simple usage::
  >>> # perform Pipek-Mezey localization for occupied PI orbitals of benzene
  >>> # a file named benzene_rhf_LMO.fch will be created
  >>> from mokit.lib.gaussian import loc
  >>> loc(fchname='benzene_rhf.fch',idx=range(6,21))
  '''
  from pyscf.lo.boys import dipole_integral

  if alpha is True:
    spin = 'a'
  else:
    spin = 'b'
  fchname1 = fchname[0:fchname.rindex('.fch')]+'_LMO.fch'
  mol = load_mol_from_fch(fchname)
  nbf, nif = read_nbf_and_nif_from_fch(fchname)
  mo_coeff = fch2py(fchname, nbf, nif, spin)
  nmo = len(idx)

  if method == 'pm':
    S = mol.intor_symmetric('int1e_ovlp')
    loc_orb = pm(mol.nbas, mol._bas[:,0],mol._bas[:,1],mol._bas[:,3], mol.cart,
                 nbf, nmo, mo_coeff[:,idx], S, 'mulliken')
  elif method == 'boys':
    mo_dipole = dipole_integral(mol, mo_coeff[:,idx])
    loc_orb = boys(nbf, nmo, mo_coeff[:,idx], mo_dipole)
  else:
    raise ValueError('Localization method cannot be recognized.')

  mo_coeff[:,idx] = loc_orb.copy()
  noon = np.zeros(nif)
  shutil.copyfile(fchname, fchname1)
  py2fch(fchname1, nbf, nif, mo_coeff, spin, noon, False, False)
  print('Localized orbitals exported to file '+fchname1)


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
  import mokit.lib.uno as pyuno
  from mokit.lib.construct_vir import construct_vir
  from mokit.lib.rwwfn import read_na_and_nb_from_fch

  os.system('fch_u2r '+fchname)
  fchname0 = fchname[0:fchname.rindex('.fch')]+'_r.fch'
  fchname1 = fchname[0:fchname.rindex('.fch')]+'_UNO.fch'
  os.rename(fchname0, fchname1)
  nbf, nif = read_nbf_and_nif_from_fch(fchname)
  na, nb = read_na_and_nb_from_fch(fchname)
  alpha_mo = fch2py(fchname, nbf, nif, 'a')
  beta_mo  = fch2py(fchname, nbf, nif, 'b')
  mol = load_mol_from_fch(fchname)
  S = mol.intor_symmetric('int1e_ovlp')
  idx, noon, alpha_coeff = pyuno.uno(nbf, nif, na, nb, alpha_mo, beta_mo, S, 1e-5)
  alpha_coeff = construct_vir(nbf, nif, idx[1], alpha_coeff, S)
  os.remove('uno.out')
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
  cons = math.sqrt(2e0)
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


def gen_fcidump(fchname, nacto, nacte, mem=4000, np=None):
  '''
  generate a FCIDUMP file using the provided .fch(k) file
  nacto: the number of active orbitals
  nacte: the number of active electrons
  mem: total memory, in MB
  np: the number of OpenMP threads
  '''
  from pyscf import scf, mcscf, ao2mo, lib
  from pyscf.tools.fcidump import from_integrals

  if (np):
    lib.num_threads(np)

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
  from mokit.lib.rwwfn import read_na_and_nb_from_fch
  from mokit.lib.rwgeom import read_natom_from_fch, read_coor_from_fch
  from mokit.lib.mirror_wfn import rotate_atoms_wfn2
  from mokit.lib.mo_svd import orb_resemble
  from mokit.lib.qchem import read_hf_type_from_fch

  if nmo is None:
    nmo, nb = read_na_and_nb_from_fch(ref_fch)
    # set nmo as the number of alpha occupied orbitals

  if align is True:
    natom = read_natom_from_fch(target_fch)
    coor = read_coor_from_fch(target_fch, natom)
    ref_fch1 = ref_fch[0:ref_fch.rindex('.fch')]+'_rot.fch'
    rotate_atoms_wfn2(ref_fch, natom, coor, ref_fch1)
  else:
    ref_fch1 = ref_fch

  mol1 = load_mol_from_fch(target_fch)
  mol2 = load_mol_from_fch(ref_fch1)
  cross_S = gto.intor_cross('int1e_ovlp', mol1, mol2)

  nbf1, nif1 = read_nbf_and_nif_from_fch(target_fch)
  mo1 = fch2py(target_fch, nbf1, nif1, 'a')
  nbf2, nif2 = read_nbf_and_nif_from_fch(ref_fch1)
  mo2 = fch2py(ref_fch1, nbf2, nif2, 'a')

  # rotate alpha MOs of target molecule at target basis to resemble known orbitals
  mo3 = orb_resemble(nbf1, nif1, mo1, nbf2, nmo, mo2[:,0:nmo], cross_S)
  noon = np.zeros(nif1)
  py2fch(target_fch, nbf1, nif1, mo3, 'a', noon, False, False)

  # If UHF, rotate beta MOs, too
  ihf = read_hf_type_from_fch(ref_fch1)
  if ihf==1 or ihf==101: # real R(O)HF
    if align is True:
      os.remove(ref_fch1)
  elif ihf == 2: # UHF
    if nmo is None:
      na, nmo = read_na_and_nb_from_fch(ref_fch1)
      # set nmo as the number of beta occupied orbitals
    mo1 = fch2py(target_fch, nbf1, nif1, 'b')
    mo2 = fch2py(ref_fch1, nbf2, nif2, 'b')
    if align is True:
      os.remove(ref_fch1)
    mo3 = orb_resemble(nbf1, nif1, mo1, nbf2, nmo, mo2[:,0:nmo], cross_S)
    py2fch(target_fch, nbf1, nif1, mo3, 'b', noon, False, False)
  else:
    raise NotImplementedError('make_orb_resemble supports only R(O)HF/UHF currently.')


def proj2target_basis(fchname, target_basis='cc-pVTZ', nmo=None, cart=False):
  '''
  Project MOs of the original basis set onto the target basis set.
  cart: True/False for Cartesian-type or spherical harmonic type functions
  '''
  from pyscf import scf
  from mokit.lib.qchem import read_hf_type_from_fch
  from mokit.lib.rwwfn import get_no_from_density_and_ao_ovlp
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
  nif = nbf
  # we cannot check linear dependency of basis functions here, so set nif to nbf

  dm0 = mf.get_init_guess(mol, '1e')
  if ihf == 1:   # real RHF
    mf.mo_energy, mf.mo_coeff = get_no_from_density_and_ao_ovlp(nbf, nif, dm0, S)
  elif ihf == 2: # UHF
    dm0 = dm0[0] + dm0[1]
    mo_e_a, alpha_mo = get_no_from_density_and_ao_ovlp(nbf, nif, dm0, S)
    mf.mo_energy = (mo_e_a, mo_e_a)
    mf.mo_coeff = (alpha_mo, alpha_mo)
  if ihf == 101: # real ROHF
    dm0 = dm0[0] + dm0[1]
    mf.mo_energy, mf.mo_coeff = get_no_from_density_and_ao_ovlp(nbf, nif, dm0, S)

  target_fch = fchname[0:fchname.rindex('.fch')]+'_proj.fch'
  fchk(mf, target_fch)
  make_orb_resemble(target_fch, fchname, nmo=nmo, align=False)


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
  from mokit.lib.qchem import read_hf_type_from_fch
  from mokit.lib.rwwfn import read_na_and_nb_from_fch
  from mokit.lib.mirror_wfn import mo_grassmann_intrplt
  from mokit.lib.rwgeom import replace_coor_in_fch_by_gjf
  from mokit.lib.construct_vir import construct_vir

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

