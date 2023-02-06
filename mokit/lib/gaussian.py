#!/usr/bin/env python
# written by jxzou at 20210129: subroutines involving Gaussian files

import random, os, shutil
import numpy as np
from pyscf import gto
from pyscf.lo.boys import dipole_integral
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
  proname = 'gau'+str(random.randint(1,10000))
  tmp_fch = proname+'.fch'
  tmp_py  = proname+'.py'
  shutil.copyfile(fchname, tmp_fch)
  os.system('bas_fch2py '+tmp_fch)
  os.remove(tmp_fch)

  os.system("sed -i '1,3d' "+tmp_py)
  os.system("sed -i '1i from pyscf import gto' "+tmp_py)
  os.system("sed -i '/build/{p;:a;N;$!ba;d}' "+tmp_py)
  molpy = __import__(proname)
  os.remove(tmp_py)
  shutil.rmtree('__pycache__')
  return molpy.mol

def mo_fch2py(fchname):
  '''
  Read MOs from a given Gaussian .fch(k) file, and convert MOs for usage in PySCF

  Simple usage::
  >>> from mokit.lib.gaussian import mo_fch2py
  >>> mo = mo_fch2py('h2o.fch')
  '''
  from mokit.lib.rwwfn import read_nbf_and_nif_from_fch
  from mokit.lib.excited import check_uhf_in_fch
  from mokit.lib.fch2py import fch2py

  nbf, nif = read_nbf_and_nif_from_fch(fchname)
  stat = check_uhf_in_fch(fchname)

  if stat == 0:
    mo = fch2py(fchname, nbf, nif, 'a')
  elif stat == -1:
    mo_a = fch2py(fchname, nbf, nif, 'a')
    mo_b = fch2py(fchname, nbf, nif, 'b')
    mo = (mo_a, mo_b)
  else:
    raise ValueError("Neither R(O)HF nor UHF. Confused.")
  return mo

def loc(fchname, idx, method=None):
  '''
  Perform orbital localization for a specified set of orbitals in a given
  Gaussian .fch(k) file.
  (The following 1e AO-basis integrals are computed using PySCF:
   1) overlap integrals for Pipek-Mezey localization;
   2) dipole integrals for Boys localization.)
  The method can be 'pm' or 'boys'.

  Simple usage::
  >>> # perform Pipek-Mezey localization for occupied PI orbitals of benzene
  >>> # a file named benzene_rhf_LMO.fch will be created
  >>> from mokit.lib.gaussian import loc
  >>> loc(fchname='benzene_rhf.fch',idx=range(6,21))
  '''

  if method is None:
    method = 'pm'
  fchname1 = fchname[0:fchname.rindex('.fch')]+'_LMO.fch'
  mol = load_mol_from_fch(fchname)
  nbf, nif = read_nbf_and_nif_from_fch(fchname)
  mo_coeff = fch2py(fchname, nbf, nif, 'a')
  nmo = len(idx)

  if method == 'pm':
    S = mol.intor_symmetric('int1e_ovlp')
    loc_orb = pm(mol.nbas,mol._bas[:,0],mol._bas[:,1],mol._bas[:,3],mol.cart,nbf,nmo,mo_coeff[:,idx],S,'mulliken')
  elif method == 'boys':
    mo_dipole = dipole_integral(mol, mo_coeff[:,idx])
    loc_orb = boys(nbf, nmo, mo_coeff[:,idx], mo_dipole)
  else:
    raise ValueError("Localization method cannot be recognized.")

  mo_coeff[:,idx] = loc_orb.copy()
  noon = np.zeros(nif)
  shutil.copyfile(fchname, fchname1)
  py2fch(fchname1, nbf, nif, mo_coeff, 'a', noon, False)
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
  py2fch(fchname1, nbf, nif, alpha_coeff, 'a', noon, True)
  print('UNOs exported to file '+fchname1)

def permute_orb(fchname, orb1, orb2):
  '''
  Permute two orbitals in a given Gaussian .fch(k) file
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

  # call Gaussian to calculate dipole integrals and calculate electronic dipole
  e_dip = get_e_dipole_using_density_in_fch(fchname, itype)
  print(' Dipole moment from electrons (a.u.):', e_dip)

  # (total) electric dipole moment
  dipole = e_dip + n_dip
  print(' Dipole moment (a.u.):', dipole)
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

def make_orb_resemble(target_fch, ref_fch, nmo=None):
  '''
  make a set of target MOs resembles the reference MOs
  (different geometries and different basis set in two .fch files are allowed)
  target_fch: the .fch file which holds MOs to be updated
  ref_fch: the .fch file which holds reference MOs
  nmo: indices 1~nmo MOs in ref_fch will be set as reference MOs
  If nmo is not given, it will be set as na (number of alpha electrons)
  '''
  from mokit.lib.rwwfn import read_na_and_nb_from_fch
  from mokit.lib.mo_svd import orb_resemble
  if nmo is None:
    nmo, nb = read_na_and_nb_from_fch(ref_fch)

  mol1 = load_mol_from_fch(target_fch)
  mol2 = load_mol_from_fch(ref_fch)
  nbf1, nif1 = read_nbf_and_nif_from_fch(target_fch)
  nbf2, nif2 = read_nbf_and_nif_from_fch(ref_fch)
  # rotate MOs of target molecule at target basis to resemble known orbitals
  cross_S = gto.intor_cross('int1e_ovlp', mol1, mol2)
  mo1 = fch2py(target_fch, nbf1, nif1, 'a')
  mo2 = fch2py(ref_fch, nbf2, nif2, 'a')
  mo3 = orb_resemble(nbf1, nif1, mo1, nbf2, nmo, mo2[:,0:nmo], cross_S)
  noon = np.zeros(nif1)
  py2fch(target_fch, nbf1, nif1, mo3, 'a', noon, False)

