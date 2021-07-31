#!/usr/bin/env python
# written by jxzou at 20210129: subroutines involving Gaussian files

import random, os, shutil
import numpy as np
from pyscf import gto
from pyscf.lo.boys import dipole_integral
from fch2py import fch2py
from py2fch import py2fch
from rwwfn import read_nbf_and_nif_from_fch, read_na_and_nb_from_fch
import uno as pyuno
from construct_vir import construct_vir
from lo import boys, pm

def load_mol_from_fch(fchname):
  '''
  Load the PySCF mol object from a given Gaussian .fch(k) file

  Simple usage::
  >>> from pyscf import scf
  >>> from gaussian import load_mol_from_fch
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

def loc(fchname, method, idx):
  '''
  Perform orbital localization for a specified set of orbitals in a given
  Gaussian .fch(k) file.
  (The following 1e AO-basis integrals are computed using PySCF:
   1) overlap integrals for Pipek-Mezey localization;
   2) dipole integrals for Boys localization.)

  Simple usage::
  >>> # perform Pipek-Mezey localization for valence occupied orbitals of benzene
  >>> # a file named benzene_rhf_LMO.fch will be created
  >>> from gaussian import loc
  >>> loc(fchname='benzene_rhf.fch',method='pm',idx=range(6,21))
  '''
  fchname1 = fchname[0:fchname.rindex('.fch')]+'_LMO.fch'
  mol = load_mol_from_fch(fchname)
  nbf, nif = read_nbf_and_nif_from_fch(fchname)
  mo_coeff = fch2py(fchname, nbf, nif, 'a')
  nmo = len(idx)

  if(method == 'boys'):
    mo_dipole = dipole_integral(mol, mo_coeff[:,idx])
    loc_orb = boys(nbf, nmo, mo_coeff[:,idx], mo_dipole)
  elif(method == 'pm'):
    S = mol.intor_symmetric('int1e_ovlp')
    loc_orb = pm(mol.nbas,mol._bas[:,0],mol._bas[:,1],mol._bas[:,3],mol.cart,nbf,nmo,mo_coeff[:,idx],S,'mulliken')
  else:
    raise ValueError("Localization method cannot be recognized or supported.")

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
  >>> from gaussian import uno
  >>> uno(fchname='benzene_uhf.fch')
  '''
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
  idx, noon, alpha_coeff = pyuno.uno(nbf, nif, na, nb, alpha_mo, beta_mo, S, 0.99999E+00)
  alpha_coeff = construct_vir(nbf, nif, idx[1], alpha_coeff, S)
  os.remove('uno.out')
  py2fch(fchname1, nbf, nif, alpha_coeff, 'a', noon, True)
  print('UNOs exported to file '+fchname1)
