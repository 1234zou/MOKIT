#!/usr/bin/env python
# written by jxzou at 20210129: subroutines involving Gaussian files

import random, os, shutil
import numpy as np
from pyscf import gto
from pyscf.lo.boys import dipole_integral
from fch2py import fch2py
from py2fch import py2fch
from rwwfn import read_nbf_and_nif_from_fch
from lo import boys, pm

def load_mol_from_fch(fchname):
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

