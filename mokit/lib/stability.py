#!/usr/bin/env python

# This file is firstly created at 20240701.
# The reasons to create this file:
# 1) PySCF use tol=1e-4 by default in the Davidson step for wave function stability
#    check, which is insufficient for some difficult systems (1e-7 required). The
#    tol value cannot be modified in the master branch currently. So jxzou copies
#    and modifies some content from pyscf/scf/stability.py to this file.
# 2) The `stable=opt` functionality is moved from mokit/src/do_hf.f90 to this file.

# This file is for temporarily usage, which may be removed in ~2 year (when most
# people are used to new PySCF, say >= v2.7). Forcing all users to update to the
# newest PySCF is not realistic currently.

# Related pull request (https://github.com/pyscf/pyscf/pull/2299).
# If any PySCF developers think this action is inappropriate, please contact
# MOKIT developers to delete this file.

from pyscf import lib
from pyscf.scf import hf
from pyscf.soscf import newton_ah
import numpy as np

def _rotate_mo(mo_coeff, mo_occ, dx):
  dr = hf.unpack_uniq_var(dx, mo_occ)
  u = newton_ah.expmat(dr)
  return np.dot(mo_coeff, u)

def rohf_internal(mf, tol=1e-7, verbose=None):
  log = lib.logger.new_logger(mf, verbose)
  g, hop, hdiag = newton_ah.gen_g_hop_rohf(mf, mf.mo_coeff, mf.mo_occ,
                                           with_symmetry=True)
  hdiag *= 2
  stable = True

  def precond(dx, e, x0):
    hdiagd = hdiag - e
    hdiagd[abs(hdiagd)<1e-8] = 1e-8
    return dx/hdiagd

  def hessian_x(x): # See comments in function rhf_internal
    return hop(x).real * 2

  x0 = np.zeros_like(g)
  x0[g!=0] = 1. / hdiag[g!=0]
  e, v = lib.davidson(hessian_x, x0, precond, tol=tol, verbose=log)
  if e < -1e-5:
    log.note(f'{mf.__class__} wavefunction has an internal instability.')
    mo = _rotate_mo(mf.mo_coeff, mf.mo_occ, v)
    stable = False
  else:
    log.note(f'{mf.__class__} wavefunction is stable in the internal '
             'stability analysis')
    mo = mf.mo_coeff
  return mo, stable

def uhf_internal(mf, tol=1e-7, verbose=None):
  log = lib.logger.new_logger(mf, verbose)
  g, hop, hdiag = newton_ah.gen_g_hop_uhf(mf, mf.mo_coeff, mf.mo_occ,
                                          with_symmetry=True)
  hdiag *= 2
  stable = True

  def precond(dx, e, x0):
    hdiagd = hdiag - e
    hdiagd[abs(hdiagd)<1e-8] = 1e-8
    return dx/hdiagd

  def hessian_x(x):
    return hop(x).real * 2

  x0 = np.zeros_like(g)
  x0[g!=0] = 1. / hdiag[g!=0]
  e, v = lib.davidson(hessian_x, x0, precond, tol=tol, verbose=log)
  if e < -1e-5:
    log.note(f'{mf.__class__} wavefunction has an internal instability.')
    nocca = np.count_nonzero(mf.mo_occ[0]> 0)
    nvira = np.count_nonzero(mf.mo_occ[0]==0)
    mo = (_rotate_mo(mf.mo_coeff[0], mf.mo_occ[0], v[:nocca*nvira]),
          _rotate_mo(mf.mo_coeff[1], mf.mo_occ[1], v[nocca*nvira:]))
    stable = False
  else:
    log.note(f'{mf.__class__} wavefunction is stable in the internal '
             'stability analysis')
    mo = mf.mo_coeff
  return mo, stable

def rohf_stable_opt_internal(mf):
  i = 0
  while (i < 10):
    i += 1
    mo, stable = rohf_internal(mf, verbose=5)
    if (stable):
      break
    else:
      dm = mf.make_rdm1(mo, mf.mo_occ)
      mf = mf.newton()
      mf.kernel(dm0=dm)
  if not stable:
    raise OSError('PySCF ROHF stable=opt failed after 10 attempts.')
  return mf

def uhf_stable_opt_internal(mf):
  i = 0
  while (i < 10):
    i += 1
    mo, stable = uhf_internal(mf, verbose=5)
    if (stable):
      break
    else:
      dm = mf.make_rdm1(mo, mf.mo_occ)
      mf = mf.newton()
      mf.kernel(dm0=dm)
  if not stable:
    raise OSError('PySCF UHF stable=opt failed after 10 attempts.')
  return mf

