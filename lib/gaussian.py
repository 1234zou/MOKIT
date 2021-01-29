#!/usr/bin/env python
# written by jxzou at 20210129: subroutines involving Gaussian files

import os, shutil
from pyscf import gto

def load_mol_from_fch(fchname):
    proname = 'gau'+str(os.getpid())
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

