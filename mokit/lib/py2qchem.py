# Transfer MOs from PySCF to Q-Chem

def py2qchem(mf, inpname, npair=None):
    from mokit.lib.py2fch_direct import fchk
    from os import system, remove
    fchname = inpname[0:inpname.rindex('.in')]+'.fch'
    fchk(mf, fchname, overwrite_mol=True)
    if npair is None:
        system('fch2qchem '+fchname)
    else:
        system('fch2qchem '+fchname+' -gvb '+str(npair))
    remove(fchname)

