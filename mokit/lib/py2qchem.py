# Transfer MOs from PySCF to Q-Chem

def py2qchem(mf, inpname, npair=None):
    from py2fch_direct import fchk
    from os import system, remove
    fchname = inpname[0:inpname.rindex('.in')]+'.fch'
    fchk(mf, fchname)
    if npair is None:
        system('fch2qchem '+fchname)
    else:
        system('fch2qchem '+fchname+' -gvb '+str(npair))
    remove(fchname)

