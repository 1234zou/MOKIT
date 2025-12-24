# Transfer MOs from PySCF -> Q-Chem

def py2qchem(mf, inpname, npair=None, sf=False, sasf=False, sfcis=False, sasfcis=False):
    from mokit.lib.py2fch_direct import fchk
    from os import system, remove

    fchname = inpname[0:inpname.rindex('.in')]+'.fch'
    fchk(mf, fchname, overwrite_mol=True)
    if npair is None:
        if sf:
            system('fch2qchem '+fchname+' -sf')
        elif sasf:
            system('fch2qchem '+fchname+' -sasf')
        elif sfcis:
            system('fch2qchem '+fchname+' -sfcis')
        elif sasfcis:
            system('fch2qchem '+fchname+' -sasfcis')
        else:
            system('fch2qchem '+fchname)
    else:
        system('fch2qchem '+fchname+' -gvb '+str(npair))
    remove(fchname)

