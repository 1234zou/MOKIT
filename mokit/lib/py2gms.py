# Transfer MOs from PySCF to GAMESS

def py2gms(mf, inpname, npair=None, nopen=None):
    from mokit.lib.py2fch_direct import fchk
    from os import system, remove
    fchname = inpname[0:inpname.rindex('.inp')]+'.fch'
    fchk(mf, fchname)
    if npair is None:
        system('fch2inp '+fchname)
    else:
        if nopen is None:
            system('fch2inp '+fchname+' -gvb '+str(npair))
        else:
            system('fch2inp '+fchname+' -gvb '+str(npair)+' -open '+str(nopen))
    remove(fchname)

