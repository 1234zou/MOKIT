# Transfer MOs from PySCF to GAMESS

def py2gms(mf, inpname, npair=None, nopen=None, sf=False, mrsf=False):
    from mokit.lib.py2fch_direct import fchk
    from os import system, remove
    fchname = inpname[0:inpname.rindex('.inp')]+'.fch'
    fchk(mf, fchname, overwrite_mol=True)
    if npair is None:
        if sf:
            system('fch2inp '+fchname+' -sf')
        elif mrsf:
            system('fch2inp '+fchname+' -mrsf')
        else:
            system('fch2inp '+fchname)
    else:
        if sf or mrsf:
            raise ValueError('npair cannot be used with sf=True or mrsf=True.')
        if nopen is None:
            system('fch2inp '+fchname+' -gvb '+str(npair))
        else:
            system('fch2inp '+fchname+' -gvb '+str(npair)+' -open '+str(nopen))
    remove(fchname)

