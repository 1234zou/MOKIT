# Transfer MOs from PySCF -> Dalton

def py2dalton(mf, inpname, xc=None):
    from mokit.lib.py2fch_direct import fchk
    from os import system, remove
    fchname = inpname[0:inpname.rindex('.dal')]+'.fch'
    fchk(mf, fchname, overwrite_mol=True)
    if xc is None:
        system('fch2dal '+fchname)
    else:
        system('fch2dal '+fchname+' -dft "'+xc+'"')
    remove(fchname)

