# Transfer MOs from PySCF -> REST

def py2rest(mf, inpname, xc=None):
    from mokit.lib.py2fch_direct import fchk
    from os import system, remove

    fchname = inpname[0:inpname.rindex('.in')]+'.fch'
    fchk(mf, fchname, overwrite_mol=True)
    if xc is None:
        system('fch2rest '+fchname)
    else:
        system('fch2rest '+fchname+" -dft '"+xc+"'")
    remove(fchname)

