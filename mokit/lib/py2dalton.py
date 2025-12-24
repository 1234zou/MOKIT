# Transfer MOs from PySCF -> Dalton

def py2dalton(mf, inpname):
    from mokit.lib.py2fch_direct import fchk
    from os import system, remove
    fchname = inpname[0:inpname.rindex('.dal')]+'.fch'
    fchk(mf, fchname, overwrite_mol=True)
    system('fch2dal '+fchname)
    remove(fchname)

