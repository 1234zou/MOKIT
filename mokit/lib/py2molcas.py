# Transfer MOs from PySCF -> (Open)Molcas

def py2molcas(mf, inpname, xc=None):
    from mokit.lib.py2fch_direct import fchk
    from os import system, remove
    fchname = inpname[0:inpname.rindex('.input')]+'.fch'
    fchk(mf, fchname, overwrite_mol=True)
    if xc is None:
        system('fch2inporb '+fchname)
    else:
        system('fch2inporb '+fchname+' -dft "'+xc+'"')
    remove(fchname)

