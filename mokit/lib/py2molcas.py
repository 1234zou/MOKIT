# Transfer MOs from PySCF to (Open)Molcas

def py2molcas(mf, inpname):
    from mokit.lib.py2fch_direct import fchk
    from os import system, remove
    fchname = inpname[0:inpname.rindex('.input')]+'.fch'
    fchk(mf, fchname)
    system('fch2inporb '+fchname)
    remove(fchname)

