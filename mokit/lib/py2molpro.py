# Transfer MOs from PySCF to Molpro

def py2molpro(mf, inpname):
    from py2fch_direct import fchk
    from os import system, remove
    fchname = inpname[0:inpname.rindex('.com')]+'.fch'
    fchk(mf, fchname)
    system('fch2com '+fchname)
    remove(fchname)

