# Transfer MOs from PySCF to AMESP

def py2amesp(mf, aipname):
    from mokit.lib.py2fch_direct import fchk
    from os import system, remove
    fchname = aipname[0:aipname.rindex('.aip')]+'.fch'
    fchk(mf, fchname)
    system('fch2amo '+fchname)
    remove(fchname)

