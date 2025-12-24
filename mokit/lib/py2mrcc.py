# Transfer MOs from PySCF -> MRCC

def py2mrcc(mf):
    '''
    generate MINP, GENBAS and MOCOEF
    '''
    from mokit.lib.py2fch_direct import fchk
    from os import system, remove
    from random import randint
    fchname = str(randint(1,100000))+'.fch'
    fchk(mf, fchname, overwrite_mol=True)
    system('fch2mrcc '+fchname)
    remove(fchname)

