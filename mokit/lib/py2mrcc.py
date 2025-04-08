# Transfer MOs from PySCF to MRCC

def py2mrcc(mf):
    '''
    generate MINP, GENBAS and MOCOEF
    '''
    from mokit.lib.py2fch_direct import fchk
    from os import system, remove
    from random import randint
    fchname = str(randint(1,100000))+'.fch'
    fchk(mf, fchname)
    system('fch2mrcc '+fchname)
    remove(fchname)

