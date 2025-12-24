# Transfer MOs from PySCF -> CFOUR

def py2cfour(mf):
    '''
    generate ZMAT, GENBAS and OLDMOS
    '''
    from mokit.lib.py2fch_direct import fchk
    from os import system, remove
    from random import randint
    fchname = str(randint(1,10000))+'.fch'
    fchk(mf, fchname, overwrite_mol=True)
    system('fch2cfour '+fchname)
    remove(fchname)

