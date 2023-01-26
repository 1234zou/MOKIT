# Transfer MOs from PySCF to BDF

def py2bdf(mf, inpname, write_no=None):
    from py2fch_direct import fchk
    from os import system, remove, rename
    proname = inpname[0:inpname.rindex('.inp')]
    fchname = proname+'.fch'
    inpname1 = proname+'_bdf.inp'
    fchk(mf, fchname)
    if (write_no is None) or (write_no is False):
        orbname = proname+'.scforb'
        orbname1 = proname+'_bdf.scforb'
        system('fch2bdf '+fchname)
    elif write_no is True:
        orbname = proname+'.inporb'
        orbname1 = proname+'_bdf.inporb'
        system('fch2bdf '+fchname+' -no')
    else:
        raise AttributeError('write_no can only be None, True or False.')
    remove(fchname)
    rename(inpname1, inpname)
    rename(orbname1, orbname)

