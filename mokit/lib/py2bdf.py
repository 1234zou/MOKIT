# Transfer MOs from PySCF -> BDF

def py2bdf(mf, inpname, write_no=False, xtda=False, socxtda=False, xtddft=False):
    from mokit.lib.py2fch_direct import fchk
    from os import system, remove, rename
    proname = inpname[0:inpname.rindex('.inp')]
    fchname = proname+'.fch'
    inpname1 = proname+'_bdf.inp'
    fchk(mf, fchname, overwrite_mol=True)
    if (write_no):
        orbname = proname+'.inporb'
        orbname1 = proname+'_bdf.inporb'
        system('fch2bdf '+fchname+' -no')
    else:
        orbname = proname+'.scforb'
        orbname1 = proname+'_bdf.scforb'
        if (xtda):
            system('fch2bdf '+fchname+' -xtda')
        elif (socxtda):
            system('fch2bdf '+fchname+' -socxtda')
        elif (xtddft):
            system('fch2bdf '+fchname+' -xtddft')
        else:
            system('fch2bdf '+fchname)
    remove(fchname)
    rename(inpname1, inpname)
    rename(orbname1, orbname)

