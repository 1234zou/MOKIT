# Transfer MOs from PySCF -> PSI4

def py2psi(mf, inpname, xc=None):
    from mokit.lib.py2fch_direct import fchk
    from os import system, remove, rename
    proname = inpname[0:inpname.rindex('.inp')]
    fchname = proname+'.fch'
    inpname1 = proname+'_psi.inp'
    inpname2 = proname+'.inp'
    fchk(mf, fchname, overwrite_mol=True)
    if xc is None:
        system('fch2psi '+fchname)
    else:
        system('fch2psi '+fchname+' -dft "'+xc+'"')
    remove(fchname)
    rename(inpname1, inpname2)

