# Transfer MOs from PySCF to PSI4

def py2psi(mf, inpname):
    from mokit.lib.py2fch_direct import fchk
    from os import system, remove, rename
    proname = inpname[0:inpname.rindex('.inp')]
    fchname = proname+'.fch'
    inpname1 = proname+'_psi.inp'
    inpname2 = proname+'.inp'
    fchk(mf, fchname, overwrite_mol=True)
    system('fch2psi '+fchname)
    remove(fchname)
    rename(inpname1, inpname2)

