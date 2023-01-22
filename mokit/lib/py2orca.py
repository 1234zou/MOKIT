# Transfer MOs from PySCF to ORCA

def py2orca(mf, inpname):
    from py2fch_direct import fchk
    from os import system, remove, rename
    proname = inpname[0:inpname.rindex('.inp')]
    fchname = proname+'.fch'
    inpname1 = proname+'_o.inp'
    mklname1 = proname+'_o.mkl'
    inpname2 = proname+'.inp'
    mklname2 = proname+'.mkl'
    fchk(mf, fchname)
    system('fch2mkl '+fchname)
    remove(fchname)
    rename(inpname1, inpname2)
    rename(mklname1, mklname2)
    print('\n .mkl and .inp files are generated. Now do mkl->gbw...')
    i = system('orca_2mkl '+proname+' -gbw')
    if i != 0:
        raise OSError('Failed to call utility orca_2mkl')

