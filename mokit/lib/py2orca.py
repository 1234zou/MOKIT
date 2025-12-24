# Transfer MOs from PySCF -> ORCA
# Note: the Bohr constant is different between PySCF/MOKIT v.s. ORCA, so the
# difference of the nuclear repulsion energy may be as large as 1e-5 ~ 1e-6 a.u.
# which is usually larger than the difference of electronic energies. If you
# observe a 1e-5 ~ 1e-6 a.u. total energy difference between PySCF/ORCA, it is
# likely due to the Bohr constant. You can find the nuclear repulsion energy
# both in PySCF and ORCA output file, and make a comparison.

def py2orca(mf, inpname):
    from mokit.lib.py2fch_direct import fchk
    from os import system, remove, rename
    proname = inpname[0:inpname.rindex('.inp')]
    fchname = proname+'.fch'
    inpname1 = proname+'_o.inp'
    mklname1 = proname+'_o.mkl'
    inpname2 = proname+'.inp'
    mklname2 = proname+'.mkl'
    fchk(mf, fchname, overwrite_mol=True)
    system('fch2mkl '+fchname)
    remove(fchname)
    rename(inpname1, inpname2)
    rename(mklname1, mklname2)
    print('\n .mkl and .inp files are generated. Now do mkl->gbw...')
    i = system('orca_2mkl '+proname+' -gbw')
    if i != 0:
        raise OSError('Failed to call utility orca_2mkl')

