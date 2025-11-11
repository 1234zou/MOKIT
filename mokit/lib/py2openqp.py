# Transfer MOs from PySCF to OpenQP

def py2openqp(mf, inpname, sfcis=False, sf=False, mrsfcis=False, mrsf=False):
    from mokit.lib.py2fch_direct import fchk
    from os import system, rename, remove

    fchname = inpname[0:inpname.rindex('.inp')]+'.fch'
    fchk(mf, fchname, overwrite_mol=True)

    if mf.mol.cart is False:
        system('fch_sph2cart '+fchname)
        remove(fchname)
        cart_fch = fchname[0:fchname.rindex('.fch')]+'_c.fch'
        rename(cart_fch, fchname)

    if sfcis:
        system('fch2openqp '+fchname+' -sfcis')
    elif sf:
        system('fch2openqp '+fchname+' -sf')
    elif mrsfcis:
        system('fch2openqp '+fchname+' -mrsfcis')
    elif mrsf:
        system('fch2openqp '+fchname+' -mrsf')
    else:
        system('fch2openqp '+fchname)

    remove(fchname)

