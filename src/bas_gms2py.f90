! written by jxzou at 20180814
! This is a subroutine used to transform the GAMESS format of basis sets to the
!  NWChem format of basis sets, which can be used as input file for PySCF

! updated by jxzou at 20180815: Effective Core Potential (ECP) supported
! updated by jxzou at 20181212: fix the bug in ECP
! updated by jxzou at 20200324: simplify code
! updated by jxzou at 20200506: fix the bug of user-defined basis set; add content of fch2py
! updated by jxzou at 20210207: simplify code
! updated by jxzou at 20210527: remove string Sdiag

! Note: 1) Only basis set data and ECP data in GAMESS format are supported!
!       2) Isotopes are not supported so far!
program main
 implicit none
 integer :: i
 character(len=4) :: buf
 character(len=240) :: fname = ' '
 logical :: cart

 i = iargc()
 if(.not. (i==1 .or. i==2)) then
  write(6,'(/,A)') ' ERROR in subroutine bas_gms2py: wrong command line arguments!'
  write(6,'(A)')   ' Example 1: bas_gms2py a.inp '
  write(6,'(A,/)') ' Example 2: bas_gms2py a.inp -sph'
  stop
 end if

 call getarg(1, fname)
 call require_file_exist(fname)

 cart = .true.
 if(i == 2) then
  call getarg(2, buf)
  if(buf /= '-sph') then
   write(6,'(/,A)') ' ERROR in subroutine bas_gms2py: wrong command line arguments!'
   write(6,'(A)') "The second argument must be '-sph'."
   stop
  else
   cart = .false.
  end if
 end if

 call bas_gms2py(fname, cart)
end program main

! transform the GAMESS format of basis sets to the NWChem format,
! which can be used as input file for PySCF
subroutine bas_gms2py(inpname, cart)
 implicit none
 integer :: i, j, k, m, p, lmax, charge, mult
 integer :: nline, ncol, natom, nif, nbf
 integer :: inpid, pyid
 integer, allocatable :: ntimes(:) ! number of times of an atom appears
 integer, allocatable :: nuc(:)
 integer, external :: detect_ncol_in_buf
 real(kind=8) :: so_coeff
 real(kind=8), allocatable :: coor(:,:), prim_gau(:,:)
 character(len=1) :: bastype
 character(len=2) :: new_bastype
 character(len=2), allocatable :: elem(:)
 character(len=240) :: buf, pyname
 character(len=240), intent(in) :: inpname
 character(len=1), parameter :: am_type(0:6) = ['S','P','D','F','G','H','I']
 logical :: ecp, uhf, ghf, lin_dep
 logical, intent(in) :: cart

 buf = ' '
 pyname = ' '
 bastype = ' '
 new_bastype = ' '
 i = INDEX(inpname, '.', back=.true.)
 pyname = inpname(1:i-1)//'.py'

 call read_charge_and_mult_from_gms_inp(inpname, charge, mult, uhf, ghf, ecp)
 call read_nbf_and_nif_from_gms_inp(inpname, nbf, nif)
 if(nbf > nif) then
  lin_dep = .true.
 else if(nbf == nif) then
  lin_dep = .false.
 else
  write(6,'(A)') 'ERROR in ERROR in subroutine bas_gms2py: nbf<nif.'
  write(6,'(2(A,I0))') 'nbf=', nbf, ', nif=', nif
  stop
 end if

 call read_natom_from_gms_inp(inpname, natom)
 allocate(elem(natom), nuc(natom), coor(3,natom))
 call read_elem_nuc_coor_from_gms_inp(inpname, natom, elem, nuc, coor)
 deallocate(nuc)

 ! find the number of times of each atom occurs
 allocate(ntimes(natom))
 call calc_ntimes(natom, elem, ntimes)

 ! print the coordinates (optional) into .nwc file
 open(newunit=pyid,file=TRIM(pyname),status='replace')
 write(pyid,'(A)') 'from pyscf import gto, scf'
 write(pyid,'(A)') 'from mokit.lib.fch2py import fch2py'
 if(ghf) then
  write(pyid,'(A)') 'from mokit.lib.ortho import check_cghf_orthonormal'
 else
  write(pyid,'(A)') 'from mokit.lib.ortho import check_orthonormal'
 end if
 write(pyid,'(/,A)') 'mol = gto.M()'
 write(pyid,'(A,I0,A)') '# ',natom,' atom(s)'
 write(pyid,'(A)') "mol.atom = '''"

 do i = 1, natom, 1
  write(buf(1:6),'(A,I0)') TRIM(elem(i)), ntimes(i)
  write(pyid,'(A6,2X,3(1X,F17.8))') buf(1:6), coor(1:3,i)
 end do ! for i
 write(pyid,'(A)') "'''"
 deallocate(coor)
 ! print coordinates done

 open(newunit=inpid,file=TRIM(inpname),status='old',position='rewind')
 do while(.true.)
  read(inpid,'(A)') buf
  if(buf(2:2) == '$') then
   call upper(buf(3:6))
   if(buf(3:6)=='DATA') exit
  end if
 end do ! for while

 ! skip 2 lines: the Title Card line and the Point Group line
 read(inpid,'(A)') buf
 read(inpid,'(A)') buf

 ! read and print basis sets
 write(pyid,'(/,A)') 'mol.basis = {'
 do k = 1, natom, 1
  read(inpid,'(A)') buf
  ! 'buf' contains the element, relative atomic mass and coordinates
  if(buf(2:2) == '$') then
   call upper(buf(3:5))
   if(buf(3:5)=='END') exit
  end if

  ! deal with the basis data
  write(pyid,'(A,I0,A)') "'"//TRIM(elem(k)),ntimes(k),"': gto.basis.parse('''"

  do while(.true.)
   read(inpid,'(A)') buf
   if(LEN_TRIM(buf) == 0) exit

   read(buf,*) bastype, nline
   if(bastype == 'L') then
    ncol = 3
    new_bastype = 'SP'
   else
    ncol = 2
    new_bastype = bastype//' '
   end if
   allocate(prim_gau(ncol,nline), source=0d0)
   write(pyid,'(A2,4X,A2)') elem(k), new_bastype

   do i = 1, nline, 1
    read(inpid,*) j, prim_gau(1:ncol,i)
    if(nline==1 .and. prim_gau(2,1)<1d-6) prim_gau(2,1) = 1d0
    write(pyid,'(3(1X,E18.9))') prim_gau(:,i)
   end do ! for i
   deallocate(prim_gau)
  end do ! for while

  if(k < natom) then
   write(pyid,'(A)') "'''),"
  else
   write(pyid,'(A)') "''')}"
  end if
 end do ! k

 ! Note: it is required in .inp file that ECP is defined for every atom.
 !  For no ECP atom in .inp file, it must be assigned as 'ECP-NONE'
 if(ecp) then
  rewind(inpid)
  do while(.true.)
   read(inpid,'(A)') buf
   if(buf(2:2) == '$') then
    call upper(buf(3:5))
    if(buf(3:5) == 'ECP') exit
   end if
  end do ! for while

  write(pyid,'(/,A)') 'mol.ecp = {'

  do m = 1, natom, 1
   read(inpid,'(A)') buf
   if(index(buf,'NONE') /= 0) cycle

   write(pyid,'(A,I0,A)') "'"//TRIM(elem(m)),ntimes(m),"': gto.basis.parse_ecp('''"
   i = index(buf,'GEN')
   read(buf(i+3:),*) k, lmax
   write(pyid,'(A,I3)') TRIM(elem(m))//' nelec ', k

   allocate(prim_gau(2,1), source=0d0)
   do k = 1, lmax+1, 1
    read(inpid,'(A)') buf
    if(k == 1) then
     new_bastype = 'ul'
    else
     new_bastype = am_type(k-2)
    end if
    write(pyid,'(A)') TRIM(elem(m))//' '//TRIM(new_bastype)
    read(buf,*) nline
    do i = 1, nline, 1
     read(inpid,'(A)') buf
     p = detect_ncol_in_buf(buf)
     select case(p)
     case(3)
      read(buf,*) prim_gau(1,1), j, prim_gau(2,1)
      write(pyid,'(I1,1X,F17.10,3X,F17.10)') j, prim_gau(2,1), prim_gau(1,1)
     case(4)
      read(buf,*) prim_gau(1,1), j, prim_gau(2,1), so_coeff
      write(pyid,'(I1,1X,F17.10,2(3X,F17.10))') j, prim_gau(2,1), prim_gau(1,1),&
                                                so_coeff
     case default
      write(6,'(A,I0)') 'ERROR in subroutine bas_gms2py: invalid p=', p
      stop
     end select
    end do ! for i
   end do ! for k

   deallocate(prim_gau)
   write(pyid,'(A)') "'''),"
  end do ! for m

  BACKSPACE(pyid)
  write(pyid,'(A)') "''')}"
 end if ! for if(ecp)

 close(inpid)
 deallocate(elem, ntimes)

 write(pyid,'(/,A)') '# Remember to check the charge and spin'
 write(pyid,'(A,I0)') 'mol.charge = ', charge
 write(pyid,'(A,I0)') 'mol.spin = ', mult-1
 write(pyid,'(A)') 'mol.verbose = 4'
 if(cart) write(pyid,'(A)') 'mol.cart = True'
 write(pyid,'(A,/)') 'mol.build()'

 i = INDEX(inpname, '.', back=.true.)

 if(uhf) then
  if(lin_dep) then
   write(pyid,'(A)') 'old_mf = scf.UHF(mol)'
   write(pyid,'(A)') 'mf = scf.remove_linear_dep_(old_mf, threshold=1.1e-6, lindep=1.1e-6)'
  else
   write(pyid,'(A)') 'mf = scf.UHF(mol)'
  end if
  write(pyid,'(A)') 'mf.max_cycle = 1'
  write(pyid,'(A)') 'mf.kernel()'
  write(pyid,'(/,A)') '# read MOs from .fch(k) file'
  write(pyid,'(A)') 'nbf = mf.mo_coeff[0].shape[0]'
  write(pyid,'(A)') 'nif = mf.mo_coeff[0].shape[1]'
  write(pyid,'(A)') "alpha_coeff = fch2py('"//inpname(1:i-1)//".fch', nbf, nif, 'a')"
  write(pyid,'(A)') "beta_coeff = fch2py('"//inpname(1:i-1)//".fch', nbf, nif, 'b')"
  write(pyid,'(A)') 'mf.mo_coeff = (alpha_coeff, beta_coeff)'
  write(pyid,'(A)') '# read done'
  write(pyid,'(/,A)') '# check if input MOs are orthonormal'
  write(pyid,'(A)') "S = mol.intor_symmetric('int1e_ovlp')"
  write(pyid,'(A)') 'check_orthonormal(nbf, nif, mf.mo_coeff[0], S)'
  write(pyid,'(A)') 'check_orthonormal(nbf, nif, mf.mo_coeff[1], S)'

 else if(ghf) then ! complex GHF
  if(lin_dep) then
   write(pyid,'(A)') 'old_mf = scf.GHF(mol)'
   write(pyid,'(A)') 'mf = scf.remove_linear_dep_(old_mf, threshold=1.1e-6, lindep=1.1e-6)'
  else
   write(pyid,'(A)') 'mf = scf.GHF(mol)'
  end if
  if(ecp) then
   write(pyid,'(A)') 'mf.with_soc = True'
  else
   write(pyid,'(A)') '#mf.with_soc = True'
  end if
  write(pyid,'(A)') "dm = mf.get_init_guess(key='1e') + 0j"
  ! For super-heavy atoms, the initial guess 'minao' or 'atom' does not work, so
  !  I use '1e', which means initial guess from Hcore.
  ! Anyway, here the initial guess is not important since we will use fch2py to
  !  read MOs from a given .fch(k) file.
  write(pyid,'(A)') 'mf.max_cycle = 1'
  write(pyid,'(A)') 'mf.kernel(dm0=dm)'
  write(pyid,'(/,A)') '# read MOs from .fch(k) file'
  write(pyid,'(A)') 'nbf = mf.mo_coeff.shape[0]'
  write(pyid,'(A)') 'nif = mf.mo_coeff.shape[1]'
  write(pyid,'(A)') "mf.mo_coeff.real = fch2py('"//inpname(1:i-1)//".fch', nbf, nif, 'r')"
  write(pyid,'(A)') "mf.mo_coeff.imag = fch2py('"//inpname(1:i-1)//".fch', nbf, nif, 'i')"
  write(pyid,'(A)') '# read done'
  write(pyid,'(/,A)') '# check if input MOs are orthonormal'
  write(pyid,'(A)') "S = mol.intor_symmetric('int1e_ovlp')"
  write(pyid,'(A)') 'check_cghf_orthonormal(nbf, nif, mf.mo_coeff, S)'

 else ! using R(O)HF format
  if(mult == 1) then
   if(lin_dep) then
    write(pyid,'(A)') 'old_mf = scf.RHF(mol)'
    write(pyid,'(A)') 'mf = scf.remove_linear_dep_(old_mf, threshold=1.1e-6, lindep=1.1e-6)'
   else
    write(pyid,'(A)') 'mf = scf.RHF(mol)'
   end if
  else
   if(lin_dep) then
    write(pyid,'(A)') 'old_mf = scf.ROHF(mol)'
    write(pyid,'(A)') 'mf = scf.remove_linear_dep_(old_mf, threshold=1.1e-6, lindep=1.1e-6)'
   else
    write(pyid,'(A)') 'mf = scf.ROHF(mol)'
   end if
  end if
  write(pyid,'(A)') 'mf.max_cycle = 1'
  write(pyid,'(A)') 'mf.kernel()'
  write(pyid,'(/,A)') '# read MOs from .fch(k) file'
  write(pyid,'(A)') 'nbf = mf.mo_coeff.shape[0]'
  write(pyid,'(A)') 'nif = mf.mo_coeff.shape[1]'
  write(pyid,'(A)') "mf.mo_coeff = fch2py('"//inpname(1:i-1)//".fch', nbf, nif, 'a')"
  write(pyid,'(A)') '# read done'
  write(pyid,'(/,A)') '# check if input MOs are orthonormal'
  write(pyid,'(A)') "S = mol.intor_symmetric('int1e_ovlp')"
  write(pyid,'(A)') 'check_orthonormal(nbf, nif, mf.mo_coeff, S)'
 end if

 write(pyid,'(/,A)') '#dm = mf.make_rdm1()'
 write(pyid,'(A)')   '#mf.max_cycle = 10'
 write(pyid,'(A,/)') '#mf.kernel(dm)'
 close(pyid)
end subroutine bas_gms2py

