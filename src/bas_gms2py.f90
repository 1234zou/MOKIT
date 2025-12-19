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
 integer :: i, j
 character(len=5) :: buf
 character(len=240) :: inpname = ' '
 logical :: cart, pbc, obj_only

 i = iargc()
 if(i<1 .or. i>4) then
  write(6,'(/,A)') ' ERROR in program bas_gms2py: wrong command line arguments!'
  write(6,'(A)')   ' Example 1: bas_gms2py a.inp'
  write(6,'(A)')   ' Example 2: bas_gms2py a.inp -sph'
  write(6,'(A)')   ' Example 3: bas_gms2py a.inp -sph -pbc'
  write(6,'(A)')   ' Example 4: bas_gms2py a.inp -sph -obj'
  write(6,'(A)')   ' Example 5: bas_gms2py a.inp -obj'
  write(6,'(A,/)') ' Example 6: bas_gms2py a.inp -sph -pbc -obj'
  stop
 end if

 buf = ' '; inpname = ' '
 call getarg(1, inpname)
 call require_file_exist(inpname)

 cart = .true.; pbc = .false.; obj_only = .false.

 if(i > 1) then
  do j = 2, i, 1
   call getarg(j, buf)
   buf = ADJUSTL(buf)
   select case(TRIM(buf))
   case('-sph')
    cart = .false.
   case('-pbc')
    pbc = .true.
   case('-obj')
    obj_only = .true.
   case('-rest')
    write(6,'(/,A)') "ERROR in program bas_gms2py: '-rest' is no longer suppor&
                     &ted in bas_gms2py."
    stop
   case default
    write(6,'(/,A)') 'ERROR in program bas_gms2py: wrong command line arguments!'
    write(6,'(A)') "The arguments can only be -sph/-pbc/-obj"
    stop
   end select
  end do ! for j
 end if

 call bas_gms2py(inpname, cart, pbc, obj_only)
end program main

! transform the GAMESS format of basis sets to the NWChem format,
! which can be used as input file for PySCF
subroutine bas_gms2py(inpname, cart, pbc, obj_only)
 implicit none
 integer :: i, j, k, m, p, lmax, charge, mult, isph, inpid, pyid, order
 integer :: nline, ncol, natom, nif, nbf
 integer, allocatable :: ntimes(:) ! number of times of an atom appears
 integer, allocatable :: nuc(:)
 integer, external :: detect_ncol_in_buf
 real(kind=8) :: so_coeff
 real(kind=8), allocatable :: coor(:,:), prim_gau(:,:)
 character(len=1) :: bastype
 character(len=2) :: new_bastype
 character(len=2), allocatable :: elem(:)
 character(len=4) :: str4
 character(len=64), parameter :: lin_dep_str = &
  'mf = scf.remove_linear_dep_(old_mf, threshold=1e-6, lindep=1e-6)'
 character(len=240) :: buf, pyname, hf_fch
 character(len=240), intent(in) :: inpname
 character(len=1), parameter :: am_type(0:6) = ['S','P','D','F','G','H','I']
 logical :: ecp, uhf, ghf, X2C, lin_dep
 logical, intent(in) :: cart, pbc, obj_only
 logical, allocatable :: ghost(:) ! size natom

 str4 = 'mol '
 if(pbc) str4 = 'cell'
 buf = ' '; pyname = ' '; bastype = ' '; new_bastype = ' '

 i = INDEX(inpname, '.', back=.true.)
 if(i == 0) then
  write(6,'(/,A)') "ERROR in subroutine bas_gms2py: '.' symbol is required in f&
                   &ilename."
  write(6,'(A)') 'But the provided filename='//TRIM(inpname)
  stop
 end if
 pyname = inpname(1:i-1)//'.py'
 hf_fch = inpname(1:i-1)//'.fch'

 call read_charge_mult_isph_from_gms_inp(inpname, charge, mult, isph, uhf, ghf,&
                                         ecp)
 call read_nbf_and_nif_from_gms_inp(inpname, nbf, nif)

 ! We have no AO overlap here to detect the basis set linear dependency, so use
 ! nbf and nif.
 ! Some .fch(k) file only includes occupied MOs and thus nif<nbf does not neces-
 ! sarily means basis set linear dependency. Luckily, setting lin_dep=.true.
 ! would not lead to anything wrong.
 if(nbf > nif) then
  lin_dep = .true.
 else if(nbf == nif) then
  lin_dep = .false.
 else
  write(6,'(/,A)') 'ERROR in subroutine bas_gms2py: nbf<nif, impossible.'
  write(6,'(2(A,I0))') 'nbf=', nbf, ', nif=', nif
  stop
 end if

 call check_X2C_in_gms_inp(inpname, X2C)
 if(.not. X2C) then
  call check_DKH_in_gms_inp(inpname, order)
  if(order > -2) then
   write(6,'(/,A)') 'Warning from subroutine check_DKH_in_gms_inp: this type of &
                    &relativistic'
   write(6,'(A)') 'hamiltonian is not supported by PySCF. Anyway, file conversio&
                  &n will be proceeded.'
  end if
 end if
 call read_natom_from_gms_inp(inpname, natom)
 allocate(elem(natom), nuc(natom), coor(3,natom), ghost(natom))
 call read_elem_nuc_coor_from_gms_inp(inpname, natom, elem, nuc, coor, ghost)
 deallocate(nuc)

 ! find the number of times of each atom occurs
 allocate(ntimes(natom))
 call calc_ntimes(natom, elem, ntimes)

 open(newunit=pyid,file=TRIM(pyname),status='replace')
 if(obj_only) then
  write(pyid,'(A)') 'from pyscf import gto'
 else
  if(pbc) then
   write(pyid,'(A)') 'from pyscf import gto, lib'
   write(pyid,'(A)') 'from pyscf.pbc import scf'
  else
   write(pyid,'(A)') 'from pyscf import gto, scf, lib'
  end if
 end if

 if(pbc) then
  write(pyid,'(A)') 'from pyscf.pbc import gto as pgto'
  write(pyid,'(A)') 'import numpy as np'
 end if

 if(.not. obj_only) then
  write(pyid,'(A)') 'from mokit.lib.fch2py import fch2py'
  if(ghf) then
   write(pyid,'(A)') 'from mokit.lib.ortho import check_cghf_orthonormal'
  else
   write(pyid,'(A)') 'from mokit.lib.ortho import check_orthonormal'
  end if
  write(pyid,'(/,A)') 'lib.num_threads(2) # 2 CPU cores'
 else
  write(pyid,'(/)',advance='no')
 end if

 if(pbc) then
  write(pyid,'(A)') 'cell = pgto.Cell()'
 else
  write(pyid,'(A)') 'mol = gto.M()'
 end if
 write(pyid,'(A,I0,A)') '# ',natom,' atom(s)'
 write(pyid,'(A)') TRIM(str4)//".atom = '''"

 do i = 1, natom, 1
  if(ghost(i)) then
   write(buf,'(A,I0)') 'X-'//TRIM(elem(i)), ntimes(i)
  else
   write(buf,'(A,I0)') TRIM(elem(i)), ntimes(i)
  end if
  write(pyid,'(A,2X,3(1X,F19.9))') TRIM(buf), coor(1:3,i)
 end do ! for i
 write(pyid,'(A)') "'''"
 ! print coordinates done

 call goto_data_section_in_gms_inp(inpname, inpid)
 read(inpid,'(A)') buf
 read(inpid,'(A)') buf

 ! read and print basis sets
 write(pyid,'(/,A)') TRIM(str4)//'.basis = {'
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
   write(pyid,'(A2,4X,A)') elem(k), TRIM(new_bastype)

   do i = 1, nline, 1
    read(inpid,*) j, prim_gau(1:ncol,i)
    !if(nline==1 .and. prim_gau(2,1)<1d-6) prim_gau(2,1) = 1d0
    !if(nline == 1) then
    ! rtmp = prim_gau(2,1)
    ! if(rtmp>-1d-5 .and. rtmp<0d0) then
    !  prim_gau(2,1) = -1d0
    ! else if(rtmp>0d0 .and. rtmp<1d-5) then
    !  prim_gau(2,1) = 1d0
    ! end if
    !end if
    if(nline == 1) prim_gau(2,1) = DSIGN(1d0, prim_gau(2,1))
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

  write(pyid,'(/,A)') TRIM(str4)//'.ecp = {'

  do m = 1, natom, 1
   read(inpid,'(A)') buf
   if(index(buf,'NONE') /= 0) cycle

   write(pyid,'(A,I0,A)') "'"//TRIM(elem(m)),ntimes(m),"': gto.basis.parse_ecp('''"
   i = INDEX(buf,'GEN')
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
      write(6,'(/,A,I0)') 'ERROR in subroutine bas_gms2py: invalid p=', p
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
 deallocate(elem, ntimes, coor, ghost)

 write(pyid,'(/,A)') '# Remember to check the charge and spin'
 if(pbc) then
  write(pyid,'(A,I0)') TRIM(str4)//'.charge = 0'
 else
  write(pyid,'(A,I0)') TRIM(str4)//'.charge = ', charge
 end if
 write(pyid,'(A,I0)') TRIM(str4)//'.spin = ', mult-1
 if(pbc) then
  write(pyid,'(A)') "cell.pseudo = 'gth-pbe'"
  write(pyid,'(A)') 'cell.verbose = 3'
  write(pyid,'(A)') 'cell.a = np.eye(3)*100.0'
 else
  write(pyid,'(A)') 'mol.verbose = 4'
 end if
 if(cart) write(pyid,'(A)') TRIM(str4)//'.cart = True'

 if(.not. (pbc .and. obj_only)) then
  write(pyid,'(A,/)') TRIM(str4)//'.build(parse_arg=False)'
 end if
 if(obj_only) return ! no more to print

 if(uhf) then
  if(lin_dep) then
   write(pyid,'(A)') 'old_mf = scf.UHF('//TRIM(str4)//')'
   write(pyid,'(A)') lin_dep_str
  else
   write(pyid,'(A)') 'mf = scf.UHF('//TRIM(str4)//')'
  end if
  write(pyid,'(A)') 'mf.max_cycle = 1'
  write(pyid,'(A)') "mf.init_guess = '1e'"
  write(pyid,'(A)') 'mf.max_memory = 4000 # MB'
  if(X2C) write(pyid,'(A)') 'mf = mf.x2c1e()'
  write(pyid,'(A)') 'mf.kernel()'
  write(pyid,'(/,A)') '# read MOs from .fch(k) file'
  write(pyid,'(A)') "hf_fch = '"//TRIM(hf_fch)//"'"
  write(pyid,'(A)') 'nbf = mf.mo_coeff[0].shape[0]'
  write(pyid,'(A)') 'nif = mf.mo_coeff[0].shape[1]'
  write(pyid,'(A)') "alpha_coeff = fch2py(hf_fch, nbf, nif, 'a')"
  write(pyid,'(A)') "beta_coeff = fch2py(hf_fch, nbf, nif, 'b')"
  write(pyid,'(A)') 'mf.mo_coeff = (alpha_coeff, beta_coeff)'
  write(pyid,'(A)') '# read done'
  write(pyid,'(/,A)') '# check if input MOs are orthonormal'
  if(pbc) then
   write(pyid,'(A)') "S = cell.pbc_intor('int1e_ovlp',hermi=1)"
  else
   write(pyid,'(A)') "S = mol.intor_symmetric('int1e_ovlp')"
  end if
  write(pyid,'(A)') 'check_orthonormal(nbf, nif, mf.mo_coeff[0], S)'
  write(pyid,'(A)') 'check_orthonormal(nbf, nif, mf.mo_coeff[1], S)'

 else if(ghf) then ! complex GHF
  if(lin_dep) then
   write(pyid,'(A)') 'old_mf = scf.GHF('//TRIM(str4)//')'
   write(pyid,'(A)') lin_dep_str
  else
   write(pyid,'(A)') 'mf = scf.GHF('//TRIM(str4)//')'
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
  write(pyid,'(A)') 'mf.max_memory = 4000 # MB'
  if(X2C) write(pyid,'(A)') 'mf = mf.x2c1e()'
  write(pyid,'(A)') 'mf.kernel(dm0=dm)'
  write(pyid,'(/,A)') '# read MOs from .fch(k) file'
  write(pyid,'(A)') "hf_fch = '"//TRIM(hf_fch)//"'"
  write(pyid,'(A)') 'nbf = mf.mo_coeff.shape[0]'
  write(pyid,'(A)') 'nif = mf.mo_coeff.shape[1]'
  write(pyid,'(A)') "mf.mo_coeff.real = fch2py(hf_fch, nbf, nif, 'r')"
  write(pyid,'(A)') "mf.mo_coeff.imag = fch2py(hf_fch, nbf, nif, 'i')"
  write(pyid,'(A)') '# read done'
  write(pyid,'(/,A)') '# check if input MOs are orthonormal'
  if(pbc) then
   write(pyid,'(A)') "S = cell.pbc_intor('int1e_ovlp',hermi=1)"
  else
   write(pyid,'(A)') "S = mol.intor_symmetric('int1e_ovlp')"
  end if
  write(pyid,'(A)') 'check_cghf_orthonormal(nbf, nif, mf.mo_coeff, S)'

 else ! using R(O)HF format
  if(mult == 1) then
   if(lin_dep) then
    write(pyid,'(A)') 'old_mf = scf.RHF('//TRIM(str4)//')'
    write(pyid,'(A)') lin_dep_str
   else
    write(pyid,'(A)') 'mf = scf.RHF('//TRIM(str4)//')'
   end if
  else
   if(lin_dep) then
    write(pyid,'(A)') 'old_mf = scf.ROHF('//TRIM(str4)//')'
    write(pyid,'(A)') lin_dep_str
   else
    write(pyid,'(A)') 'mf = scf.ROHF('//TRIM(str4)//')'
   end if
  end if
  write(pyid,'(A)') 'mf.max_cycle = 1'
  write(pyid,'(A)') "mf.init_guess = '1e'"
  write(pyid,'(A)') 'mf.max_memory = 4000 # MB'
  if(X2C) write(pyid,'(A)') 'mf = mf.x2c1e()'
  write(pyid,'(A)') 'mf.kernel()'

  write(pyid,'(/,A)') '# read MOs from .fch(k) file'
  write(pyid,'(A)') "hf_fch = '"//TRIM(hf_fch)//"'"
  write(pyid,'(A)') 'nbf = mf.mo_coeff.shape[0]'
  write(pyid,'(A)') 'nif = mf.mo_coeff.shape[1]'
  write(pyid,'(A)') "mf.mo_coeff = fch2py(hf_fch, nbf, nif, 'a')"
  write(pyid,'(A)') '# read done'
  write(pyid,'(/,A)') '# check if input MOs are orthonormal'
  if(pbc) then
   write(pyid,'(A)') "S = cell.pbc_intor('int1e_ovlp',hermi=1)"
  else
   write(pyid,'(A)') "S = mol.intor_symmetric('int1e_ovlp')"
  end if
  write(pyid,'(A)') 'check_orthonormal(nbf, nif, mf.mo_coeff, S)'
 end if

 write(pyid,'(/,A)') '#dm = mf.make_rdm1()'
 write(pyid,'(A)')   '#mf.max_cycle = 10'
 write(pyid,'(A,/)') '#mf.kernel(dm0=dm)'
 close(pyid)
end subroutine bas_gms2py

