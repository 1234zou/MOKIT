! written by jxzou at 20210404: generate BDF .BAS, .inp, .scforb/.inporb files
!  from Gaussian .fch(k) file
! updated by jxzou at 20210114: enlarge coeff(nbf,nif) to coeff(nbf,nbf) for
!  linear dependence usage
! updated by jxzou at 20210407: remove '-uhf' and add automatic determination

! TODO: Remove subroutine bas_gms2bdf. Generating GAMESS .inp file is not an
!  efficient way.

program main
 use util_wrapper, only: formchk, fch2inp_wrap
 implicit none
 integer :: i, narg, itype
 ! itype = -5/-4/-3/-2/-1/0/1 for
 !    X-TDDFT/SOC-X-TDA/X-TDA/TD-DFT/DFT/HF/CAS
 character(len=8) :: str8
 character(len=26), parameter :: error_warn = 'ERROR in program fch2bdf: '
 character(len=30) :: dftname
 character(len=240) :: fchname, inpname

 narg = iargc()
 if(narg<1 .or. narg>3) then
  write(6,'(/,1X,A)') error_warn//'wrong command line arguments!'
  write(6,'(A)')   ' Example 1 (R(O)HF/UHF): fch2bdf a.fch'
  write(6,'(A)')   ' Example 2 (KS-DFT)    : fch2bdf a.fch -dft ''B3LYP'''
  write(6,'(A)')   ' Example 3 (KS-DFT)    : fch2bdf a.fch -dft ''GB3LYP'''
  write(6,'(A)')   ' Example 4 (BDF BHHLYP): fch2bdf a.fch -dft ''BHandHLYP'''
  write(6,'(A)')   ' Example 5 (DFT-D3)    : fch2bdf a.fch -dft ''B3LYP D3zero'''
  write(6,'(A)')   ' Example 6 (DFT-D3BJ)  : fch2bdf a.fch -dft ''B3LYP D3BJ'''
  write(6,'(A)')   ' Example 7 (X-TDA)     : fch2bdf a.fch -xtda'
  write(6,'(A)')   ' Example 8 (SOC-X-TDA) : fch2bdf a.fch -socxtda'
  write(6,'(A)')   ' Example 9 (X-TDDFT)   : fch2bdf a.fch -xtddft'
  write(6,'(A,/)') ' Example10 (CAS NO)    : fch2bdf a.fch -no'
  stop
 end if

 itype = 0; str8 = ' '; dftname = ' '; fchname = ' '
 call getarg(1, fchname)
 call require_file_exist(fchname)

 ! if .chk file provided, convert into .fch file automatically
 i = LEN_TRIM(fchname)
 if(fchname(i-3:i) == '.chk') then
  call formchk(fchname)
  fchname = fchname(1:i-3)//'fch'
 end if

 if(narg > 1) then
  call getarg(2, str8)
  select case(TRIM(str8))
  case('-no')
   itype = 1
  case('-dft')
   itype = -1
   if(narg /= 3) then
    write(6,'(/,A)') error_warn//'three command line arguments are required'
    write(6,'(A)') 'when `-dft` is specified.'
    stop
   end if
   call getarg(3, dftname)
   if(INDEX(dftname,'-') > 0) then
    write(6,'(/,A)') error_warn//'''-'' symbol is not allowed in dftname. For'
    write(6,'(A)') 'example: Gaussian B3LYP-D3(BJ) should be specified as `-dft&
                   & ''GB3LYP D3BJ''`'
    stop
   end if
  case('-xtda','-socxtda','-xtddft')
   if(narg /= 2) then
    write(6,'(/,A)') error_warn//'no argument is accepted after'
    write(6,'(A)') '-xtda/-socxtda/-xtddft'
    stop
   end if
   dftname = 'BHHLYP'; itype = -3
   if(TRIM(str8) == '-socxtda') itype = -4
   if(TRIM(str8) == '-xtddft') itype = -5
  case default
   write(6,'(/,A)') error_warn//'the 2nd argument can only be `-no` or `-dft`.'
   stop
  end select
 end if

 call fch2inp_wrap(fchname, .false., 0, 0, .true., .false., .false.)
 call find_specified_suffix(fchname, '.fch', i)
 inpname = fchname(1:i-1)//'.inp'
 call bas_gms2bdf(inpname, itype, dftname)
 call delete_file(inpname)

 call fch2bdf(fchname, itype)
end program main

! read the MOs in .fch(k) file and adjust its p,d,f,g, etc. functions order
!  of Gaussian to that of BDF
subroutine fch2bdf(fchname, itype)
 use fch_content, only: Bohr_const
 implicit none
 integer :: i, j, k, m, length, natom, orbid
 integer :: na, nb, nif, nbf, nbf0, nbf1, mult
 integer :: n3pmark, n5dmark, n7fmark, n9gmark, n11hmark
 integer, intent(in) :: itype
 integer, allocatable :: shell_type(:), shell2atom_map(:)
 ! mark the index where d, f, g, h functions begin
 integer, allocatable :: p_mark(:), d_mark(:), f_mark(:), g_mark(:), h_mark(:)
 integer, allocatable :: nuc(:), ntimes(:)
 real(kind=8), allocatable :: coeff(:,:), occ_num(:), coor(:,:)
 character(len=2), allocatable :: elem(:)
 character(len=29), parameter :: error_warn = 'ERROR in subroutine fch2bdf: '
 character(len=240) :: inpname, orbfile
 character(len=240), intent(in) :: fchname
 logical :: uhf

 ! transform basis set data and MOs
 call check_nobasistransform_in_fch(fchname)
 call check_nosymm_in_fch(fchname)
 call read_ncontr_from_fch(fchname, k)

 allocate(shell_type(2*k), source=0)
 allocate(shell2atom_map(2*k), source=0)
 call read_shltyp_and_shl2atm_from_fch(fchname, k, shell_type, shell2atom_map)

 if( ANY(shell_type>1) ) then ! whether Cartesian/spherical harmonic
  deallocate(shell_type, shell2atom_map)
  write(6,'(/,A)') error_warn//'Cartesian-type basis functions are not supported'
  write(6,'(A)') 'in BDF. You can use spherical harmonic basis functions.'
  stop
 end if

 call read_na_and_nb_from_fch(fchname, na, nb)
 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 nbf0 = nbf   ! save a copy of the original nbf

 ! read MO Coefficients
 ! The array size of occ_num is actually 2*nif or nif, and coeff is (nbf,nif)
 ! or (nbf,2*nif). But BDF requires them to use nbf. If linear dependence
 ! occurs, the extra elements are zero and must be printed into .scforb/.inporb.
 call check_uhf_in_fch(fchname, uhf) ! determine whether UHF
 if(uhf) then
  allocate(occ_num(2*nbf), coeff(nbf,2*nbf))
  call read_eigenvalues_from_fch(fchname,nif,'a',occ_num(1:nif))
  call read_eigenvalues_from_fch(fchname,nif,'b',occ_num(nif+1:2*nif))
  call read_mo_from_fch(fchname, nbf, nif, 'a', coeff(:,1:nif))
  call read_mo_from_fch(fchname, nbf, nif, 'b', coeff(:,nif+1:2*nif))
  nif = 2*nif   ! double the size
 else
  allocate(occ_num(nbf), coeff(nbf,nbf))
  call read_eigenvalues_from_fch(fchname, nif, 'a', occ_num)
  call read_mo_from_fch(fchname, nbf, nif, 'a', coeff(1:nbf,1:nif))
 end if

! first we adjust the basis functions in each MO according to the Shell to atom map
 ! 1) split the 'L' into 'S' and 'P', this is to ensure that D comes after L functions
 call split_L_func(k, shell_type, shell2atom_map, length)

 ! 2) sort the shell_type and shell2atom_map by ascending order
 ! MOs will be adjusted accordingly
 call sort_shell_and_mo(length, shell_type, shell2atom_map, nbf, nif, coeff)
! adjust done

! then we adjust the basis functions in each MO according to the type of basis functions
 k = length  ! update k
 n3pmark = 0
 n5dmark = 0
 n7fmark = 0
 n9gmark = 0
 n11hmark = 0
 allocate(p_mark(k), d_mark(k), f_mark(k), g_mark(k), h_mark(k))
 p_mark = 0
 d_mark = 0
 f_mark = 0
 g_mark = 0
 h_mark = 0
 nbf = 0
 do i = 1, k, 1
  select case(shell_type(i))
  case( 0)   ! S
   nbf = nbf + 1
  case( 1)   ! 3P
   n3pmark = n3pmark + 1
   p_mark(n3pmark) = nbf + 1
   nbf = nbf + 3
  case(-1)   ! SP or L
   n3pmark = n3pmark + 1
   p_mark(n3pmark) = nbf + 2
   nbf = nbf + 4
  case(-2)   ! 5D
   n5dmark = n5dmark + 1
   d_mark(n5dmark) = nbf + 1
   nbf = nbf + 5
  case(-3)   ! 7F
   n7fmark = n7fmark + 1
   f_mark(n7fmark) = nbf + 1
   nbf = nbf + 7
  case(-4)   ! 9G
   n9gmark = n9gmark + 1
   g_mark(n9gmark) = nbf + 1
   nbf = nbf + 9
  case(-5)   ! 11H
   n11hmark = n11hmark + 1
   h_mark(n11hmark) = nbf + 1
   nbf = nbf + 11
  case default
   write(6,'(/,A)') 'ERROR in subroutine fch2bdf: shell_type(i) out of range.'
   write(6,'(2(A,I0))') 'i=', i, ', k=', k
   stop
  end select
 end do ! for i

 ! adjust the order of d, f, etc. functions
 do i = 1, n3pmark, 1
  call fch2bdf_permute_3p(nif,coeff(p_mark(i):p_mark(i)+2,:))
 end do
 do i = 1, n5dmark, 1
  call fch2bdf_permute_5d(nif,coeff(d_mark(i):d_mark(i)+4,:))
 end do
 do i = 1, n7fmark, 1
  call fch2bdf_permute_7f(nif,coeff(f_mark(i):f_mark(i)+6,:))
 end do
 do i = 1, n9gmark, 1
  call fch2bdf_permute_9g(nif,coeff(g_mark(i):g_mark(i)+8,:))
 end do
 do i = 1, n11hmark, 1
  call fch2bdf_permute_11h(nif,coeff(h_mark(i):h_mark(i)+10,:))
 end do
! adjustment finished

 deallocate(p_mark, d_mark, f_mark, g_mark, h_mark)

! move the 2nd, 3rd, ... Zeta basis functions forward
 i = 0
 nbf = 0
 do while(i < k)
  i = i + 1
  j = shell2atom_map(i)
  m = shell_type(i)
  nbf1 = nbf
  select case(m)
  case( 0)   ! S
   nbf = nbf + 1
  case( 1)   ! 3P
   nbf = nbf + 3
  case(-1)   ! SP or L
   nbf = nbf + 4
  case(-2)   ! 5D
   nbf = nbf + 5
  case(-3)   ! 7F
   nbf = nbf + 7
  case(-4)   ! 9G
   nbf = nbf + 9
  case(-5)   ! 11H
   nbf = nbf + 11
  end select
  if(m == 0) cycle

  length = 1
  do while(i < k)
   i = i + 1
   if(shell_type(i) /= m) exit
   if(shell2atom_map(i) /= j) exit
   length = length + 1
  end do
  if(i < k) i = i - 1
  if(length > 1) then
   call zeta_mv_forwd(nbf1, m, length, nbf0, nif, coeff)
   nbf = nbf1 + length*(nbf-nbf1)
  end if
 end do
 deallocate(shell_type, shell2atom_map)
! move done

 ! print MOs into BDF .scforb/.inporb
 i = INDEX(fchname, '.fch', back=.true.)
 inpname = fchname(1:i-1)//'_bdf.inp'
 if(itype == 1) then
  orbfile = fchname(1:i-1)//'_bdf.inporb'
 else
  orbfile = fchname(1:i-1)//'_bdf.scforb'
 end if

 open(newunit=orbid,file=TRIM(orbfile),status='replace')
 write(orbid,'(A)') 'TITLE - transformed MOs'

 call read_charge_and_mult_from_fch(fchname, i, mult)

 if((mult==1 .and. (.not.uhf)) .or. itype==1) then
  write(orbid,'(A)') '$MOCOEF  1'  ! RHF/CAS
 else
  write(orbid,'(A)') '$MOCOEF  2'  ! ROHF, UHF
 end if

 if(uhf) nif = nif/2
 write(orbid,'(A,I9,A)') 'SYM=  1 NORB=', nbf0, ' ALPHA'
 do i = 1, nbf0, 1
  write(orbid,'(5(2X,E23.16))') (coeff(j,i),j=1,nbf0)
 end do ! for i

 if(uhf) then    ! UHF
  write(orbid,'(A,I9,A)') 'SYM=  1 NORB=', nbf0, ' BETA'
  do i = 1, nbf0, 1
   write(orbid,'(5(2X,E23.16))') (coeff(j,i+nbf0),j=1,nbf0)
  end do ! for i
 else if(mult/=1 .and. itype/=1) then ! ROHF
  write(orbid,'(A,I9,A)') 'SYM=  1 NORB=', nbf0, ' BETA'
  do i = 1, nbf0, 1
   write(orbid,'(5(2X,E23.16))') (coeff(j,i),j=1,nbf0)
  end do ! for i
 end if
 deallocate(coeff)

 write(orbid,'(A)') 'ORBITAL ENERGY'
 write(orbid,'(5(2X,E23.16))') (occ_num(i),i=1,nbf0)

 if(uhf) then
  write(orbid,'(5(2X,E23.16))') (occ_num(i),i=nbf0+1,2*nbf0)
 else if(mult/=1 .and. itype/=1) then
  write(orbid,'(5(2X,E23.16))') (occ_num(i),i=1,nbf0)
 end if

 write(orbid,'(A)') 'OCCUPATION'
 if(itype == 1) then
  write(orbid,'(10(1X,E11.5))') (occ_num(i),i=1,nbf0)
  write(orbid,'(A)') '$END'
  close(orbid)
  return
 end if

 occ_num = 0d0
 occ_num(1:na) = 1d0
 write(orbid,'(10(1X,E11.5))') (occ_num(i),i=1,nbf0)

 if(mult/=1 .or. uhf) then
  if(nb < na) occ_num(nb+1:na) = 0d0
  write(orbid,'(A)') 'OCCUPATION'
  write(orbid,'(10(1X,E11.5))') (occ_num(i),i=1,nbf0)
 end if
 deallocate(occ_num)
 write(orbid,'(A)') '$END'

 call read_natom_from_fch(fchname, natom)
 write(orbid,'(A,1X,I8)') '$COORD', natom

 allocate(elem(natom), coor(3,natom), nuc(natom), ntimes(natom))
 call read_elem_and_coor_from_fch(fchname, natom, elem, nuc, coor, i, j)
 deallocate(nuc)
 call calc_ntimes(natom, elem, ntimes)
 coor = coor/Bohr_const

 do i = 1, natom, 1
  write(orbid,'(A,I0,3(1X,F19.12))') TRIM(elem(i)),ntimes(i),coor(1:3,i)
 end do ! for i
 deallocate(elem, ntimes, coor)

 write(orbid,'(A)') '$END'

 if(mult /= 1) then
  write(orbid,'(A,2(1X,I8),A)') '$NBF', nbf0, nbf0, '   2'
 else
  write(orbid,'(A,2(1X,I8),A)') '$NBF', nbf0, nbf0, '   1'
 end if

 write(orbid,'(A,2X,I7,1X,I7)') 'NOCC', na, nb
 write(orbid,'(A)') '$END'
 close(orbid)
 ! done print MOs
end subroutine fch2bdf

! print primitive gaussians
subroutine prt_prim_gau_bdf(iatom, fid)
 use pg, only: prim_gau, nuc, highest, all_ecp, elem, ntimes, ecp_exist
 implicit none
 integer :: i, j, k, m, n, nline, ncol
 integer, intent(in) :: iatom, fid
 integer, allocatable :: list(:)
 character(len=1), parameter :: am(0:5) = ['s','p','d','f','g','h']

 call get_highest_am()
 write(fid,'(A)') '****'
 write(fid,'(A,I0,3X,I0,3X,I1)') TRIM(elem(iatom)),ntimes(iatom),nuc(iatom),highest

 do i = 1, 7, 1
  if(.not. allocated(prim_gau(i)%coeff)) cycle
  write(fid,'(A)',advance='no') prim_gau(i)%stype
  nline = prim_gau(i)%nline
  ncol = prim_gau(i)%ncol
  write(fid,'(2(1X,I4))') nline, ncol-1
  do j = 1, nline, 1
   write(fid,'(3X,E16.9)') prim_gau(i)%coeff(j,1)
  end do ! for j
  do j = 1, nline, 1
   write(fid,'(10(E16.9,3X))') (prim_gau(i)%coeff(j,k), k=2,ncol)
  end do ! for j
 end do ! for i

 if(.not. ecp_exist) return

 if(all_ecp(iatom)%ecp) then
  write(fid,'(A)') 'ECP'
  m = all_ecp(iatom)%highest
  write(fid,'(A,2(I0,1X),I0)') TRIM(elem(iatom)),ntimes(iatom),all_ecp(iatom)%core_e,m
  allocate(list(m+1))
  list(1) = m
  forall(i=2:m+1) list(i) = i-2
  do i = 1, m+1, 1
   n = all_ecp(iatom)%potential(i)%n
   if(i == 1) then
    write(fid,'(A,I0)') am(list(i))//' potential ', n
   else ! i > 1
    write(fid,'(A,I0)') am(list(i))//' potential ', n
   end if

   do j = 1, n, 1
    write(fid,'(I0,2(1X,F16.8))') all_ecp(iatom)%potential(i)%col2(j), all_ecp(iatom)%potential(i)%col3(j), &
                                  all_ecp(iatom)%potential(i)%col1(j)
   end do ! for j
  end do ! for i
 end if
end subroutine prt_prim_gau_bdf

! Transform the basis sets in GAMESS format to those in BDF format
subroutine bas_gms2bdf(fort7, itype, dftname)
 use pg, only: natom, nuc, ntimes, coor, elem, all_ecp, ecp_exist
 implicit none
 integer :: i, k, nline, rc, rel, nbf, nif, fid1, fid2
 integer :: charge, mult, isph
 integer, intent(in) :: itype
 character(len=7) :: str
 character(len=240), intent(in) :: fort7
 character(len=240) :: buf, input, basfile
 ! input is the BDF input file
 ! if you do not like the suffix .bdf, you can change it into .inp
 character(len=1) :: stype
 character(len=30) :: dftname1
 character(len=30), intent(in) :: dftname
 character(len=33), parameter :: error_str = 'ERROR in subroutine bas_gms2bdf: '
 logical :: X2C, uhf, ghf, ks_dft, d3zero, d3bj
 logical, allocatable :: ghost(:)

 ! detect and parse dftname
 ks_dft = .false.; d3zero = .false.; d3bj = .false.
 dftname1 = ADJUSTL(dftname)
 call upper(dftname1)
 k = LEN_TRIM(dftname1)

 if(k > 0) then
  ks_dft = .true.
  if(dftname1(1:5) == 'B3LYP') then
   write(6,'(/,A)') 'Warning from subroutine fch2bdf: if you want to use Gaussi&
                    &an B3LYP functional'
   write(6,'(A)') 'in BDF, you need to specify `GB3LYP` (which uses VWN3). An e&
                  &xample is:'
   write(6,'(A)') ' fch2bdf h2o_rks.fch -dft ''GB3LYP D3BJ'''
   write(6,'(/,A)') 'If you simply want to use BDF default B3LYP (which uses VW&
                    &N5), you can ignore'
   write(6,'(A)') 'this warning.'
  end if
  i = INDEX(dftname1(1:k), ' ', back=.true.)
  if(i > 0) then
   select case(TRIM(dftname1(i+1:k)))
   case('D3')
    write(6,'(/,A)') error_str//'D3 is not accepted in this utility. Please'
    write(6,'(A)') 'specify D3zero or D3BJ.'
    stop
   case('D3ZERO')
    d3zero = .true.
   case('D3BJ')
    d3bj = .true.
   case default
    write(6,'(/,A)') error_str//'dftname cannot be parsed correctly.'
    write(6,'(A)') 'dftname='//TRIM(dftname)
    stop
   end select
   dftname1(i+1:k) = ' '
  end if
  if(TRIM(dftname1) == 'BHANDHLYP') dftname1 = 'BHHLYP'
 end if

 buf = ' '; input = ' '; basfile = ' '
 call read_nbf_and_nif_from_gms_inp(fort7, nbf, nif)

 k = INDEX(fort7, '.', back=.true.)
 input = fort7(1:k-1)//'_bdf.inp'
 call read_natom_from_gms_inp(fort7, natom)
 allocate(elem(natom), nuc(natom), coor(3,natom), ntimes(natom), ghost(natom))
 call read_elem_nuc_coor_from_gms_inp(fort7, natom, elem, nuc, coor, ghost)
 deallocate(ghost)
 ! ram cannot be deallocated here since subroutine prt_prim_gau_bdf will use it

 call calc_ntimes(natom, elem, ntimes)
 call read_charge_mult_isph_from_gms_inp(fort7, charge, mult, isph, uhf, ghf, &
                                         ecp_exist)
 call read_all_ecp_from_gms_inp(fort7)

 if(itype<-2 .and. mult<2) then
  write(6,'(/,A)') error_str//'X-TD must be based on ROKS. But got'
  write(6,'(A,I0)') 'mult=', mult
  call delete_file(TRIM(fort7))
  stop
 end if

 if(ecp_exist) then
  basfile = 'ECP.'//fort7(1:k-1)//'.BAS'
 else
  basfile = fort7(1:k-1)//'.BAS'
 end if
 ! read ECP/PP done

 call upper(basfile)
 open(newunit=fid2,file=TRIM(input),status='replace')
 write(fid2,'(A)') '$COMPASS'
 write(fid2,'(A)') 'Title'
 write(fid2,'(A)') 'generated by fch2bdf of MOKIT'
 write(fid2,'(A)') 'Geometry'

 call check_X2C_in_gms_inp(fort7, X2C)
 call check_DKH_in_gms_inp(fort7, rel)

 do i = 1, natom, 1
  str = ' '
  write(str,'(A,I0)') TRIM(elem(i)), ntimes(i)
  write(fid2,'(A7,1X,3F18.8)') str, coor(1:3,i)
 end do ! for i
 deallocate(coor)

 write(fid2,'(A)') 'End Geometry'
 write(fid2,'(A,/,A)') 'Basis', TRIM(basfile)
 write(fid2,'(A)') 'Nosymm'
 if(itype == 1) write(fid2,'(A)') 'Saorb'
 write(fid2,'(A)') '$END'

 write(fid2,'(/,A)') '$XUANYUAN'
 if(rel>-1 .or. X2C) then
  write(fid2,'(A)') '# for very old BDF, you may need to replace the following &
                    &two lines by ''Scalar'''
  write(fid2,'(A,/,1X,I0)') 'heff', 3
  if(rel>-1 .and. (.not.X2C)) then
   write(6,'(/,A)') 'Warning in subroutine bas_gms2bdf: BDF program does not su&
                    &pport DKH Hamiltonian.'
   write(6,'(A)') "Spin-free X2C keyword 'Scalar' are written in $XUANYUAN inst&
                  &ead."
  end if
 end if
 if(itype == -4) write(fid2,'(A,/,1X,I0)') 'hsoc', 2
 write(fid2,'(A)') '$END'

 write(fid2,'(/,A)') '$SCF'

 if(ks_dft) then ! KS-DFT
  if(uhf) then
   write(fid2,'(A)') 'UKS'
  else
   if(mult == 1) then
    write(fid2,'(A)') 'RKS'
   else
    write(fid2,'(A)') 'ROKS'
   end if
  end if
  write(fid2,'(A)') 'DFT'
  write(fid2,'(1X,A)') TRIM(dftname1)
  if(d3zero) write(fid2,'(A)') 'D3zero'
  if(d3bj) write(fid2,'(A)') 'D3' ! be careful
  write(fid2,'(A,/,1X,A)') 'Grid', 'ultra fine'
  write(fid2,'(A,/,1X,A)') 'GridTol', '1e-8'
 else            ! HF
  if(uhf) then
   write(fid2,'(A)') 'UHF'
  else
   if(mult == 1) then
    write(fid2,'(A)') 'RHF'
   else
    write(fid2,'(A)') 'ROHF'
   end if
  end if
 end if

 write(fid2,'(A,/,1X,I0)') 'Charge', charge
 write(fid2,'(A,/,1X,I0)') 'Spin', mult
 write(fid2,'(A,/,1X,A)') 'Guess','read'
 ! By default, BDF uses 1d-5 for checking basis set linear dependency, which
 ! is not equal to that of Gaussian/GAMESS. So let's always specify 1d-6.
 write(fid2,'(2(A,/),1X,A)') 'CheckLin','TolLin','1d-6'
 ! Other options which might be used in the future:
 ! Vshift=x: equivalent to scf(vshift=1000*x) of Gaussian
 ! NoDeltaP: turn off the incremental Fock technique, usually useless
 ! NoGridSwitch: equivalent to `NoVarAcc` of Gaussian, which seems useless when
 !               `ultra fine` is specified
 write(fid2,'(A)') '$END'

 if(itype == -5) then ! X-TDDFT
  write(fid2,'(/,A)') '$TDDFT'
  write(fid2,'(A,/,1X,I0)') 'iroot', 4
  write(fid2,'(A,/,1X,A)') 'Grid', 'ultra fine'
  write(fid2,'(A)') '$END'
 else if(itype==-3 .or. itype==-4) then ! (SOC-)X-TDA
  write(fid2,'(/,A)') '$TDDFT'
  write(fid2,'(A,/,1X,I0)') 'itda', 1
  write(fid2,'(A,/,1X,I0)') 'iroot', 4
  if(itype == -4) then
   write(fid2,'(A,/,1X,I0)') 'istore', 1
   write(fid2,'(A,/,1X,I0)') 'itrans', 1
  end if
  write(fid2,'(A,/,1X,A)') 'Grid', 'ultra fine'
  write(fid2,'(A)') '$END'
  write(fid2,'(/,A)') '$TDDFT'
  write(fid2,'(A,/,1X,I0)') 'itda', 1
  write(fid2,'(A,/,1X,I0)') 'isf', 1
  ! Important: ALDA0 noncollinear XC kernel, affect excited states much
  write(fid2,'(A,/,1X,I0)') 'ialda', 2
  write(fid2,'(A,/,1X,I0)') 'iroot', 4
  if(itype == -4) then
   write(fid2,'(A,/,1X,I0)') 'istore', 2
   write(fid2,'(A,/,1X,I0)') 'itrans', 1
  end if
  write(fid2,'(A,/,1X,A)') 'Grid', 'ultra fine'
  write(fid2,'(A)') '$END'
  ! so far, SOC calculation is supported using X-TDA, but not X-TDDFT
  if(itype == -4) then
   write(fid2,'(/,A)') '$TDDFT'
   write(fid2,'(A,/,1X,I0)') 'isoc', 2
   write(fid2,'(A,/,1X,I0)') 'nfiles', 2
   write(fid2,'(A,/,1X,I0)') 'ifgs', 1
   write(fid2,'(A,/,1X,I0)') 'imatsoc', -1
   write(fid2,'(A,/,1X,I0)') 'imatrso', -2
   write(fid2,'(A)') '$END'
  end if
 end if
 close(fid2)

 ! open the .inp/.dat file and find the $DATA section
 call goto_data_section_in_gms_inp(fort7, fid1)
 read(fid1,'(A)') buf
 read(fid1,'(A)') buf

 ! initialization: clear all primitive gaussians
 call clear_prim_gau()
 open(newunit=fid2,file=TRIM(basfile),status='replace')

 ! read the element, relative atomic mass (ram) and coordinates
 do i = 1, natom, 1
  read(fid1,'(A)',iostat=rc) buf
  ! 'buf' contains the element, ram and coordinates
  if(rc /= 0) exit

  ! deal with primitive gaussians
  do while(.true.)
   read(fid1,'(A)') buf
   if(LEN_TRIM(buf) == 0) exit

   read(buf,*) stype, nline
   call read_prim_gau(stype, nline, fid1)
  end do ! for while

  ! print basis sets and ECP/PP (if any) of this atom in BDF format
  call prt_prim_gau_bdf(i, fid2)

  ! clear all primitive gaussians for next cycle
  call clear_prim_gau()
 end do ! for i

 deallocate(nuc, ntimes, elem, all_ecp)
 close(fid1)

 if(rc /= 0) then
  write(6,'(/,A)') error_str//'it seems the "$DATA" has no corresponding "$END".'
  write(6,'(A)') 'Incomplete file '//TRIM(fort7)
  close(fid2,status='delete')
  stop
 end if

 close(fid2)
end subroutine bas_gms2bdf

! move the 2nd, 3rd, ... Zeta basis functions forward
subroutine zeta_mv_forwd(i0, shell_type, length, nbf, nif, coeff2)
 implicit none
 integer i, j, k
 integer, intent(in) :: i0, shell_type, length, nbf, nif
 integer, parameter :: num0(-5:5) = [11, 9, 7, 5, 0, 0, 3, 6, 10, 15, 21]
 !                                   11H 9G 7F 5D L  S 3P 6D 10F 15G 21H
 real(kind=8), intent(inout) :: coeff2(nbf,nif)
 real(kind=8), allocatable :: coeff(:,:)

 if(length == 1) return

 if(shell_type==0 .or. shell_type==-1) then
  write(6,'(/,A)') 'ERROR in subroutine zeta_mv_forwd: this element of shell_ty&
                   &pe is 0 or -1.'
  write(6,'(2(A,I0))') 'shell_type=', shell_type, ', length=', length
  stop
 end if

 coeff = coeff2
 k = num0(shell_type)
 do i = 1, k, 1
  do j = 1, length, 1
   coeff(i0+j+(i-1)*length,:) = coeff2(i0+i+(j-1)*k,:)
  end do ! for j
 end do ! for i

 coeff2 = coeff
 deallocate(coeff)
end subroutine zeta_mv_forwd

subroutine fch2bdf_permute_3p(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(3) = [2,3,1]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(3,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical p functions in Gaussian
! To: the order of spherical p functions in BDF
! 1    2    3
! Px,  Py,  Pz
! Py,  Pz,  Px

 allocate(coeff2(3,nif), source=0d0)
 forall(i = 1:3) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
end subroutine fch2bdf_permute_3p

subroutine fch2bdf_permute_5d(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(5) = [5, 3, 1, 2, 4]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(5,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical d functions in Gaussian
! To: the order of spherical d functions in BDF
! 1    2    3    4    5
! d0 , d+1, d-1, d+2, d-2
! d-2, d-1, d0 , d+1, d+2

 allocate(coeff2(5,nif), source=0d0)
 forall(i = 1:5) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
end subroutine fch2bdf_permute_5d

subroutine fch2bdf_permute_7f(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(7) = [7, 5, 3, 1, 2, 4, 6]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(7,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical f functions in Gaussian
! To: the order of spherical f functions in BDF
! 1    2    3    4    5    6    7
! f0 , f+1, f-1, f+2, f-2, f+3, f-3
! f-3, f-2, f-1, f0 , f+1, f+2, f+3

 allocate(coeff2(7,nif), source=0d0)
 forall(i = 1:7) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
end subroutine fch2bdf_permute_7f

subroutine fch2bdf_permute_9g(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(9) = [9, 7, 5, 3, 1, 2, 4, 6, 8]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(9,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical g functions in Gaussian
! To: the order of spherical g functions in BDF
! 1    2    3    4    5    6    7    8    9
! g0 , g+1, g-1, g+2, g-2, g+3, g-3, g+4, g-4
! g-4, g-3, g-2, g-1, g0 , g+1, g+2, g+3, g+4

 allocate(coeff2(9,nif), source=0d0)
 forall(i = 1:9) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
end subroutine fch2bdf_permute_9g

subroutine fch2bdf_permute_11h(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(11) = [11, 9, 7, 5, 3, 1, 2, 4, 6, 8, 10]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(11,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical h functions in Gaussian
! To: the order of spherical h functions in BDF
! 1    2    3    4    5    6    7    8    9    10   11
! h0 , h+1, h-1, h+2, h-2, h+3, h-3, h+4, h-4, h+5, h-5
! h-5, h-4, h-3, h-2, h-1, h0 , h+1, h+2, h+3, h+4, h+5

 allocate(coeff2(11,nif), source=0d0)
 forall(i = 1:11) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
end subroutine fch2bdf_permute_11h

