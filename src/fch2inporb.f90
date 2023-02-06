! written by jxzou at 20180817: adjust the orders of d,f,g, etc functions in
!  .fch(k) file, to the orders in Molcas and OpenMolcas

! updated by jxzou at 20180826: support Cartesian H functions
! updated by jxzou at 20190226: change data to parameter when declare constant arrays
! updated by jxzou at 20190624: able to transform UHF orbitals
! updated by jxzou at 20200504: combined with rwwfn.f90
! updated by jxzou at 20200809: combined with util_wrapper.f90
! updated by jxzou at 20210407: remove '-uhf', add automatic determination

program main
 use util_wrapper, only: formchk, fch2inp_wrap
 implicit none
 integer :: i, k, system
 character(len=3) :: str
 character(len=240) :: fname, inpname
 logical :: sph, prt_no

 i = iargc()
 if(i<1 .or. i>2) then
  write(6,'(/,A)') ' ERROR in subroutine fch2inporb: wrong command line arguments!'
  write(6,'(A)')   ' Example 1 (R(O)HF, UHF, CAS): fch2inporb a.fch'
  write(6,'(A,/)') ' Example 2 (for CAS NO)      : fch2inporb a.fch -no'
  stop
 end if

 fname = ' '
 call getarg(1, fname)
 call require_file_exist(fname)

 ! if .chk file provided, convert into .fch file automatically
 k = LEN_TRIM(fname)
 if(fname(k-3:k) == '.chk') then
  call formchk(fname)
  fname = fname(1:k-3)//'fch'
 end if

 call check_nobasistransform_in_fch(fname)
 call check_nosymm_in_fch(fname)

 if(i == 2) then
  str = ' '
  call getarg(2, str)
  if(str /= '-no') then
   write(6,'(/,A)') " ERROR in subroutine fch2inporb: the 2nd argument is&
                    & wrong! Only '-no' is accepted."
   stop
  else ! str = '-no'
   prt_no = .true.
  end if
 else ! i = 1
  prt_no = .false.
 end if

 ! ->prt_no, <-sph
 call fch2inporb(fname, prt_no, sph)
 call fch2inp_wrap(fname, .false., 0, 0)

 k = index(fname,'.fch', back=.true.)
 inpname = fname(1:k-1)//'.inp'

 if(sph) then
  i = system('bas_gms2molcas '//TRIM(inpname)//' -sph')
 else
  i = system('bas_gms2molcas '//TRIM(inpname))
 end if

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine fch2inporb: call utility bas_gms2molcas failed.'
  write(6,'(A)')   'Three possible reasons:'
  write(6,'(A)')   '(1) You forget to compile the utility bas_gms2molcas.'
  write(6,'(A)')   '(2) The file '//TRIM(fname)//' may be incomplete.'
  write(6,'(A,/)') '(3) This is a bug of the utility bas_gms2molcas.'
  stop
 end if

 call delete_file(inpname)
end program main

! nbf: the number of basis functions
! nif: the number of independent functions, i.e., the number of MOs

! read the MOs in .fch(k) file and adjust its d,f,g, etc. functions order
!  of Gaussian to that of Molcas
subroutine fch2inporb(fchname, prt_no, sph)
 use fch_content, only: check_uhf_in_fch
 implicit none
 integer :: i, j, k, m, length, orbid
 integer :: nalpha, nbeta, nbf, nif
 integer :: nbf0, nbf1
 integer :: n6dmark,n10fmark,n15gmark,n21hmark
 integer :: n5dmark,n7fmark, n9gmark, n11hmark
 integer(kind=4) :: getpid, hostnm
 integer, allocatable :: shell_type(:), shell2atom_map(:)
 ! mark the index where d, f, g, h functions begin
 integer, allocatable :: d_mark(:), f_mark(:), g_mark(:), h_mark(:)
 character(len=8) :: hostname
 character(len=24) :: data_string
 character(len=240), intent(in) :: fchname
 character(len=240) :: orbfile   ! INPORB file of (Open)Molcas
 real(kind=8), allocatable :: coeff(:,:), occ_num(:)
 logical :: uhf
 logical, intent(in) :: prt_no
 logical, intent(out) :: sph

 sph = .true. ! default value
 orbfile = ' '

 call read_na_and_nb_from_fch(fchname, nalpha, nbeta)
 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 nbf0 = nbf

 call check_uhf_in_fch(fchname, uhf) ! determine whether UHF

 ! read MO Coefficients
 if(uhf) then
  allocate(coeff(nbf,2*nif))
  call read_mo_from_fch(fchname, nbf, nif, 'a', coeff(:,1:nif))
  call read_mo_from_fch(fchname, nbf, nif, 'b', coeff(:,nif+1:2*nif))
  nif = 2*nif   ! double the size
 else   ! not UHF
  if(prt_no) then
   allocate(occ_num(nif), coeff(nbf,nif))
   call read_eigenvalues_from_fch(fchname, nif, 'a', occ_num)
   call read_mo_from_fch(fchname, nbf, nif, 'a', coeff)
  else
   allocate(coeff(nbf,nif))
   call read_mo_from_fch(fchname, nbf, nif, 'a', coeff)
  end if
 end if

 call read_ncontr_from_fch(fchname, k)
 allocate(shell_type(2*k), source=0)
 allocate(shell2atom_map(2*k), source=0)
 call read_shltyp_and_shl2atm_from_fch(fchname, k, shell_type, shell2atom_map)

 ! check if any spherical functions
 if(ANY(shell_type<-1) .and. ANY(shell_type>1)) then
  write(6,'(A)') 'ERROR in subroutine fch2inporb: mixed spherical harmonic/&
                 &Cartesian functions detected.'
  write(6,'(A)') 'You probably used a basis set like 6-31G(d) in Gaussian. Its&
                 & default setting is (6D,7F).'
  write(6,'(A)') "You need to add '5D 7F' or '6D 10F' keywords in Gaussian."
  stop
 else if( ANY(shell_type<-1) ) then
  sph = .true.
 else
  sph = .false.
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
 n6dmark = 0
 n10fmark = 0
 n15gmark = 0
 n21hmark = 0
 n5dmark = 0
 n7fmark = 0
 n9gmark = 0
 n11hmark = 0
 allocate(d_mark(k), f_mark(k), g_mark(k), h_mark(k))
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
   nbf = nbf + 3
  case(-1)   ! SP or L
   nbf = nbf + 4
  case(-2)   ! 5D
   n5dmark = n5dmark + 1
   d_mark(n5dmark) = nbf + 1
   nbf = nbf + 5
  case( 2)   ! 6D
   n6dmark = n6dmark + 1
   d_mark(n6dmark) = nbf + 1
   nbf = nbf + 6
  case(-3)   ! 7F
   n7fmark = n7fmark + 1
   f_mark(n7fmark) = nbf + 1
   nbf = nbf + 7
  case( 3)   ! 10F
   n10fmark = n10fmark + 1
   f_mark(n10fmark) = nbf + 1
   nbf = nbf + 10
  case(-4)   ! 9G
   n9gmark = n9gmark + 1
   g_mark(n9gmark) = nbf + 1
   nbf = nbf + 9
  case( 4)   ! 15G
   n15gmark = n15gmark + 1
   g_mark(n15gmark) = nbf + 1
   nbf = nbf + 15
  case(-5)   ! 11H
   n11hmark = n11hmark + 1
   h_mark(n11hmark) = nbf + 1
   nbf = nbf + 11
  case( 5)   ! 21H
   n21hmark = n21hmark + 1
   h_mark(n21hmark) = nbf + 1
   nbf = nbf + 21
  case default
   write(6,'(A)') 'ERROR in subroutine fch2inporb: shell_type(i) out of range.'
   write(6,'(2(A,I0))') 'k=', k, ', i=', i
   stop
  end select
 end do ! for i

 ! adjust the order of d, f, etc. functions
 call fch2inporb_permute_sph(n5dmark, n7fmark, n9gmark, n11hmark, k, d_mark, &
                             f_mark, g_mark, h_mark, nbf, nif, coeff)
 call fch2inporb_permute_cart(n6dmark, n10fmark, n15gmark, n21hmark, k, d_mark, &
                              f_mark, g_mark, h_mark, nbf, nif, coeff)
 !do i = 1, n6dmark, 1
 ! call fch2inporb_permute_6d(nif,coeff(d_mark(i):d_mark(i)+5,:))
 !end do
 !do i = 1, n10fmark, 1
 ! call fch2inporb_permute_10f(nif,coeff(f_mark(i):f_mark(i)+9,:))
 !end do
 !do i = 1, n15gmark, 1
 ! call fch2inporb_permute_15g(nif,coeff(g_mark(i):g_mark(i)+14,:))
 !end do
 !do i = 1, n21hmark, 1
 ! call fch2inporb_permute_21h(nif,coeff(h_mark(i):h_mark(i)+20,:))
 !end do
! adjustment finished
 deallocate(d_mark, f_mark, g_mark, h_mark)

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
  case( 2)   ! 6D
   nbf = nbf + 6
  case(-3)   ! 7F
   nbf = nbf + 7
  case( 3)   ! 10F
   nbf = nbf + 10
  case(-4)   ! 9G
   nbf = nbf + 9
  case( 4)   ! 15G
   nbf = nbf + 15
  case(-5)   ! 11H
   nbf = nbf + 11
  case( 5)   ! 21H
   nbf = nbf + 21
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

! print MOs into INPORB
 i = SCAN(fchname, '.', BACK=.TRUE.)
 orbfile = fchname(1:i-1)//'.INPORB'
 open(newunit=orbid,file=TRIM(orbfile),status='replace')

 write(orbid,'(A,/,A)') '#INPORB 2.2', '#INFO'

 if(uhf) then
  write(orbid,'(A)') '* UHF orbitals'
  write(orbid,'(A)') '       1       1       4'
  nif = nif/2   ! change to original size
 else
  if(prt_no) then
   write(orbid,'(A)') '* natural orbitals'
   write(orbid,'(A)') '       0       1       0'
  else
   write(orbid,'(A)') '* SCF orbitals'
   write(orbid,'(A)') '       0       1       2'
  end if
 end if

 write(orbid,'(2X,I6,/,2X,I6)') nbf0, nif

 hostname = ' '
 data_string = ' '
 i = getpid()
 j = hostnm(hostname)
 call fdate(data_string)
 write(orbid,'(A,I6,A)') '*BC:HOST '//TRIM(hostname)//' PID', i, ' DATE '//&
                         TRIM(data_string)//' Generated by MOKIT'

 write(orbid,'(A)') '#ORB'
 do i = 1, nif, 1
  write(orbid,'(A14,I5)') '* ORBITAL    1', i
  write(orbid,'(5(1X,ES21.14))') (coeff(j,i),j=1,nbf0)
 end do

 if(uhf) then
  write(orbid,'(A)') '#UORB'
  do i = 1, nif, 1
   write(orbid,'(A14,I5)') '* ORBITAL    1', i
   write(orbid,'(5(1X,ES21.14))') (coeff(j,i+nif),j=1,nbf0)
  end do
 end if
 deallocate(coeff)

 if(uhf) then
  allocate(occ_num(2*nif), source=0d0)
  occ_num(1:nalpha) = 1d0
  occ_num(nif+1:nif+nbeta) = 1d0
  write(orbid,'(A,/,A)') '#OCC', '* OCCUPATION NUMBERS'
  write(orbid,'(5(1X,ES21.14))') (occ_num(i),i=1,nif)
  write(orbid,'(A,/,A)') '#UOCC', '* Beta OCCUPATION NUMBERS'
  write(orbid,'(5(1X,ES21.14))') (occ_num(i),i=nif+1,2*nif)
 else
  if(prt_no) then
   write(orbid,'(A,/,A)') '#OCC', '* OCCUPATION NUMBERS'
   write(orbid,'(5(1X,ES21.14))') (occ_num(i),i=1,nif)
  else
   allocate(occ_num(nif), source=0d0)
   if(nbeta < nalpha) then
    occ_num(1:nbeta) = 2d0
    occ_num(nbeta+1:nalpha) = 1d0
   else if(nbeta == nalpha) then
    occ_num(1:nalpha) = 2d0
   end if
   write(orbid,'(A,/,A)') '#OCC', '* OCCUPATION NUMBERS'
   write(orbid,'(5(1X,ES21.14))') (occ_num(i),i=1,nif)
  end if
 end if

 write(orbid,'(A,/,A)') '#OCHR', '* OCCUPATION NUMBERS (HUMAN-READABLE)'
 write(orbid,'(10(1X,F7.4))') (occ_num(i),i=1,nif)
 if(uhf) then
  write(orbid,'(A,/,A)') '#UOCHR', '* Beta OCCUPATION NUMBERS (HUMAN-READABLE)'
  write(orbid,'(10(1X,F7.4))') (occ_num(i),i=nif+1,2*nif)
 end if

 deallocate(occ_num)
 close(orbid)
! print done
end subroutine fch2inporb

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
  write(6,'(A)') 'ERROR in subroutine zeta_mv_forwd: this element of&
                 & shell_type is 0 or -1. Impossible.'
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

