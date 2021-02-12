! written by jxzou at 20180817: adjust the orders of d,f,g, etc functions in
!  .fch(k) file, to the orders in Molcas and OpenMolcas

! updated by jxzou at 20180826: support Cartesian H functions
! updated by jxzou at 20190226: change data to parameter when declare constant arrays
! updated by jxzou at 20190624: able to transform UHF orbitals
! updated by jxzou at 20200504: combined with rwwfn.f90
! updated by jxzou at 20200809: combined with util_wrapper.f90

module root_parameter
 implicit none
 real(kind=8), parameter :: root3   = DSQRT(3.0d0)     ! SQRT(3)
 real(kind=8), parameter :: root9   = 3.0d0            ! SQRT(9)
 real(kind=8), parameter :: root15  = DSQRT(15.0d0)    ! SQRT(15)
 real(kind=8), parameter :: root45  = DSQRT(45.0d0)    ! SQRT(45)
 real(kind=8), parameter :: root105 = DSQRT(105.0d0)   ! SQRT(105)
 real(kind=8), parameter :: root945 = DSQRT(945.0d0)   ! SQRT(945)
end module

program main
 use util_wrapper, only: fch2inp_wrap
 implicit none
 integer :: i, k, system
 integer, parameter :: iout = 6
 character(len=4) :: ab
 character(len=240) :: fname, inpname
 logical :: sph, uhf

 i = iargc()
 if(.not. (i==1 .or. i==2)) then
  write(iout,'(/,1X,A)') 'ERROR in subroutine fch2inporb: wrong command line arguments!'
  write(iout,'(1X,A)') 'Example 1 (for RHF, CAS): fch2inporb a.fch'
  write(iout,'(1X,A)') 'Example 2 (for CAS NO)  : fch2inporb a.fch -no'
  write(iout,'(1X,A,/)') 'Example 3 (for UHF)     : fch2inporb a.fch -uhf'
  stop
 end if

 ab = ' '; fname = ' '
 call getarg(1,fname)

 if(i == 2) then
  call getarg(2, ab)
  ab = ADJUSTL(ab)
  if(ab/='-uhf' .and. ab/='-no') then
   write(iout,'(/,1X,A)') "ERROR in subroutine fch2inporb: the 2nd argument is&
                         & wrong! Only '-uhf' or '-no' is accepted."
   stop
  end if
 end if

 call fch2inporb(fname, ab, sph, uhf)
 call fch2inp_wrap(fname, uhf, .false., 0, 0)

 k = index(fname,'.fch', back=.true.)
 inpname = fname(1:k-1)//'.inp'

 if(sph) then
  i = system('bas_gms2molcas '//TRIM(inpname)//' -sph')
 else
  i = system('bas_gms2molcas '//TRIM(inpname))
 end if

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine fch2inporb: call utility bas_gms2molcas failed.'
  write(iout,'(A)') 'Did you forget to compile bas_gms2molcas? Or the file '//&
                     TRIM(fname)//' may be incomplete.'
  stop
 end if

 open(newunit=i,file=TRIM(inpname),status='old')
 close(unit=i,status='delete')
 stop
end program main

! nbf: the number of basis functions
! nif: the number of independent functions, i.e., the number of MOs

! read the MOs in .fch(k) file and adjust its d,f,g, etc. functions order
!  of Gaussian to that of Molcas
subroutine fch2inporb(fchname, ab, sph, uhf)
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
 character(len=4), intent(in) :: ab
 character(len=8) :: hostname
 character(len=24) :: data_string
 character(len=240), intent(in) :: fchname
 character(len=240) :: orbfile
 ! orbfile is the INPORB file of Molcas/OpenMolcas

 real(kind=8), allocatable :: coeff(:,:), occ_num(:)
 logical, intent(out) :: sph, uhf

 sph = .true. ! default value
 uhf = .false.
 orbfile = ' '

 call read_na_and_nb_from_fch(fchname, nalpha, nbeta)
 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 nbf0 = nbf

 ! read MO Coefficients
 select case(ab)
 case('-uhf')
  uhf = .true.
  allocate(coeff(nbf,2*nif))
  call read_mo_from_fch(fchname, nbf, nif, 'a', coeff(:,1:nif))
  call read_mo_from_fch(fchname, nbf, nif, 'b', coeff(:,nif+1:2*nif))
  nif = 2*nif   ! double the size
 case('-no')
  allocate(occ_num(nif), coeff(nbf,nif))
  call read_eigenvalues_from_fch(fchname, nif, 'a', occ_num)
  call read_mo_from_fch(fchname, nbf, nif, 'a', coeff)
 case default
  allocate(coeff(nbf,nif))
  call read_mo_from_fch(fchname, nbf, nif, 'a', coeff)
 end select

 call read_ncontr_from_fch(fchname, k)
 allocate(shell_type(2*k), source=0)
 allocate(shell2atom_map(2*k), source=0)
 call read_shltyp_and_shl2atm_from_fch(fchname, k, shell_type, shell2atom_map)
 if( ANY(shell_type>1) ) then ! whether Cartesian/spherical harmonic
  sph = .false.
 else
  sph = .true.
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
  end select
 end do

 ! adjust the order of d, f, etc. functions
 do i = 1, n5dmark, 1
  call fch2inporb_permute_5d(nif,coeff(d_mark(i):d_mark(i)+4,:))
 end do
 do i = 1, n6dmark, 1
  call fch2inporb_permute_6d(nif,coeff(d_mark(i):d_mark(i)+5,:))
 end do
 do i = 1, n7fmark, 1
  call fch2inporb_permute_7f(nif,coeff(f_mark(i):f_mark(i)+6,:))
 end do
 do i = 1, n10fmark, 1
  call fch2inporb_permute_10f(nif,coeff(f_mark(i):f_mark(i)+9,:))
 end do
 do i = 1, n9gmark, 1
  call fch2inporb_permute_9g(nif,coeff(g_mark(i):g_mark(i)+8,:))
 end do
 do i = 1, n15gmark, 1
  call fch2inporb_permute_15g(nif,coeff(g_mark(i):g_mark(i)+14,:))
 end do
 do i = 1, n11hmark, 1
  call fch2inporb_permute_11h(nif,coeff(h_mark(i):h_mark(i)+10,:))
 end do
 do i = 1, n21hmark, 1
  call fch2inporb_permute_21h(nif,coeff(h_mark(i):h_mark(i)+20,:))
 end do
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

 select case(ab)
 case('-uhf')
  write(orbid,'(A)') '* UHF orbitals'
  write(orbid,'(A)') '       1       1       4'
  nif = nif/2   ! change to original size
 case('-no')
  write(orbid,'(A)') '* natural orbitals'
  write(orbid,'(A)') '       0       1       0'
 case default ! '-a '
  write(orbid,'(A)') '* SCF orbitals'
  write(orbid,'(A)') '       0       1       2'
 end select

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

 if(ab == '-uhf') then
  write(orbid,'(A)') '#UORB'
  do i = 1, nif, 1
   write(orbid,'(A14,I5)') '* ORBITAL    1', i
   write(orbid,'(5(1X,ES21.14))') (coeff(j,i+nif),j=1,nbf0)
  end do
 end if
 deallocate(coeff)

 select case(ab)
 case('-uhf')
  allocate(occ_num(2*nif), source=0.0d0)
  occ_num(1:nalpha) = 1.0d0
  occ_num(nif+1:nif+nbeta) = 1.0d0
  write(orbid,'(A,/,A)') '#OCC', '* OCCUPATION NUMBERS'
  write(orbid,'(5(1X,ES21.14))') (occ_num(i),i=1,nif)
  write(orbid,'(A,/,A)') '#UOCC', '* Beta OCCUPATION NUMBERS'
  write(orbid,'(5(1X,ES21.14))') (occ_num(i),i=nif+1,2*nif)

 case('-no')
  write(orbid,'(A,/,A)') '#OCC', '* OCCUPATION NUMBERS'
  write(orbid,'(5(1X,ES21.14))') (occ_num(i),i=1,nif)

 case default
  allocate(occ_num(nif), source=0.0d0)
  if(nbeta < nalpha) then
   occ_num(1:nbeta) = 2.0d0
   occ_num(nbeta+1:nalpha) = 1.0d0
  else if(nbeta == nalpha) then
   occ_num(1:nalpha) = 2.0d0
  end if
  write(orbid,'(A,/,A)') '#OCC', '* OCCUPATION NUMBERS'
  write(orbid,'(5(1X,ES21.14))') (occ_num(i),i=1,nif)
 end select

 write(orbid,'(A,/,A)') '#OCHR', '* OCCUPATION NUMBERS (HUMAN-READABLE)'
 write(orbid,'(10(1X,F7.4))') (occ_num(i),i=1,nif)
 if(ab == '-uhf') then
  write(orbid,'(A,/,A)') '#UOCHR', '* Beta OCCUPATION NUMBERS (HUMAN-READABLE)'
  write(orbid,'(10(1X,F7.4))') (occ_num(i),i=nif+1,2*nif)
 end if

 deallocate(occ_num)
 close(orbid)
! print done
 return
end subroutine fch2inporb

! move the 2nd, 3rd, ... Zeta basis functions forward
subroutine zeta_mv_forwd(i0, shell_type, length, nbf, nif, coeff2)
 implicit none
 integer i, j, k
 integer, intent(in) :: i0, shell_type, length, nbf, nif
 integer, parameter :: iout = 6
 integer, parameter :: num0(-5:5) = [11, 9, 7, 5, 0, 0, 3, 6, 10, 15, 21]
 !                                   11H 9G 7F 5D L  S 3P 6D 10F 15G 21H
 real(kind=8), intent(inout) :: coeff2(nbf,nif)
 real(kind=8), allocatable :: coeff(:,:)

 if(length == 1) return

 if(shell_type==0 .or. shell_type==-1) then
  write(iout,'(A)') 'ERROR in subroutine zeta_mv_forwd: this element of&
                   & shell_type is 0 or -1. Impossible.'
  write(iout,'(2(A,I0))') 'shell_type=', shell_type, ', length=', length
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
 return
end subroutine zeta_mv_forwd

! split the 'L' into 'S' and 'P'
subroutine split_L_func(k, shell_type, shell2atom_map, length)
 implicit none
 integer i, k0
 integer,intent(in) :: k
 integer,intent(inout) :: shell_type(2*k), shell2atom_map(2*k)
 integer,intent(out) :: length
 integer,allocatable :: temp1(:), temp2(:)

 k0 = 2*k
 length = k
 ! set initial values for arrays shell_type, assume 15 has not be used
 shell_type(k+1:k0) = 15
 i = 1
 do while(shell_type(i) /= 15)
  if(shell_type(i) /= -1) then
   i = i + 1
   cycle
  end if
  shell_type(i) = 0
  allocate( temp1(i+1 : k0-1), temp2(i+1 : k0-1) )
  temp1(i+1 : k0-1) = shell_type(i+1 : k0-1)
  shell_type(i+2 : k0) = temp1(i+1 : k0-1)
  temp2(i+1 : k0-1) = shell2atom_map(i+1 : k0-1)
  shell2atom_map(i+2 : k0) = temp2(i+1 : k0-1)
  deallocate(temp1, temp2)
  shell_type(i+1) = 1
  shell2atom_map(i+1) = shell2atom_map(i)
  i = i + 2
 end do

 length = i - 1
 shell_type(i : k0) = 0
 return
end subroutine split_L_func

! sort the shell_type, shell2atom_map by ascending order
! MOs will be adjusted accordingly
subroutine sort_shell_and_mo(ilen, shell_type, shell2atom_map, nbf, nif, coeff2)
 implicit none
 integer :: i, j, k
 integer :: ibegin, iend, natom
 integer :: jbegin, jend
 integer, parameter :: ntype = 10
 integer, parameter :: num0(ntype) = [0, 1, -2, 2, -3, 3, -4, 4, -5, 5]
 integer, parameter :: num1(ntype) = [1, 3, 5, 6, 7, 10, 9, 15, 11, 21]
 !                                    S  P  5D 6D 7F 10F 9G 15G 11H 21H
 integer num(ntype)

 integer, intent(in) :: ilen, nbf, nif
 integer, intent(inout) :: shell_type(ilen), shell2atom_map(ilen)
 integer, allocatable :: ith(:), ith_bas(:), tmp_type(:)
 real(kind=8), intent(inout) :: coeff2(nbf,nif)

 ! find the number of atoms
 natom = shell2atom_map(ilen)

 allocate(ith(0:natom), ith_bas(0:natom))
 ith = 0
 ith_bas = 0

 ! find the end position of each atom in array shell2atom_map
 do i = 1, natom, 1
  ith(i) = count(shell2atom_map==i) + ith(i-1)
 end do

 ! find the end position of basis functions between two atoms
 do i = 1, natom, 1
  ibegin = ith(i-1) + 1
  iend = ith(i)
  allocate(tmp_type(ibegin:iend))
  tmp_type = 0
  tmp_type = shell_type(ibegin:iend)
  num = 0
  do j = 1, ntype, 1
   num(j) = count(tmp_type == num0(j))
  end do
  k = 0
  do j = 1, ntype, 1
   k = k + num(j)*num1(j)
  end do
  ith_bas(i) = ith_bas(i-1) + k
  deallocate(tmp_type)
 end do

 ! adjust the MOs in each atom
 do i = 1, natom, 1
  ibegin = ith(i-1) + 1
  iend = ith(i)
  jbegin = ith_bas(i-1) + 1
  jend = ith_bas(i)
  call sort_shell_and_mo_in_each_atom(iend-ibegin+1, shell_type(ibegin:iend), &
  & jend-jbegin+1, nif, coeff2(jbegin:jend,1:nif))
 end do
 ! adjust the MOs in each atom done
 deallocate(ith, ith_bas)
 return
end subroutine sort_shell_and_mo

subroutine sort_shell_and_mo_in_each_atom(ilen1, shell_type, ilen2, nif, coeff2)
 implicit none
 integer i, tmp_type
 integer ibegin, iend, jbegin, jend
 integer,parameter :: ntype = 10
 integer,parameter :: num0(ntype) = [0, 1, -2, 2, -3, 3, -4, 4, -5, 5]
 integer,parameter :: num1(ntype) = [1, 3, 5, 6, 7, 10, 9, 15, 11, 21]
 !                                   S  P  5D 6D 7F 10F 9G 15G 11H 21H
 integer,parameter :: rnum(-5:5) = [9, 7, 5, 3, 0, 1, 2, 4, 6, 8, 10]

 integer,intent(in) :: nif, ilen1, ilen2
 integer,intent(inout) :: shell_type(ilen1)
 integer,allocatable :: ith_bas(:)
 real(kind=8),intent(inout) :: coeff2(ilen2,nif)
 real(kind=8),allocatable :: tmp_coeff1(:,:), tmp_coeff2(:,:)
 logical sort_done

 tmp_type = 0

 ! find the end position of basis functions within an atom
 allocate(ith_bas(0:ilen1))
 ith_bas = 0
 do i = 1, ilen1, 1
   ith_bas(i) = ith_bas(i-1) + num1(rnum(shell_type(i)))
 end do

 sort_done = .false.
 do while(.not. sort_done)
  sort_done = .true.
  do i = 1, ilen1-1, 1
    if(shell_type(i) == 0) cycle
    if(ABS(shell_type(i+1)) >= ABS(shell_type(i))) cycle
    sort_done = .false.
    tmp_type = shell_type(i+1)
    shell_type(i+1) = shell_type(i)
    shell_type(i) = tmp_type
    ibegin = ith_bas(i-1) + 1
    iend = ith_bas(i)
    jbegin = ith_bas(i) + 1
    jend = ith_bas(i+1)
    allocate(tmp_coeff1(ibegin:iend,nif), tmp_coeff2(jbegin:jend,nif))
    tmp_coeff1 = 0.0d0
    tmp_coeff2 = 0.0d0
    tmp_coeff1 = coeff2(ibegin:iend,:)
    tmp_coeff2 = coeff2(jbegin:jend,:)
    ith_bas(i) = ibegin+jend-jbegin
    coeff2(ibegin: ith_bas(i),:) = tmp_coeff2
    ith_bas(i+1) = jend+iend-jbegin+1
    coeff2(ibegin+jend-jbegin+1: ith_bas(i+1),:) = tmp_coeff1
    deallocate(tmp_coeff1, tmp_coeff2)
  end do
 end do
 deallocate(ith_bas)
 return
end subroutine sort_shell_and_mo_in_each_atom

subroutine fch2inporb_permute_5d(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(5) = [5, 3, 1, 2, 4]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(5,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical d functions in Gaussian
! To: the order of spherical d functions in Molcas
! 1    2    3    4    5
! d0 , d+1, d-1, d+2, d-2
! d-2, d-1, d0 , d+1, d+2

 allocate(coeff2(5,nif), source=0.0d0)
 forall(i = 1:5) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
 return
end subroutine fch2inporb_permute_5d

subroutine fch2inporb_permute_6d(nif,coeff)
 use root_parameter, only: root3
 implicit none
 integer :: i, j
 integer, parameter :: order(6) = [1, 4, 5, 2, 6, 3]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(6,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian d functions in Gaussian
! To: the order of Cartesian d functions in Molcas
! 1  2  3  4  5  6
! XX,YY,ZZ,XY,XZ,YZ
! XX,XY,XZ,YY,YZ,ZZ

 forall(j=1:3 , i=1:nif) coeff(j,i) = coeff(j,i)/root3

 allocate(coeff2(6,nif), source=0.0d0)
 forall(i = 1:6) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
 return
end subroutine fch2inporb_permute_6d

subroutine fch2inporb_permute_7f(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(7) = [7, 5, 3, 1, 2, 4, 6]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(7,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical f functions in Gaussian
! To: the order of spherical f functions in Molcas
! 1    2    3    4    5    6    7
! f0 , f+1, f-1, f+2, f-2, f+3, f-3
! f-3, f-2, f-1, f0 , f+1, f+2, f+3

 allocate(coeff2(7,nif), source=0.0d0)
 forall(i = 1:7) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
 return
end subroutine fch2inporb_permute_7f

subroutine fch2inporb_permute_10f(nif,coeff)
 use root_parameter, only: root3, root15
 implicit none
 integer :: i, j
 integer, parameter :: order(10) = [1, 5, 6, 4, 10, 7, 2, 9, 8, 3]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(10,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian f functions in Gaussian
! To: the order of Cartesian f functions in Molcas
! 1   2   3   4   5   6   7   8   9   10
! XXX,YYY,ZZZ,XYY,XXY,XXZ,XZZ,YZZ,YYZ,XYZ
! XXX,XXY,XXZ,XYY,XYZ,XZZ,YYY,YYZ,YZZ,ZZZ

 forall(j=1:3 , i=1:nif) coeff(j,i) = coeff(j,i)/root15
 forall(j=4:9 , i=1:nif) coeff(j,i) = coeff(j,i)/root3

 allocate(coeff2(10,nif), source=0.0d0)
 forall(i = 1:10) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
 return
end subroutine fch2inporb_permute_10f

subroutine fch2inporb_permute_9g(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(9) = [9, 7, 5, 3, 1, 2, 4, 6, 8]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(9,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical g functions in Gaussian
! To: the order of spherical g functions in Molcas
! 1    2    3    4    5    6    7    8    9
! g0 , g+1, g-1, g+2, g-2, g+3, g-3, g+4, g-4
! g-4, g-3, g-2, g-1, g0 , g+1, g+2, g+3, g+4

 allocate(coeff2(9,nif), source=0.0d0)
 forall(i = 1:9) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
 return
end subroutine fch2inporb_permute_9g

subroutine fch2inporb_permute_15g(nif,coeff)
 use root_parameter, only: root3, root9, root15, root105
 implicit none
 integer :: i, j
 integer, intent(in) :: nif
 real(kind=8), parameter :: ratio(15) = [root105, root15, root9, root15, root105, &
  root15, root3, root3, root15, root9, root3, root9, root15, root15, root105]
 real(kind=8), intent(inout) :: coeff(15,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian g functions in Gaussian
! To: the order of Cartesian g functions in Molcas
! 1    2    3    4    5    6    7    8    9    10   11   12   13   14   15
! ZZZZ,YZZZ,YYZZ,YYYZ,YYYY,XZZZ,XYZZ,XYYZ,XYYY,XXZZ,XXYZ,XXYY,XXXZ,XXXY,XXXX
! xxxx,xxxy,xxxz,xxyy,xxyz,xxzz,xyyy,xyyz,xyzz,xzzz,yyyy,yyyz,yyzz,yzzz,zzzz

 forall(j=1:15 , i=1:nif) coeff(j,i) = coeff(j,i)/ratio(j)

 allocate(coeff2(15,nif), source=0.0d0)
 forall(i = 1:15) coeff2(i,:) = coeff(16-i,:)
 coeff = coeff2
 deallocate(coeff2)
 return
end subroutine fch2inporb_permute_15g

subroutine fch2inporb_permute_11h(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(11) = [11, 9, 7, 5, 3, 1, 2, 4, 6, 8, 10]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(11,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical h functions in Gaussian
! To: the order of spherical h functions in Molcas
! 1    2    3    4    5    6    7    8    9    10   11
! h0 , h+1, h-1, h+2, h-2, h+3, h-3, h+4, h-4, h+5, h-5
! h-5, h-4, h-3, h-2, h-1, h0 , h+1, h+2, h+3, h+4, h+5

 allocate(coeff2(11,nif), source=0.0d0)
 forall(i = 1:11) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
 return
end subroutine fch2inporb_permute_11h

subroutine fch2inporb_permute_21h(nif,coeff)
 use root_parameter
 implicit none
 integer :: i, j
 integer, intent(in) :: nif
 real(kind=8), parameter :: ratio(21) = [root945, root105, root45, root45, root105, &
  root945, root105, root15, root9, root15, root105, root45, root9, root9, root45, &
  root45, root15, root45, root105, root105, root945]
 real(kind=8), intent(inout) :: coeff(21,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian h functions in Gaussian
! To: the order of Cartesian h functions in Molcas
! 1     2     3     4     5     6     7     8     9     10    11    12    13    14    15    16    17    18    19    20    21
! ZZZZZ,YZZZZ,YYZZZ,YYYZZ,YYYYZ,YYYYY,XZZZZ,XYZZZ,XYYZZ,XYYYZ,XYYYY,XXZZZ,XXYZZ,XXYYZ,XXYYY,XXXZZ,XXXYZ,XXXYY,XXXXZ,XXXXY,XXXXX
! xxxxx,xxxxy,xxxxz,xxxyy,xxxyz,xxxzz,xxyyy,xxyyz,xxyzz,xxzzz,xyyyy,xyyyz,xyyzz,xyzzz,xzzzz,yyyyy,yyyyz,yyyzz,yyzzz,yzzzz,zzzzz

 forall(j=1:21 , i=1:nif) coeff(j,i) = coeff(j,i)/ratio(j)

 allocate(coeff2(21,nif), source=0.0d0)
 forall(i = 1:21) coeff2(i,:) = coeff(22-i,:)
 coeff = coeff2
 deallocate(coeff2)
 return
end subroutine fch2inporb_permute_21h

