! written by jxzou at 20200907: adjust the orders of d,f,g, etc functions in
!  OpenMolcas, to that in Gaussian .fch(k) file
! Originally copied from fch2inporb.f90, some modifications are made
! updated by jxzou at 20210413: remove '-uhf', add automatic determination

module root_parameter
 implicit none
 real(kind=8), parameter :: root3   = DSQRT(3d0)   ! SQRT(3)
 real(kind=8), parameter :: root9   = 3d0          ! SQRT(9)
 real(kind=8), parameter :: root15  = DSQRT(15d0)  ! SQRT(15)
 real(kind=8), parameter :: root45  = DSQRT(45d0)  ! SQRT(45)
 real(kind=8), parameter :: root105 = DSQRT(105d0) ! SQRT(105)
 real(kind=8), parameter :: root945 = DSQRT(945d0) ! SQRT(945)
end module

program main
 implicit none
 integer :: i
 integer, parameter :: iout = 6
 character(len=3) :: str
 character(len=240) :: fchname, orbname
 logical :: prt_no

 i = iargc()
 if(.not. (i==2 .or. i==3)) then
  write(iout,'(/,A)') ' ERROR in subroutine orb2fch: wrong command line arguments!'
  write(iout,'(A)')   ' Example 1 (for RHF)   : orb2fch a.ScfOrb a.fch'
  write(iout,'(A)')   ' Example 2 (for UHF)   : orb2fch a.UhfOrb a.fch'
  write(iout,'(A)')   ' Example 3 (for CAS)   : orb2fch a.RasOrb a.fch'
  write(iout,'(A)')   ' Example 4 (for UNO)   : orb2fch a.UnaOrb a.fch -no'
  write(iout,'(A,/)') ' Example 5 (for CAS NO): orb2fch a.RasOrb.1 a.fch -no'
  stop
 end if

 orbname = ' '
 fchname = ' '
 call getarg(1,orbname)
 call require_file_exist(orbname)

 call getarg(2,fchname)
 call require_file_exist(fchname)
 prt_no = .false.

 if(i == 3) then
  call getarg(3, str)
  if(str /= '-no') then
   write(iout,'(/,A)') "ERROR in subroutine orb2fch: the 3rd argument is&
                      & wrong! Only '-no' is accepted."
   stop
  else
   prt_no = .true.
  end if
 end if

 call orb2fch(orbname, fchname, prt_no)
 stop
end program main

! read the MOs in orbital file of OpenMolcas and adjust its d,f,g,h functions
!  order to that of Gaussian
subroutine orb2fch(orbname, fchname, prt_no)
 use fch_content, only: iout, check_uhf_in_fch
 implicit none
 integer :: i, j, k, m, length
 integer :: na, nb, nbf, nif, nbf0, nbf1
 integer :: n6dmark, n10fmark, n15gmark, n21hmark
 integer :: n5dmark, n7fmark, n9gmark, n11hmark
 integer, allocatable :: shell_type(:), shell2atom_map(:)
 integer, allocatable :: idx(:), idx2(:)
 ! mark the index where d, f, g, h functions begin
 integer, allocatable :: d_mark(:), f_mark(:), g_mark(:), h_mark(:)
 character(len=240), intent(in) :: orbname, fchname
 ! orbname is one of .ScfOrb, .RasOrb, .RasOrb.1, .UnaOrb, .UhfOrb file of OpenMolcas
 real(kind=8), allocatable :: coeff(:,:), coeff2(:,:), occ_num(:), norm(:)
 logical :: uhf
 logical, intent(in) :: prt_no

 call check_uhf_in_fch(fchname, uhf)
 call read_na_and_nb_from_fch(fchname, na, nb)
 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 nbf0 = nbf ! make a copy of nbf

 ! read MO Coefficients from .ScfOrb, .RasOrb, .RasOrb.1, .UnaOrb, or .UhfOrb
 ! file of OpenMolcas
 if(uhf) then
  allocate(coeff(nbf,2*nif))
  call read_mo_from_orb(orbname, nbf, nif, 'a', coeff(:,1:nif))
  call read_mo_from_orb(orbname, nbf, nif, 'b', coeff(:,nif+1:2*nif))
  nif = 2*nif   ! double the size
 else
  if(prt_no) then
   allocate(occ_num(nif), coeff(nbf,nif))
   call read_on_from_orb(orbname, nif, 'a', occ_num)
   call read_mo_from_orb(orbname, nbf, nif, 'a', coeff)
  else
   allocate(coeff(nbf,nif))
   call read_mo_from_orb(orbname, nbf, nif, 'a', coeff)
  end if
 end if

 call read_ncontr_from_fch(fchname, k)
 allocate(shell_type(2*k), source=0)
 allocate(shell2atom_map(2*k), source=0)
 call read_shltyp_and_shl2atm_from_fch(fchname, k, shell_type, shell2atom_map)

 if(ANY(shell_type>1) .and. ANY(shell_type<-2)) then
  write(iout,'(A)') 'ERROR in subroutine orb2fch: mixed Cartesian/spherical harmonic&
                   & functions detected. Cannot deal with that.'
  write(iout,'(A)') 'One possible reason is that you used a Pople-type basis set in Gaussian.'
  write(iout,'(A)') "Its default setting is '6D 7F', you should add keywords '5D 7F' or '6D 10F'."
  stop
 end if

! first we adjust the basis functions in each MO according to the Shell to atom map
 ! 1) split the 'L' into 'S' and 'P', this is to ensure that D comes after L functions
 call split_L_func(k, shell_type, shell2atom_map, length)

 allocate(idx(nbf))
 forall(i = 1:nbf) idx(i) = i
 allocate(norm(nbf), source=1.0d0)

 ! 2) sort the shell_type and shell2atom_map by ascending order,
 !  the indices of MOs will be adjusted accordingly
 call sort_shell_and_mo(length, shell_type, shell2atom_map, nbf, idx)
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
  call orb2fch_permute_5d(idx(d_mark(i):d_mark(i)+4))
 end do
 do i = 1, n6dmark, 1
  call orb2fch_permute_6d(idx(d_mark(i):d_mark(i)+5), norm(d_mark(i):d_mark(i)+5))
 end do
 do i = 1, n7fmark, 1
  call orb2fch_permute_7f(idx(f_mark(i):f_mark(i)+6))
 end do
 do i = 1, n10fmark, 1
  call orb2fch_permute_10f(idx(f_mark(i):f_mark(i)+9), norm(f_mark(i):f_mark(i)+9))
 end do
 do i = 1, n9gmark, 1
  call orb2fch_permute_9g(idx(g_mark(i):g_mark(i)+8))
 end do
 do i = 1, n15gmark, 1
  call orb2fch_permute_15g(idx(g_mark(i):g_mark(i)+14), norm(g_mark(i):g_mark(i)+14))
 end do
 do i = 1, n11hmark, 1
  call orb2fch_permute_11h(idx(h_mark(i):h_mark(i)+10))
 end do
 do i = 1, n21hmark, 1
  call orb2fch_permute_21h(idx(h_mark(i):h_mark(i)+20), norm(h_mark(i):h_mark(i)+20))
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
  end do ! for while

  if(i < k) i = i - 1
  if(length > 1) then
   call zeta_mv_forwd(nbf1, m, length, nbf0, idx, norm)
   nbf = nbf1 + length*(nbf-nbf1)
  end if
 end do ! for while

 deallocate(shell_type, shell2atom_map)
! move done

 nbf = nbf0
 allocate(idx2(nbf), coeff2(nbf,nif))
 forall(i = 1:nbf) idx2(idx(i)) = i
 forall(i=1:nbf, j=1:nif) coeff2(i,j) = coeff(idx2(i),j)*norm(idx2(i))
 deallocate(idx, idx2, norm, coeff)

! print MOs into .fch(k) file
 if(uhf) then
  nif = nif/2
  call write_mo_into_fch(fchname, nbf, nif, 'a', coeff2(:,1:nif))
  call write_mo_into_fch(fchname, nbf, nif, 'b', coeff2(:,nif+1:2*nif))
 else
  if(prt_no) then
   call write_eigenvalues_to_fch(fchname, nif, 'a', occ_num, .true.)
   call write_mo_into_fch(fchname, nbf, nif, 'a', coeff2)
  else
   call write_mo_into_fch(fchname, nbf, nif, 'a', coeff2)
  end if
 end if
! print done

 deallocate(coeff2)
 if(allocated(occ_num)) deallocate(occ_num)
 return
end subroutine orb2fch

! move the 2nd, 3rd, ... Zeta basis functions forward
subroutine zeta_mv_forwd(i0, shell_type, length, nbf, idx2, norm1)
 implicit none
 integer :: i, j, k
 integer, intent(in) :: i0, shell_type, length, nbf
 integer, parameter :: iout = 6
 integer, parameter :: num0(-5:5) = [11, 9, 7, 5, 0, 0, 3, 6, 10, 15, 21]
 !                                   11H 9G 7F 5D L  S 3P 6D 10F 15G 21H
 integer, intent(inout) :: idx2(nbf)
 integer, allocatable :: idx(:)
 real(kind=8), allocatable :: norm(:)
 real(kind=8), intent(inout) :: norm1(nbf)

 if(length == 1) return

 if(shell_type==0 .or. shell_type==-1) then
  write(iout,'(A)') 'ERROR in subroutine zeta_mv_forwd: this element of&
                   & shell_type is 0 or -1. Impossible.'
  write(iout,'(2(A,I0))') 'shell_type=', shell_type, ', length=', length
  stop
 end if

 idx = idx2
 norm = norm1
 k = num0(shell_type)

 do i = 1, k, 1
  do j = 1, length, 1
   idx(i0+j+(i-1)*length) = idx2(i0+i+(j-1)*k)
   norm(i0+j+(i-1)*length) = norm1(i0+i+(j-1)*k)
  end do ! for j
 end do ! for i

 idx2 = idx
 norm1 = norm
 deallocate(idx, norm)
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

! sort the shell_type, shell2atom_map by ascending order,
! indices of MOs will be adjusted accordingly
subroutine sort_shell_and_mo(ilen, shell_type, shell2atom_map, nbf, idx)
 implicit none
 integer :: i, j, k
 integer :: ibegin, iend, natom
 integer :: jbegin, jend
 integer, parameter :: ntype = 10
 integer, parameter :: num0(ntype) = [0, 1, -2, 2, -3, 3, -4, 4, -5, 5]
 integer, parameter :: num1(ntype) = [1, 3, 5, 6, 7, 10, 9, 15, 11, 21]
 !                                    S  P  5D 6D 7F 10F 9G 15G 11H 21H
 integer :: num(ntype)
 integer, intent(in) :: ilen, nbf
 integer, intent(inout) :: shell_type(ilen), shell2atom_map(ilen)
 integer, intent(inout) :: idx(nbf)
 integer, allocatable :: ith(:), ith_bas(:), tmp_type(:)

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
  allocate(tmp_type(ibegin:iend), source=0)
  tmp_type = shell_type(ibegin:iend)
  num = 0

  do j = 1, ntype, 1
   num(j) = count(tmp_type == num0(j))
  end do ! for j

  k = 0
  do j = 1, ntype, 1
   k = k + num(j)*num1(j)
  end do ! for j

  ith_bas(i) = ith_bas(i-1) + k
  deallocate(tmp_type)
 end do ! for i

 ! adjust the indices of MOs in each atom
 do i = 1, natom, 1
  ibegin = ith(i-1) + 1
  iend = ith(i)
  jbegin = ith_bas(i-1) + 1
  jend = ith_bas(i)
  call sort_shell_and_mo_in_each_atom(iend-ibegin+1, shell_type(ibegin:iend), &
  & jend-jbegin+1, idx(jbegin:jend))
 end do ! for i

 deallocate(ith, ith_bas)
 return
end subroutine sort_shell_and_mo

subroutine sort_shell_and_mo_in_each_atom(ilen1, shell_type, ilen2, idx)
 implicit none
 integer :: i, tmp_type, ibegin, iend, jbegin, jend
 integer, parameter :: ntype = 10
 integer, parameter :: num0(ntype) = [0, 1, -2, 2, -3, 3, -4, 4, -5, 5]
 integer, parameter :: num1(ntype) = [1, 3, 5, 6, 7, 10, 9, 15, 11, 21]
 !                                   S  P  5D 6D 7F 10F 9G 15G 11H 21H
 integer, parameter :: rnum(-5:5) = [9, 7, 5, 3, 0, 1, 2, 4, 6, 8, 10]
 integer, intent(in) :: ilen1, ilen2
 integer, intent(inout) :: shell_type(ilen1), idx(ilen2)
 integer, allocatable :: ith_bas(:), idx1(:), idx2(:)
 logical :: sort_done

 tmp_type = 0

 ! find the end position of basis functions within an atom
 allocate(ith_bas(0:ilen1), source=0)
 do i = 1, ilen1, 1
  ith_bas(i) = ith_bas(i-1) + num1(rnum(shell_type(i)))
 end do ! for i

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
    ibegin = ith_bas(i-1) + 1; iend = ith_bas(i)
    jbegin = ith_bas(i) + 1  ; jend = ith_bas(i+1)
    allocate(idx1(ibegin:iend), source=idx(ibegin:iend))
    allocate(idx2(jbegin:jend), source=idx(jbegin:jend))
    ith_bas(i) = ibegin + jend - jbegin
    idx(ibegin:ith_bas(i)) = idx2
    ith_bas(i+1) = jend + iend - jbegin + 1
    idx(ith_bas(i)+1: ith_bas(i+1)) = idx1
    deallocate(idx1, idx2)
  end do ! for i
 end do ! for while

 deallocate(ith_bas)
 return
end subroutine sort_shell_and_mo_in_each_atom

subroutine orb2fch_permute_5d(idx)
 implicit none
 integer :: i, idx0(5)
 integer, parameter :: order(5) = [5, 3, 1, 2, 4]
 integer, intent(inout) :: idx(5)
! From: the order of spherical d functions in Gaussian
! To: the order of spherical d functions in Molcas
! 1    2    3    4    5
! d0 , d+1, d-1, d+2, d-2
! d-2, d-1, d0 , d+1, d+2

 idx0 = idx
 forall(i = 1:5) idx(i) = idx0(order(i))
 return
end subroutine orb2fch_permute_5d

subroutine orb2fch_permute_6d(idx, norm)
 use root_parameter, only: root3
 implicit none
 integer :: i, idx0(6)
 integer, parameter :: order(6) = [1, 4, 5, 2, 6, 3]
 integer, intent(inout) :: idx(6)
 real(kind=8) :: norm0(6)
 real(kind=8), intent(inout) :: norm(6)
! From: the order of Cartesian d functions in Gaussian
! To: the order of Cartesian d functions in Molcas
! 1  2  3  4  5  6
! XX,YY,ZZ,XY,XZ,YZ
! XX,XY,XZ,YY,YZ,ZZ

 forall(i=1:3) norm(i) = norm(i)*root3
 norm0 = norm
 idx0 = idx

 forall(i = 1:6)
  idx(i) = idx0(order(i))
  norm(i) = norm0(order(i))
 end forall

 return
end subroutine orb2fch_permute_6d

subroutine orb2fch_permute_7f(idx)
 implicit none
 integer :: i, idx0(7)
 integer, parameter :: order(7) = [7, 5, 3, 1, 2, 4, 6]
 integer, intent(inout) :: idx(7)
! From: the order of spherical f functions in Gaussian
! To: the order of spherical f functions in Molcas
! 1    2    3    4    5    6    7
! f0 , f+1, f-1, f+2, f-2, f+3, f-3
! f-3, f-2, f-1, f0 , f+1, f+2, f+3

 idx0 = idx
 forall(i = 1:7) idx(i) = idx0(order(i))
 return
end subroutine orb2fch_permute_7f

subroutine orb2fch_permute_10f(idx, norm)
 use root_parameter, only: root3, root15
 implicit none
 integer :: i, idx0(10)
 integer, parameter :: order(10) = [1, 5, 6, 4, 10, 7, 2, 9, 8, 3]
 integer, intent(inout) :: idx(10)
 real(kind=8) :: norm0(10)
 real(kind=8), intent(inout) :: norm(10)
! From: the order of Cartesian f functions in Gaussian
! To: the order of Cartesian f functions in Molcas
! 1   2   3   4   5   6   7   8   9   10
! XXX,YYY,ZZZ,XYY,XXY,XXZ,XZZ,YZZ,YYZ,XYZ
! XXX,XXY,XXZ,XYY,XYZ,XZZ,YYY,YYZ,YZZ,ZZZ

 forall(i=1:3) norm(i) = norm(i)*root15
 forall(i=4:9) norm(i) = norm(i)*root3
 norm0 = norm
 idx0 = idx

 forall(i = 1:10)
  idx(i) = idx0(order(i))
  norm(i) = norm0(order(i))
 end forall
 return
end subroutine orb2fch_permute_10f

subroutine orb2fch_permute_9g(idx)
 implicit none
 integer :: i, idx0(9)
 integer, parameter :: order(9) = [9, 7, 5, 3, 1, 2, 4, 6, 8]
 integer, intent(inout) :: idx(9)
! From: the order of spherical g functions in Gaussian
! To: the order of spherical g functions in Molcas
! 1    2    3    4    5    6    7    8    9
! g0 , g+1, g-1, g+2, g-2, g+3, g-3, g+4, g-4
! g-4, g-3, g-2, g-1, g0 , g+1, g+2, g+3, g+4

 idx0 = idx
 forall(i = 1:9) idx(i) = idx0(order(i))
 return
end subroutine orb2fch_permute_9g

subroutine orb2fch_permute_15g(idx, norm)
 use root_parameter, only: root3, root9, root15, root105
 implicit none
 integer :: i, idx0(15)
 integer, intent(inout) :: idx(15)
 real(kind=8) :: norm0(15)
 real(kind=8), parameter :: ratio(15) = [root105, root15, root9, root15, root105, &
  root15, root3, root3, root15, root9, root3, root9, root15, root15, root105]
 real(kind=8), intent(inout) :: norm(15)
! From: the order of Cartesian g functions in Gaussian
! To: the order of Cartesian g functions in Molcas
! 1    2    3    4    5    6    7    8    9    10   11   12   13   14   15
! ZZZZ,YZZZ,YYZZ,YYYZ,YYYY,XZZZ,XYZZ,XYYZ,XYYY,XXZZ,XXYZ,XXYY,XXXZ,XXXY,XXXX
! xxxx,xxxy,xxxz,xxyy,xxyz,xxzz,xyyy,xyyz,xyzz,xzzz,yyyy,yyyz,yyzz,yzzz,zzzz

 forall(i=1:15) norm(i) = norm(i)*ratio(i)
 norm0 = norm
 idx0 = idx

 forall(i = 1:15)
  idx(i) = idx0(16-i)
  norm(i) = norm0(16-i)
 end forall
 return
end subroutine orb2fch_permute_15g

subroutine orb2fch_permute_11h(idx)
 implicit none
 integer :: i, idx0(11)
 integer, parameter :: order(11) = [11, 9, 7, 5, 3, 1, 2, 4, 6, 8, 10]
 integer, intent(inout) :: idx(11)
! From: the order of spherical h functions in Gaussian
! To: the order of spherical h functions in Molcas
! 1    2    3    4    5    6    7    8    9    10   11
! h0 , h+1, h-1, h+2, h-2, h+3, h-3, h+4, h-4, h+5, h-5
! h-5, h-4, h-3, h-2, h-1, h0 , h+1, h+2, h+3, h+4, h+5

 idx0 = idx
 forall(i = 1:11) idx(i) = idx0(order(i))
 return
end subroutine orb2fch_permute_11h

subroutine orb2fch_permute_21h(idx, norm)
 use root_parameter
 implicit none
 integer :: i, idx0(21)
 integer, intent(inout) :: idx(21)
 real(kind=8) :: norm0(21)
 real(kind=8), parameter :: ratio(21) = [root945, root105, root45, root45, root105, &
  root945, root105, root15, root9, root15, root105, root45, root9, root9, root45, &
  root45, root15, root45, root105, root105, root945]
 real(kind=8), intent(inout) :: norm(21)
! From: the order of Cartesian h functions in Gaussian
! To: the order of Cartesian h functions in Molcas
! 1     2     3     4     5     6     7     8     9     10    11    12    13    14    15    16    17    18    19    20    21
! ZZZZZ,YZZZZ,YYZZZ,YYYZZ,YYYYZ,YYYYY,XZZZZ,XYZZZ,XYYZZ,XYYYZ,XYYYY,XXZZZ,XXYZZ,XXYYZ,XXYYY,XXXZZ,XXXYZ,XXXYY,XXXXZ,XXXXY,XXXXX
! xxxxx,xxxxy,xxxxz,xxxyy,xxxyz,xxxzz,xxyyy,xxyyz,xxyzz,xxzzz,xyyyy,xyyyz,xyyzz,xyzzz,xzzzz,yyyyy,yyyyz,yyyzz,yyzzz,yzzzz,zzzzz

 forall(i=1:21) norm(i) = norm(i)*ratio(i)
 norm0 = norm
 idx0 = idx

 forall(i = 1:21)
  idx(i) = idx0(22-i)
  norm(i) = norm0(22-i)
 end forall
 return
end subroutine orb2fch_permute_21h

