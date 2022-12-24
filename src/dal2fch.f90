! written by jxzou at 20210615: adjust the orders of d, f, g, h functions etc.
!  in the input file of Dalton, to the oeders of .fch(k) file in Gaussian

! diagonal elements of overlap matrix using Cartesian functions (6D 10F)
module Sdiag_dalton
 implicit none
 real(kind=8), parameter :: PI = 4d0*DATAN(1d0)
 real(kind=8), parameter :: p1 = 2d0*DSQRT(PI/15d0)
 real(kind=8), parameter :: p2 = 2d0*DSQRT(PI/5d0)
 real(kind=8), parameter :: p3 = 2d0*DSQRT(PI/7d0)
 real(kind=8), parameter :: p4 = 2d0*DSQRT(PI/35d0)
 real(kind=8), parameter :: p5 = 2d0*DSQRT(PI/105d0)
 real(kind=8), parameter :: p6 = (2d0/3d0)*DSQRT(PI)
 real(kind=8), parameter :: p7 = (2d0/3d0)*DSQRT(PI/7d0)
 real(kind=8), parameter :: p8 = (2d0/3d0)*DSQRT(PI/35d0)
 real(kind=8), parameter :: p9 = 2d0*DSQRT(PI/11d0)
 real(kind=8), parameter :: p10 = (2d0/3d0)*DSQRT(PI/11d0)
 real(kind=8), parameter :: p11 = 2d0*DSQRT(PI/231d0)
 real(kind=8), parameter :: p12 = (2d0/3d0)*DSQRT(PI/77d0)
 real(kind=8), parameter :: p13 = 2d0*DSQRT(PI/1155d0)
 real(kind=8), parameter :: Sdiag_d(6)  = [p2,p1,p1,p2,p1,p2]
 real(kind=8), parameter :: Sdiag_f(10) = [p3,p4,p4,p4,p5,p4,p3,p4,p4,p3]
 real(kind=8), parameter :: Sdiag_g(15) = [p6,p7,p7,p5,p8,p5,p7,p8,p8,p7,p6,p7,p5,p7,p6]
 real(kind=8), parameter :: Sdiag_h(21) = &
  [p9,p10,p10,p11,p12,p11,p11,p13,p13,p11,p10,p12,p13,p12,p10,p9,p10,p11,p11,p10,p9]
end module Sdiag_dalton

program main
 implicit none
 integer :: i
 character(len=3) :: str
 character(len=240) :: orbname, fchname
 logical :: prt_no

 i = iargc()
 if(i<2 .or. i>3) then
  write(6,'(/,A)') ' ERROR in subroutine dal2fch: wrong command line arguments!'
  write(6,'(A)')   ' Example 1 (R(O)HF, CAS): dal2fch DALTON.MOPUN a.fch'
  write(6,'(A,/)') ' Example 2 (CAS NO)     : dal2fch DALTON.MOPUN a.fch -no'
  stop
 end if

 call getarg(1, orbname)
 call require_file_exist(orbname)
 call getarg(2, fchname)
 call require_file_exist(fchname)

 prt_no = .false.
 if(i == 3) then
  call getarg(3, str)
  if(str /= '-no') then
   write(6,'(/,1X,A)') "ERROR in subroutine dal2fch: the 3rd argument is&
                         & wrong! Only '-no' is accepted."
   stop
  else
   prt_no = .true.
  end if
 end if

 call dal2fch(orbname, fchname, prt_no)
end program main

! transfer MOs from Dalton DALTON.MOPUN file to Gaussian .fch(k) file
subroutine dal2fch(orbname, fchname, prt_no)
 implicit none
 integer :: i, j, nbf, nif
 integer, allocatable :: idx(:)
 real(kind=8), allocatable :: norm(:), coeff(:,:), coeff2(:,:)
 character(len=240), intent(in) :: orbname, fchname
 logical, intent(in) :: prt_no

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(idx(nbf), norm(nbf), coeff(nbf,nif))
 call get_permute_idx_from_fch(fchname, nbf, idx, norm)

 call read_mo_from_dalton_mopun(orbname, nbf, nif, coeff)
 allocate(coeff2(nbf,nif), source=coeff)
 forall(i=1:nbf,j=1:nif) coeff(i,j) = coeff2(idx(i),j)*norm(i)
 deallocate(coeff2, norm, idx)

 call write_mo_into_fch(fchname, nbf, nif, 'a', coeff)
 deallocate(coeff)

 if(prt_no) then
  allocate(norm(nif))
  call read_on_from_dalton_mopun(orbname, nif, norm)
  call write_eigenvalues_to_fch(fchname, nif, 'a', norm, .true.)
  deallocate(norm)
 end if
end subroutine dal2fch

! get permutation index list from a given .fch(k) file
subroutine get_permute_idx_from_fch(fchname, nbf, idx, norm)
 implicit none
 integer :: k
 integer, intent(in) :: nbf
 integer, intent(out) :: idx(nbf)
 integer, allocatable :: shell_type(:), shell_to_atom_map(:)
 real(kind=8), intent(out) :: norm(nbf)
 character(len=240), intent(in) :: fchname

 call read_ncontr_from_fch(fchname, k)
 allocate(shell_type(k), shell_to_atom_map(k))
 call read_shltyp_and_shl2atm_from_fch(fchname, k, shell_type, shell_to_atom_map)

 call get_permute_idx_from_shell(k, shell_type, shell_to_atom_map, nbf, idx, norm)
 deallocate(shell_type, shell_to_atom_map)
end subroutine get_permute_idx_from_fch

! get permutation index list from two arrays (shell_type and shell_to_atom_map)
subroutine get_permute_idx_from_shell(ncontr, shell_type0, shell_to_atom_map0, nbf0, idx, norm)
 implicit none
 integer :: i, k, nbf
 integer, intent(in) :: ncontr, nbf0
 integer, intent(out) :: idx(nbf0)
 integer :: n6dmark,n10fmark,n15gmark,n21hmark
 integer :: n5dmark,n7fmark, n9gmark, n11hmark
 ! mark the index where d, f, g, h functions begin
 integer, intent(in) :: shell_type0(ncontr), shell_to_atom_map0(ncontr)
 integer, allocatable :: shell_type(:), shell_to_atom_map(:)
 integer, allocatable :: d_mark(:), f_mark(:), g_mark(:), h_mark(:)
 real(kind=8), intent(out) :: norm(nbf0)

 norm = 1d0   ! initialization
 forall(i = 1:nbf0) idx(i) = i

 k = 2*ncontr
 allocate(shell_type(k), shell_to_atom_map(k))
 shell_type(1:ncontr) = shell_type0
 shell_to_atom_map(1:ncontr) = shell_to_atom_map0
 k = ncontr

! first we adjust the basis functions in each MO according to the Shell to atom map
! this is to ensure that the order of basis functions in PySCF is converted to be that of Gaussian
 ! unsort the shell_type, shell_to_atom_map, MOs will be adjusted accordingly
 call unsort_shell_and_mo(k, shell_type, shell_to_atom_map, nbf0, idx)
 ! Note that k will be updated
 deallocate(shell_to_atom_map)
! adjust done

! then we adjust the basis functions in each MO according to the type of basis functions
 n5dmark  = 0 ; n6dmark  = 0
 n7fmark  = 0 ; n10fmark = 0
 n9gmark  = 0 ; n15gmark = 0
 n11hmark = 0 ; n21hmark = 0
 allocate(d_mark(k), source=0)
 allocate(f_mark(k), source=0)
 allocate(g_mark(k), source=0)
 allocate(h_mark(k), source=0)
 nbf = 0

 do i = 1, k, 1
  select case(shell_type(i))
  case( 0)   !'S'
   nbf = nbf + 1
  case( 1)   !'P'
   nbf = nbf + 3
  case(-1)   !'SP' or 'L'
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
   write(6,'(A)') 'ERROR in subroutine get_permute_idx_from_shell: shell_type(i) out of range!'
   write(6,'(3(A,I0))') 'k=', k, ', i=', i, ', shell_type(i)=', shell_type(i)
   stop
  end select
 end do ! for i
 deallocate(shell_type)

 if(nbf /= nbf0) then
  write(6,'(A)') 'ERROR in subroutine get_permute_idx_from_shell: nbf /= nbf0.'
  write(6,'(2(A,I0))') 'nbf=', nbf, ', nbf0=', nbf0
  stop
 end if

 do i = 1,n5dmark,1
  call dal2fch_permute_5d(idx(d_mark(i):d_mark(i)+4))
 end do
 do i = 1,n6dmark,1
  call dal2fch_permute_6d(idx(d_mark(i):d_mark(i)+5),norm(d_mark(i):d_mark(i)+5))
 end do
 do i = 1,n7fmark,1
  call dal2fch_permute_7f(idx(f_mark(i):f_mark(i)+6))
 end do
 do i = 1,n10fmark,1
  call dal2fch_permute_10f(idx(f_mark(i):f_mark(i)+9),norm(f_mark(i):f_mark(i)+9))
 end do
 do i = 1,n9gmark,1
  call dal2fch_permute_9g(idx(g_mark(i):g_mark(i)+8))
 end do
 do i = 1,n15gmark,1
  call dal2fch_permute_15g(idx(g_mark(i):g_mark(i)+14),norm(g_mark(i):g_mark(i)+14))
 end do
 do i = 1,n11hmark,1
  call dal2fch_permute_11h(idx(h_mark(i):h_mark(i)+10))
 end do
 do i = 1,n21hmark,1
  call dal2fch_permute_21h(idx(h_mark(i):h_mark(i)+20),norm(h_mark(i):h_mark(i)+20))
 end do

 deallocate(d_mark, f_mark, g_mark, h_mark)
end subroutine get_permute_idx_from_shell

! unsort the shell_type, shell_to_atom_map according to the order in Gaussian
! MOs will be adjusted accordingly
subroutine unsort_shell_and_mo(ilen, shell_type, shell_to_atom_map, nbf, idx)
 implicit none
 integer :: i, j, k, length, natom, ibegin, iend, jbegin, jend
 integer, parameter :: ntype = 10
 integer, parameter :: num0(ntype) = [0, 1, -2, 2, -3, 3, -4, 4, -5, 5]
 integer, parameter :: num1(ntype) = [1, 3, 5, 6, 7, 10, 9, 15, 11, 21]
 !                                     S  P  5D 6D 7F 10F 9G 15G 11H 21H
 integer :: num(ntype)
 integer, intent(in) :: nbf
 integer, intent(inout) :: ilen, idx(nbf)
 ! Note that this 'ilen' is not identical to that in fchk2py.f90.
 !  The -1 in array shell_type has not been split. So below we use 2*ilen.
 integer, intent(inout) :: shell_type(2*ilen), shell_to_atom_map(2*ilen)
 integer :: new_shell_type(2*ilen), new_shell_to_atom_map(2*ilen)
 integer, allocatable :: ith(:), new_ith(:), ith_bas(:), tmp_type(:)

 ! split the 'L' into 'S' and 'P'
 new_shell_to_atom_map = shell_to_atom_map
 call split_L_func(ilen, shell_type, new_shell_to_atom_map, length)
 new_shell_type = shell_type

 ! find the number of atoms
 natom = shell_to_atom_map(ilen)

 allocate(ith(0:natom), new_ith(0:natom), ith_bas(0:natom))
 ith = 0
 new_ith = 0
 ith_bas = 0

 ! find the end position of each atom in array shell_to_atom_map
 do i = 1, natom, 1
  ith(i) = count(shell_to_atom_map==i) + ith(i-1)
  new_ith(i) = count(new_shell_to_atom_map==i) + new_ith(i-1)
 end do

 ! find the end position of basis functions of each atom
 do i = 1, natom, 1
  ibegin = new_ith(i-1) + 1
  iend = new_ith(i)
  allocate(tmp_type(ibegin:iend))
  tmp_type = 0
  tmp_type = new_shell_type(ibegin:iend)
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
  ibegin = new_ith(i-1) + 1
  iend = new_ith(i)
  jbegin = ith_bas(i-1) + 1
  jend = ith_bas(i)
  call sort_shell_type_in_each_atom2(iend-ibegin+1, new_shell_type(ibegin:iend))
  call unsort_mo_in_each_atom(iend-ibegin+1, shell_type(ibegin:iend), new_shell_type(ibegin:iend), &
       & jend-jbegin+1, idx(jbegin:jend))
 end do ! for i
 ! adjust the MOs in each atom done

 deallocate(ith, new_ith, ith_bas)
 ilen = length ! update ilen
end subroutine unsort_shell_and_mo

! sort the shell_type within each atom
subroutine sort_shell_type_in_each_atom2(ilen, shell_type)
 implicit none
 integer :: i, tmp_type
 integer, intent(in) :: ilen
 integer, intent(inout) :: shell_type(ilen)
 logical :: sort_done

 sort_done = .false.
 do while(.not. sort_done)
  sort_done = .true.
  do i = 1, ilen-1, 1
    if(shell_type(i) == 0) cycle
    if(ABS(shell_type(i+1)) >= ABS(shell_type(i))) cycle
    sort_done = .false.
    tmp_type = shell_type(i+1)
    shell_type(i+1) = shell_type(i)
    shell_type(i) = tmp_type
  end do
 end do
end subroutine sort_shell_type_in_each_atom2

! Unsort the MO within each atom. (Only for 'L' in Pople basis)
! Example:
!  PySCF:    1s, 2s, 3s, 2px, 2py, 2pz, 3px, 3py, 3pz
!  Gaussian: 1s, 2s, 2px, 2py, 2pz, 3s, 3px, 3py, 3pz
subroutine unsort_mo_in_each_atom(ilen1, shell_type, new_shell_type, ilen2, idx)
 implicit none
 integer :: i, j, k, m   ! temporary variables
 integer :: ibegin
 integer, parameter :: ntype = 10
 integer :: ith_shell(ntype)
 integer, parameter :: num0(ntype) = [0, 1, -2, 2, -3, 3, -4, 4, -5, 5]
 integer, parameter :: num1(ntype) = [1, 3, 5, 6, 7, 10, 9, 15, 11, 21]
 !                                    S  P  5D 6D 7F 10F 9G 15G 11H 21H
 integer, parameter :: rnum(-5:5) = [9, 7, 5, 3, 0, 1, 2, 4, 6, 8, 10]

 integer, intent(in) :: ilen1, ilen2
 integer, intent(in) :: shell_type(ilen1), new_shell_type(ilen1)
 integer, intent(inout) :: idx(ilen2)
 integer, allocatable :: ith_bas(:), new_idx(:)

 ! get the begin position of each type of basis functions within an atom
 ith_shell = 0
 k = rnum(new_shell_type(1))
 ith_shell(k) = 1
 do i = k+1, ntype, 1
  call get_1st_loc(num0(i), ith_shell(i), ilen1, new_shell_type)
 end do ! for i

 ! find the begin position of basis functions within an atom
 allocate(ith_bas(ilen1))
 ith_bas = 0
 ith_bas(1) = 1
 do i = 2, ilen1, 1
  ith_bas(i) = ith_bas(i-1) + num1(rnum(new_shell_type(i-1)))
 end do ! for i

 ! unsort MO indices within an atom
 allocate(new_idx(ilen2), source=0)
 j = 0
 do i = 1, ilen1, 1
  k = rnum(shell_type(i))
  ibegin = ith_bas(ith_shell(k))
  m = num1(k) - 1
  j = j + 1
  new_idx(j:j+m) = idx(ibegin:ibegin+m)
  ith_shell(k) = ith_shell(k) + 1
  j = j + m
 end do
 deallocate(ith_bas)
 ! unsort done

 ! update array idx
 idx = new_idx
 deallocate(new_idx)
end subroutine unsort_mo_in_each_atom

subroutine get_1st_loc(inum, loc, ilen, a)
 implicit none
 integer :: i
 integer, intent(out) :: loc
 integer, intent(in) :: inum, ilen
 integer, intent(in) :: a(ilen)

 do i = 1, ilen, 1
  if(a(i) == inum) exit
 end do
 loc = i
 if(i == ilen+1) loc = 0
end subroutine get_1st_loc

subroutine dal2fch_permute_5d(idx)
 implicit none
 integer :: i, idx0(5)
 integer, parameter :: order(5) = [3, 4, 2, 5, 1]
 integer, intent(inout) :: idx(5)
! From: the order of spherical d functions in PySCF
! To: the order of spherical d functions in Gaussian
! 1    2    3    4    5
! d-2, d-1, d0 , d+1, d+2
! d0 , d+1, d-1, d+2, d-2

 idx0 = idx
 forall(i = 1:5) idx(i) = idx0(order(i))
end subroutine dal2fch_permute_5d

subroutine dal2fch_permute_6d(idx, norm)
 use Sdiag_dalton, only: Sdiag_d
 implicit none
 integer :: i, idx0(6)
 integer, parameter :: order(6) = [1, 4, 6, 2, 3, 5]
 integer, intent(inout) :: idx(6)
 real(kind=8), intent(out) :: norm(6)
! From: the order of Cartesian d functions in PySCF
! To: the order of Cartesian d functions in Gaussian
! 1  2  3  4  5  6
! XX,XY,XZ,YY,YZ,ZZ
! XX,YY,ZZ,XY,XZ,YZ

 idx0 = idx
 forall(i = 1:6)
  idx(i) = idx0(order(i))
  norm(i) = Sdiag_d(order(i))
 end forall
end subroutine dal2fch_permute_6d

subroutine dal2fch_permute_7f(idx)
 implicit none
 integer :: i, idx0(7)
 integer, parameter :: order(7) = [4, 5, 3, 6, 2, 7, 1]
 integer, intent(inout) :: idx(7)
! From: the order of spherical f functions in PySCF
! To: the order of spherical f functions in Gaussian
! 1    2    3    4    5    6    7
! f-3, f-2, f-1, f0 , f+1, f+2, f+3
! f0 , f+1, f-1, f+2, f-2, f+3, f-3

 idx0 = idx
 forall(i = 1:7) idx(i) = idx0(order(i))
end subroutine dal2fch_permute_7f

subroutine dal2fch_permute_10f(idx, norm)
 use Sdiag_dalton, only: Sdiag_f
 implicit none
 integer :: i, idx0(10)
 integer, parameter :: order(10) = [1, 7, 10, 4, 2, 3, 6, 9, 8, 5]
 integer, intent(inout) :: idx(10)
 real(kind=8), intent(out) :: norm(10)
! From: the order of Cartesian f functions in PySCF
! To: the order of Cartesian f functions in Gaussian
! 1   2   3   4   5   6   7   8   9   10
! XXX,XXY,XXZ,XYY,XYZ,XZZ,YYY,YYZ,YZZ,ZZZ
! XXX,YYY,ZZZ,XYY,XXY,XXZ,XZZ,YZZ,YYZ,XYZ

 idx0 = idx
 forall(i = 1:10)
  idx(i) = idx0(order(i))
  norm(i) = Sdiag_f(order(i))
 end forall
end subroutine dal2fch_permute_10f

subroutine dal2fch_permute_9g(idx)
 implicit none
 integer :: i, idx0(9)
 integer, parameter :: order(9) = [5, 6, 4, 7, 3, 8, 2, 9, 1]
 integer, intent(inout) :: idx(9)
! From: the order of spherical g functions in PySCF
! To: the order of spherical g functions in Gaussian
! 1    2    3    4    5    6    7    8    9
! g-4, g-3, g-2, g-1, g0 , g+1, g+2, g+3, g+4
! g0 , g+1, g-1, g+2, g-2, g+3, g-3, g+4, g-4

 idx0 = idx
 forall(i = 1:9) idx(i) = idx0(order(i))
end subroutine dal2fch_permute_9g

subroutine dal2fch_permute_15g(idx, norm)
 use Sdiag_dalton, only: Sdiag_g
 implicit none
 integer :: i, idx0(15)
 integer, intent(inout) :: idx(15)
 real(kind=8), intent(out) :: norm(15)
! From: the order of Cartesian g functions in PySCF
! To: the order of Cartesian g functions in Gaussian
! 1    2    3    4    5    6    7    8    9    10   11   12   13   14   15
! xxxx,xxxy,xxxz,xxyy,xxyz,xxzz,xyyy,xyyz,xyzz,xzzz,yyyy,yyyz,yyzz,yzzz,zzzz
! ZZZZ,YZZZ,YYZZ,YYYZ,YYYY,XZZZ,XYZZ,XYYZ,XYYY,XXZZ,XXYZ,XXYY,XXXZ,XXXY,XXXX

 idx0 = idx
 forall(i = 1:15)
  idx(i) = idx0(16-i)
  norm(i) = Sdiag_g(16-i)
 end forall
end subroutine dal2fch_permute_15g

subroutine dal2fch_permute_11h(idx)
 implicit none
 integer :: i, idx0(11)
 integer, parameter :: order(11) = [6, 7, 5, 8, 4, 9, 3, 10, 2, 11, 1]
 integer, intent(inout) :: idx(11)
! From: the order of Cartesian h functions in PySCF
! To: the order of Cartesian h functions in Gaussian
! 1    2    3    4    5    6    7    8    9    10   11
! h-5, h-4, h-3, h-2, h-1, h0 , h+1, h+2, h+3, h+4, h+5
! h0 , h+1, h-1, h+2, h-2, h+3, h-3, h+4, h-4, h+5, h-5

 idx0 = idx
 forall(i = 1:11) idx(i) = idx0(order(i))
end subroutine dal2fch_permute_11h

subroutine dal2fch_permute_21h(idx, norm)
 use Sdiag_dalton, only: Sdiag_h
 implicit none
 integer :: i, idx0(21)
 integer, intent(inout) :: idx(21)
 real(kind=8), intent(out) :: norm(21)
! From: the order of Cartesian h functions in PySCF
! To: the order of Cartesian h functions in Gaussian
! 1     2     3     4     5     6     7     8     9     10    11    12    13    14    15    16    17    18    19    20    21
! xxxxx,xxxxy,xxxxz,xxxyy,xxxyz,xxxzz,xxyyy,xxyyz,xxyzz,xxzzz,xyyyy,xyyyz,xyyzz,xyzzz,xzzzz,yyyyy,yyyyz,yyyzz,yyzzz,yzzzz,zzzzz
! ZZZZZ,YZZZZ,YYZZZ,YYYZZ,YYYYZ,YYYYY,XZZZZ,XYZZZ,XYYZZ,XYYYZ,XYYYY,XXZZZ,XXYZZ,XXYYZ,XXYYY,XXXZZ,XXXYZ,XXXYY,XXXXZ,XXXXY,XXXXX

 idx0 = idx
 forall(i = 1:21)
  idx(i) = idx0(22-i)
  norm(i) = Sdiag_h(22-i)
 end forall
end subroutine dal2fch_permute_21h

