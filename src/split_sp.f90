! written by jxzou at 20220817: move 'SP'/'L' related subroutines here

! diagonal elements of overlap matrix using Cartesian functions (6D 10F)
! this module is used in fch2py, py2fch and chk2py
module Sdiag_parameter
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
end module Sdiag_parameter

module root_parameter
 implicit none
 real(kind=8), parameter :: root3   = DSQRT(3d0)     ! SQRT(3)
 real(kind=8), parameter :: root5   = DSQRT(5d0)     ! SQRT(5)
 real(kind=8), parameter :: root6   = DSQRT(6d0)     ! SQRT(6)
 real(kind=8), parameter :: root7   = DSQRT(7d0)     ! SQRT(7)
 real(kind=8), parameter :: root9   = 3d0            ! SQRT(9)
 real(kind=8), parameter :: root10  = DSQRT(10d0)    ! SQRT(10)
 real(kind=8), parameter :: root12  = DSQRT(12d0)    ! SQRT(12)
 real(kind=8), parameter :: root15  = DSQRT(15d0)    ! SQRT(15)
 real(kind=8), parameter :: root21  = DSQRT(21d0)    ! SQRT(21)
 real(kind=8), parameter :: root30  = DSQRT(30d0)    ! SQRT(30)
 real(kind=8), parameter :: root35  = DSQRT(35d0)    ! SQRT(35)
 real(kind=8), parameter :: root42  = DSQRT(42d0)    ! SQRT(42)
 real(kind=8), parameter :: root45  = DSQRT(45d0)    ! SQRT(45)
 real(kind=8), parameter :: root105 = DSQRT(105d0)   ! SQRT(105)
 real(kind=8), parameter :: root945 = DSQRT(945d0)   ! SQRT(945)
end module root_parameter

! split the 'L' into 'S' and 'P'
subroutine split_L_func(k, shell_type, shell_to_atom_map, length)
 implicit none
 integer i, k0
 integer,intent(in) :: k
 integer,intent(inout) :: shell_type(2*k), shell_to_atom_map(2*k)
 integer,intent(out) :: length
 integer,allocatable :: temp1(:), temp2(:)

 k0 = 2*k; length = k; i = 1
 ! set initial values for arrays shell_type, assume 15 will not be used
 shell_type(k+1:k0) = 15

 do while(shell_type(i) /= 15)
  if(shell_type(i) /= -1) then
   i = i + 1
   cycle
  end if
  shell_type(i) = 0
  allocate( temp1(i+1 : k0-1), temp2(i+1 : k0-1) )
  temp1(i+1 : k0-1) = shell_type(i+1 : k0-1)
  shell_type(i+2 : k0) = temp1(i+1 : k0-1)
  temp2(i+1 : k0-1) = shell_to_atom_map(i+1 : k0-1)
  shell_to_atom_map(i+2 : k0) = temp2(i+1 : k0-1)
  deallocate(temp1, temp2)
  shell_type(i+1) = 1
  shell_to_atom_map(i+1) = shell_to_atom_map(i)
  i = i + 2
 end do ! for while

 length = i - 1
 shell_type(i : k0) = 0
end subroutine split_L_func

! sort the shell_type, shell_to_atom_map by ascending order
! MOs will be adjusted accordingly
subroutine sort_shell_and_mo(ilen, shell_type, shell_to_atom_map, nbf, nif, coeff2)
 implicit none
 integer :: i, j, k, natom
 integer :: ibegin, iend, jbegin, jend
 integer, parameter :: ntype = 10
 integer, parameter :: num0(ntype) = [0, 1, -2, 2, -3, 3, -4, 4, -5, 5]
 integer, parameter :: num1(ntype) = [1, 3, 5, 6, 7, 10, 9, 15, 11, 21]
 !                                    S  P  5D 6D 7F 10F 9G 15G 11H 21H
 integer :: num(ntype)
 integer,intent(in) :: ilen, nbf, nif
 integer,intent(inout) :: shell_type(ilen), shell_to_atom_map(ilen)
 integer,allocatable :: ith(:), ith_bas(:), tmp_type(:)
 real(kind=8),intent(inout) :: coeff2(nbf,nif)

 ! find the number of atoms
 natom = shell_to_atom_map(ilen)

 allocate(ith(0:natom), ith_bas(0:natom))
 ith = 0
 ith_bas = 0

 ! find the end position of each atom in array shell_to_atom_map
 do i = 1, natom, 1
  ith(i) = count(shell_to_atom_map==i) + ith(i-1)
 end do ! for i

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
  end do ! for j
  k = 0
  do j = 1, ntype, 1
   k = k + num(j)*num1(j)
  end do ! for j
  ith_bas(i) = ith_bas(i-1) + k
  deallocate(tmp_type)
 end do ! for i

 ! adjust the MOs in each atom
 do i = 1, natom, 1
  ibegin = ith(i-1) + 1
  iend = ith(i)
  jbegin = ith_bas(i-1) + 1
  jend = ith_bas(i)
  call sort_shell_and_mo_in_each_atom(iend-ibegin+1, shell_type(ibegin:iend), &
                                      jend-jbegin+1, nif, coeff2(jbegin:jend,1:nif))
 end do ! for i
 ! adjust the MOs in each atom done

 deallocate(ith, ith_bas)
end subroutine sort_shell_and_mo

subroutine sort_shell_and_mo_in_each_atom(ilen1, shell_type, ilen2, nif, coeff2)
 implicit none
 integer :: i, tmp_type
 integer :: ibegin, iend, jbegin, jend
 integer,parameter :: ntype = 10
 integer,parameter :: num1(ntype) = [1, 3, 5, 6, 7, 10, 9, 15, 11, 21]
 !                                   S  P  5D 6D 7F 10F 9G 15G 11H 21H
 integer,parameter :: rnum(-5:5) = [9, 7, 5, 3, 0, 1, 2, 4, 6, 8, 10]

 integer,intent(in) :: nif, ilen1, ilen2
 integer,intent(inout) :: shell_type(ilen1)
 integer,allocatable :: ith_bas(:)
 real(kind=8),intent(inout) :: coeff2(ilen2,nif)
 real(kind=8),allocatable :: tmp_coeff1(:,:), tmp_coeff2(:,:)
 logical :: sort_done

 tmp_type = 0

 ! find the end position of basis functions within an atom
 allocate(ith_bas(0:ilen1))
 ith_bas = 0
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
    ibegin = ith_bas(i-1) + 1
    iend = ith_bas(i)
    jbegin = ith_bas(i) + 1
    jend = ith_bas(i+1)
    allocate(tmp_coeff1(ibegin:iend,nif), tmp_coeff2(jbegin:jend,nif))
    tmp_coeff1 = 0d0
    tmp_coeff2 = 0d0
    tmp_coeff1 = coeff2(ibegin:iend,:)
    tmp_coeff2 = coeff2(jbegin:jend,:)
    ith_bas(i) = ibegin+jend-jbegin
    coeff2(ibegin: ith_bas(i),:) = tmp_coeff2
    ith_bas(i+1) = jend+iend-jbegin+1
    coeff2(ibegin+jend-jbegin+1: ith_bas(i+1),:) = tmp_coeff1
    deallocate(tmp_coeff1, tmp_coeff2)
  end do ! for i
 end do ! for while

 deallocate(ith_bas)
end subroutine sort_shell_and_mo_in_each_atom

! sort the shell_type, shell2atom_map by ascending order,
! indices of MOs will be adjusted accordingly
subroutine sort_shell_and_mo_idx(ilen, shell_type, shell2atom_map, nbf, idx)
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
  call sort_shell_and_mo_in_each_atom_idx(iend-ibegin+1, shell_type(ibegin:iend),&
                                          jend-jbegin+1, idx(jbegin:jend))
 end do ! for i

 deallocate(ith, ith_bas)
end subroutine sort_shell_and_mo_idx

subroutine sort_shell_and_mo_in_each_atom_idx(ilen1, shell_type, ilen2, idx)
 implicit none
 integer :: i, tmp_type, ibegin, iend, jbegin, jend
 integer, parameter :: ntype = 10
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
end subroutine sort_shell_and_mo_in_each_atom_idx

subroutine fch2inporb_permute_5d(idx)
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
end subroutine fch2inporb_permute_5d

subroutine fch2inporb_permute_6d(idx, norm)
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

 idx0 = idx
 norm0(1:3) = norm(1:3)/root3
 norm0(4:6) = norm(4:6)
 forall(i = 1:6)
  idx(i) = idx0(order(i))
  norm(i) = norm0(order(i))
 end forall
end subroutine fch2inporb_permute_6d

subroutine fch2inporb_permute_7f(idx)
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
end subroutine fch2inporb_permute_7f

subroutine fch2inporb_permute_10f(idx, norm)
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

 idx0 = idx
 norm0(1:3) = norm(1:3)/root15
 norm0(4:9) = norm(4:9)/root3
 norm0(10) = norm(10)
 forall(i = 1:10)
  idx(i) = idx0(order(i))
  norm(i) = norm0(order(i))
 end forall
end subroutine fch2inporb_permute_10f

subroutine fch2inporb_permute_9g(idx)
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
end subroutine fch2inporb_permute_9g

subroutine fch2inporb_permute_15g(idx, norm)
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

 idx0 = idx
 norm0 = norm/ratio
 forall(i = 1:15)
  idx(i) = idx0(16-i)
  norm(i) = norm0(16-i)
 end forall
end subroutine fch2inporb_permute_15g

subroutine fch2inporb_permute_11h(idx)
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
end subroutine fch2inporb_permute_11h

subroutine fch2inporb_permute_21h(idx, norm)
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

 idx0 = idx
 norm0 = norm/ratio
 forall(i = 1:21)
  idx(i) = idx0(22-i)
  norm(i) = norm0(22-i)
 end forall
end subroutine fch2inporb_permute_21h

subroutine fch2inporb_permute_sph(n5dmark, n7fmark, n9gmark, n11hmark, k, &
                                  d_mark, f_mark, g_mark, h_mark, nbf, idx)
 implicit none
 integer :: i, j
 integer, intent(in) :: n5dmark, n7fmark, n9gmark, n11hmark, k, nbf
 integer, intent(in) :: d_mark(k), f_mark(k), g_mark(k), h_mark(k)
 integer, intent(inout) :: idx(nbf)

 do i = 1, n5dmark, 1
  j = d_mark(i)
  call fch2inporb_permute_5d(idx(j:j+4))
 end do ! for i

 do i = 1, n7fmark, 1
  j = f_mark(i)
  call fch2inporb_permute_7f(idx(j:j+6))
 end do ! for i

 do i = 1, n9gmark, 1
  j = g_mark(i)
  call fch2inporb_permute_9g(idx(j:j+8))
 end do ! for i

 do i = 1, n11hmark, 1
  j = h_mark(i)
  call fch2inporb_permute_11h(idx(j:j+10))
 end do ! for i
end subroutine fch2inporb_permute_sph

subroutine fch2inporb_permute_cart(n6dmark, n10fmark, n15gmark, n21hmark, k, &
                               d_mark, f_mark, g_mark, h_mark, nbf, idx, norm)
 implicit none
 integer :: i, j
 integer, intent(in) :: n6dmark, n10fmark, n15gmark, n21hmark, k, nbf
 integer, intent(in) :: d_mark(k), f_mark(k), g_mark(k), h_mark(k)
 integer, intent(inout) :: idx(nbf)
 real(kind=8), intent(inout) :: norm(nbf)

 do i = 1, n6dmark, 1
  j = d_mark(i)
  call fch2inporb_permute_6d(idx(j:j+5), norm(j:j+5))
 end do
 do i = 1, n10fmark, 1
  j = f_mark(i)
  call fch2inporb_permute_10f(idx(j:j+9), norm(j:j+9))
 end do
 do i = 1, n15gmark, 1
  j = g_mark(i)
  call fch2inporb_permute_15g(idx(j:j+14), norm(j:j+14))
 end do
 do i = 1, n21hmark, 1
  j = h_mark(i)
  call fch2inporb_permute_21h(idx(j:j+20), norm(j:j+20))
 end do
end subroutine fch2inporb_permute_cart

! move the 2nd, 3rd, ... Zeta basis functions forward
subroutine zeta_mv_forwd_idx(i0, shell_type, length, nbf, idx2, norm1)
 implicit none
 integer :: i, j, k
 integer, intent(in) :: i0, shell_type, length, nbf
 integer, parameter :: num0(-5:5) = [11, 9, 7, 5, 0, 0, 3, 6, 10, 15, 21]
 !                                   11H 9G 7F 5D L  S 3P 6D 10F 15G 21H
 integer, intent(inout) :: idx2(nbf)
 integer, allocatable :: idx(:)
 real(kind=8), allocatable :: norm(:)
 real(kind=8), intent(inout) :: norm1(nbf)

 if(length == 1) return
 if(shell_type==0 .or. shell_type==-1) then
  write(6,'(/,A)') 'ERROR in subroutine zeta_mv_forwd_idx: this element of&
                   & shell_type is'
  write(6,'(A)') '0 or -1. Impossible.'
  write(6,'(2(A,I0))') 'shell_type=', shell_type, ', length=', length
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
end subroutine zeta_mv_forwd_idx

! read the location indices of 3p from the array shell_type
! Note: can be used for spherical harmonic type or Cartesian type basis
subroutine read3pmark_from_shltyp(ncontr, shltyp, np, p_mark)
 implicit none
 integer :: i, k, nbf
 integer, intent(in) :: ncontr
 integer, intent(in) :: shltyp(ncontr)
 integer, intent(out) :: np, p_mark(ncontr)
 integer, parameter :: nbas_per_ang(-6:6) = [13,11,9,7,5,4,1,3,6,10,15,21,28]

 if(ANY(shltyp == -1)) then
  write(6,'(/,A)') 'ERROR in subroutine read3pmark_from_shltyp: there exists so&
                   &me element -1 in'
  write(6,'(A)') 'the array shltyp. The subroutine split_L_func is supposed to &
                 &be called before'
  write(6,'(A)') 'calling this subroutine. shltyp='
  write(6,'(20I5)') shltyp
  stop
 end if

 nbf = 0; np = 0; p_mark = 0

 do i = 1, ncontr, 1
  k = shltyp(i)
  if(k == 1) then
   np = np + 1
   p_mark(np) = nbf + 1
  end if
  nbf = nbf + nbas_per_ang(k)
 end do ! for i
end subroutine read3pmark_from_shltyp

! read the position marks of 5D, 7F, etc from array shell_type
! Note: only used for spherical harmonic type basis
subroutine read_mark_from_shltyp_sph(ncontr, shltyp, nd, nf, ng, nh, d_mark, &
                                     f_mark, g_mark, h_mark)
 implicit none
 integer :: i, k, nbf
 integer, intent(in) :: ncontr
 integer, intent(in) :: shltyp(ncontr)
 integer, intent(out) :: nd, nf, ng, nh, d_mark(ncontr), f_mark(ncontr), &
  g_mark(ncontr), h_mark(ncontr)
 integer, parameter :: nbas_per_ang(-6:6) = [13,11,9,7,5,4,1,3,6,10,15,21,28]

 nbf = 0; nd = 0; nf = 0; ng = 0; nh = 0
 d_mark = 0; f_mark = 0; g_mark = 0; h_mark = 0
!     Spherical      |     Cartesian
! -6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6
!  I  H  G  F  D  L  S  P  D  F  G  H  I

 do i = 1, ncontr, 1
  k = shltyp(i)
  select case(k)
  case(-1,0,1)
  case(-2) ! 5D
   nd = nd + 1
   d_mark(nd) = nbf + 1
  case(-3) ! 7F
   nf = nf + 1
   f_mark(nf) = nbf + 1
  case(-4) ! 9G
   ng = ng + 1
   g_mark(ng) = nbf + 1
  case(-5) ! 11H
   nh = nh + 1
   h_mark(nh) = nbf + 1
  case default
   write(6,'(/,A)') 'ERROR in subroutine read_mark_from_shltyp_sph:'
   write(6,'(A,I0)') 'Invalid shltyp(i)=', shltyp(i)
   stop
  end select
  nbf = nbf + nbas_per_ang(k)
 end do ! for i
end subroutine read_mark_from_shltyp_sph

! read the position marks of 6D, 10F, etc from array shell_type
! Note: only used for Cartesian type basis
subroutine read_mark_from_shltyp_cart(ncontr, shltyp, nd, nf, ng, nh, d_mark, &
                                      f_mark, g_mark, h_mark)
 implicit none
 integer :: i, k, nbf
 integer, intent(in) :: ncontr
 integer, intent(in) :: shltyp(ncontr)
 integer, intent(out) :: nd, nf, ng, nh, d_mark(ncontr), f_mark(ncontr), &
  g_mark(ncontr), h_mark(ncontr)
 integer, parameter :: nbas_per_ang(-6:6) = [13,11,9,7,5,4,1,3,6,10,15,21,28]

 nbf = 0; nd = 0; nf = 0; ng = 0; nh = 0
 d_mark = 0; f_mark = 0; g_mark = 0; h_mark = 0
!     Spherical      |     Cartesian
! -6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6
!  I  H  G  F  D  L  S  P  D  F  G  H  I

 do i = 1, ncontr, 1
  k = shltyp(i)
  select case(k)
  case(-1,0,1)
  case(2) ! 6D
   nd = nd + 1
   d_mark(nd) = nbf + 1
  case(3) ! 10F
   nf = nf + 1
   f_mark(nf) = nbf + 1
  case(4) ! 15G
   ng = ng + 1
   g_mark(ng) = nbf + 1
  case(5) ! 21H
   nh = nh + 1
   h_mark(nh) = nbf + 1
  case default
   write(6,'(/,A)') 'ERROR in subroutine read_mark_from_shltyp_cart:'
   write(6,'(A,I0)') 'Invalid shltyp(i)=', shltyp(i)
   stop
  end select
  nbf = nbf + nbas_per_ang(k)
 end do ! for i
end subroutine read_mark_from_shltyp_cart

! get Gaussian -> Q-Chem permutation indices
subroutine get_fch2qchem_permute_idx(sph, ncontr, shell_type, nbf, idx)
 implicit none
 integer :: i
 integer :: n5dmark, n7fmark, n9gmark, n11hmark
 integer :: n6dmark, n10fmark, n15gmark, n21hmark
 integer, intent(in) :: ncontr, nbf
 integer, intent(in) :: shell_type(ncontr)
 integer, intent(out) :: idx(nbf)
 integer, allocatable :: d_mark(:), f_mark(:), g_mark(:), h_mark(:)
 logical, intent(in) :: sph

 allocate(d_mark(ncontr), f_mark(ncontr), g_mark(ncontr), h_mark(ncontr))
 forall(i = 1:nbf) idx(i) = i

 if(sph) then
  call read_mark_from_shltyp_sph(ncontr, shell_type, n5dmark, n7fmark, n9gmark, &
                                 n11hmark, d_mark, f_mark, g_mark, h_mark)
  ! adjust the order of 5d, 7f functions, etc
  call fch2inporb_permute_sph(n5dmark, n7fmark, n9gmark, n11hmark, ncontr, &
                              d_mark, f_mark, g_mark, h_mark, nbf, idx)
 else
  call read_mark_from_shltyp_cart(ncontr, shell_type, n6dmark, n10fmark, n15gmark,&
                                  n21hmark, d_mark, f_mark, g_mark, h_mark)
  ! adjust the order of 6d, 10f functions, etc
  call fch2qchem_permute_cart(n6dmark, n10fmark, n15gmark, n21hmark, ncontr, &
                              d_mark, f_mark, g_mark, h_mark, nbf, idx)
 end if

 deallocate(d_mark, f_mark, g_mark, h_mark)
end subroutine get_fch2qchem_permute_idx

subroutine fch2qchem_permute_cart(n6dmark, n10fmark, n15gmark, n21hmark, k, d_mark, &
                                  f_mark, g_mark, h_mark, nbf, idx)
 implicit none
 integer :: i
 integer, intent(in) :: n6dmark, n10fmark, n15gmark, n21hmark, k, nbf
 integer, intent(in) :: d_mark(k), f_mark(k), g_mark(k), h_mark(k)
 integer, intent(inout) :: idx(nbf)

 do i = 1, n6dmark, 1
  call fch2qchem_permute_6d(idx(d_mark(i):d_mark(i)+5))
 end do
 do i = 1, n10fmark, 1
  call fch2qchem_permute_10f(idx(f_mark(i):f_mark(i)+9))
 end do
 do i = 1, n15gmark, 1
  call fch2qchem_permute_15g(idx(g_mark(i):g_mark(i)+14))
 end do
 do i = 1, n21hmark, 1
  call fch2qchem_permute_21h(idx(h_mark(i):h_mark(i)+20))
 end do
end subroutine fch2qchem_permute_cart

subroutine fch2qchem_permute_6d(idx)
 implicit none
 integer :: i, idx0(6)
 integer, parameter :: order(6) = [1, 4, 2, 5, 6, 3]
 integer, intent(inout) :: idx(6)
! From: the order of Cartesian d functions in Gaussian
! To: the order of Cartesian d functions in QChem
! 1  2  3  4  5  6
! XX,YY,ZZ,XY,XZ,YZ
! XX,XY,YY,XZ,YZ,ZZ

 idx0 = idx
 forall(i = 1:6) idx(i) = idx0(order(i))
end subroutine fch2qchem_permute_6d

subroutine fch2qchem_permute_10f(idx)
 implicit none
 integer :: i, idx0(10)
 integer, parameter :: order(10) = [1, 5, 4, 2, 6, 10, 9, 7, 8, 3]
 integer, intent(inout) :: idx(10)
! From: the order of Cartesian f functions in Gaussian
! To: the order of Cartesian f functions in Qchem
! 1   2   3   4   5   6   7   8   9   10
! XXX,YYY,ZZZ,XYY,XXY,XXZ,XZZ,YZZ,YYZ,XYZ
! XXX,XXY,XYY,YYY,XXZ,XYZ,YYZ,XZZ,YZZ,ZZZ

 idx0 = idx
 forall(i = 1:10) idx(i) = idx0(order(i))
end subroutine fch2qchem_permute_10f

subroutine fch2qchem_permute_15g(idx)
 implicit none
 integer :: i, idx0(15)
 integer, parameter :: order(15) = [15, 14, 12, 9, 5, 13, 11, 8, 4, 10, 7, 3, 6, 2, 1]
 integer, intent(inout) :: idx(15)
! From: the order of Cartesian g functions in Gaussian
! To: the order of Cartesian g functions in QChem
! 1    2    3    4    5    6    7    8    9    10   11   12   13   14   15
! ZZZZ,YZZZ,YYZZ,YYYZ,YYYY,XZZZ,XYZZ,XYYZ,XYYY,XXZZ,XXYZ,XXYY,XXXZ,XXXY,XXXX
! xxxx,xxxy,xxyy,xyyy,yyyy,xxxz,xxyz,xyyz,yyyz,xxzz,xyzz,yyzz,xzzz,yzzz,zzzz

 idx0 = idx
 forall(i = 1:15) idx(i) = idx0(order(i))
end subroutine fch2qchem_permute_15g

subroutine fch2qchem_permute_21h(idx)
 implicit none
 integer :: i, idx0(21)
 integer, parameter :: order(21) = [21, 20, 18, 15, 11, 6, 19, 17, 14, 10, 5, &
                                    16, 13, 9, 4, 12, 8, 3, 7, 2, 1]
 integer, intent(inout) :: idx(21)
! From: the order of Cartesian h functions in Gaussian
! To: the order of Cartesian h functions in QChem
! 1     2     3     4     5     6     7     8     9     10    11    12    13    14    15    16    17    18    19    20    21
! ZZZZZ,YZZZZ,YYZZZ,YYYZZ,YYYYZ,YYYYY,XZZZZ,XYZZZ,XYYZZ,XYYYZ,XYYYY,XXZZZ,XXYZZ,XXYYZ,XXYYY,XXXZZ,XXXYZ,XXXYY,XXXXZ,XXXXY,XXXXX
! xxxxx,xxxxy,xxxyy,xxyyy,xyyyy,yyyyy,xxxxz,xxxyz,xxyyz,xyyyz,yyyyz,xxxzz,xxyzz,xyyzz,yyyzz,xxzzz,xyzzz,yyzzz,xzzzz,yzzzz,zzzzz

 idx0 = idx
 forall(i = 1:21) idx(i) = idx0(order(i))
end subroutine fch2qchem_permute_21h

subroutine fch2tm_permute_sph(n5dmark, n7fmark, n9gmark, n11hmark, k, d_mark, &
                              f_mark, g_mark, h_mark, nbf, idx, norm)
 implicit none
 integer :: i, j, m
 integer, intent(in) :: n5dmark, n7fmark, n9gmark, n11hmark, k, nbf
 integer, intent(in) :: d_mark(k), f_mark(k), g_mark(k), h_mark(k)
 integer, intent(inout) :: idx(nbf)
 real(kind=8), intent(inout) :: norm(nbf)

 do i = 1, n5dmark, 1
  j = d_mark(i) + 3
  m = idx(j)
  idx(j) = idx(j+1)
  idx(j+1) = m
 end do ! for i

 do i = 1, n7fmark, 1
  j = f_mark(i) + 3
  m = idx(j)
  idx(j) = idx(j+1)
  idx(j+1) = m
  norm(j+3) = -norm(j+3)
 end do ! for i

 do i = 1, n9gmark, 1
  j = g_mark(i) + 3
  m = idx(j)
  idx(j) = idx(j+1)
  idx(j+1) = m
  norm(j+1) = -norm(j+1)
  norm(j+3) = -norm(j+3)
  m = idx(j+4)
  idx(j+4) = idx(j+5)
  idx(j+5) = m
 end do ! for i

 do i = 1, n11hmark, 1
  j = h_mark(i) + 3
  m = idx(j)
  idx(j) = idx(j+1)
  idx(j+1) = m
  m = idx(j+4)
  idx(j+4) = idx(j+5)
  idx(j+5) = m
 end do ! for i
end subroutine fch2tm_permute_sph

subroutine fch2tm_permute_cart(n6dmark, n10fmark, n15gmark, n21hmark, k, d_mark,&
                               f_mark, g_mark, h_mark, nbf, idx, norm)
 implicit none
 integer :: i, j
 integer, intent(in) :: n6dmark, n10fmark, n15gmark, n21hmark, k, nbf
 integer, intent(in) :: d_mark(k), f_mark(k), g_mark(k), h_mark(k)
 integer, intent(inout) :: idx(nbf)
 real(kind=8), intent(inout) :: norm(nbf)

 do i = 1, n6dmark, 1
  j = d_mark(i)
  call fch2tm_permute_6d(idx(j:j+5), norm(j:j+5))
 end do
 do i = 1, n10fmark, 1
  j = f_mark(i)
  call fch2tm_permute_10f(idx(j:j+9), norm(j:j+9))
 end do
 do i = 1, n15gmark, 1
  j = g_mark(i)
  call fch2tm_permute_15g(idx(j:j+14), norm(j:j+14))
 end do
 do i = 1, n21hmark, 1
  j = h_mark(i)
  call fch2tm_permute_21h(idx(j:j+20), norm(j:j+20))
 end do
end subroutine fch2tm_permute_cart

subroutine fch2tm_permute_6d(idx, norm)
 implicit none
 integer :: i, idx0(6)
 integer, parameter :: order(6) = [1, 5, 6, 4, 2, 3]
 integer, intent(inout) :: idx(6)
 real(kind=8) :: norm0(6)
 real(kind=8), intent(inout) :: norm(6)

 idx0 = idx
 norm0 = norm
 forall(i = 1:6)
  idx(i) = idx0(order(i))
  norm(i) = norm0(order(i))
 end forall
end subroutine fch2tm_permute_6d

subroutine fch2tm_permute_10f(idx, norm)
 implicit none
 integer :: i, idx0(10)
 integer, parameter :: order(10) = [1,2,3,4,5,6,7,8,9,10]
 integer, intent(inout) :: idx(10)
 real(kind=8) :: norm0(10)
 real(kind=8), intent(inout) :: norm(10)

 idx0 = idx
 norm0 = norm
 forall(i = 1:10)
  idx(i) = idx0(order(i))
  norm(i) = norm0(order(i))
 end forall
end subroutine fch2tm_permute_10f

subroutine fch2tm_permute_15g(idx, norm)
 implicit none
 integer :: i
 integer, intent(inout) :: idx(15)
 real(kind=8), intent(inout) :: norm(15)

 forall(i = 1:15) idx(i) = i
 norm = 1d0
 write(6,'(/,A)') 'ERROR in subroutine fch2tm_permute_15g: Cartesian-type g fun&
                  &ctions not'
 write(6,'(A)') 'supported currently.'
 stop
end subroutine fch2tm_permute_15g

subroutine fch2tm_permute_21h(idx, norm)
 implicit none
 integer :: i
 integer, intent(inout) :: idx(21)
 real(kind=8), intent(inout) :: norm(21)

 forall(i = 1:21) idx(i) = i
 norm = 1d0
 write(6,'(/,A)') 'ERROR in subroutine fch2tm_permute_21h: Cartesian-type h fun&
                  &ctions not'
 write(6,'(A)') 'supported currently.'
 stop
end subroutine fch2tm_permute_21h

subroutine fch2mrcc_permute_6d(idx, norm)
 use root_parameter, only: root3
 implicit none
 integer :: i, idx0(6)
 integer, parameter :: order(6) = [3,6,2,5,4,1]
 integer, intent(inout) :: idx(6)
 real(kind=8), parameter :: ratio(6) = [1d0,root3,1d0,root3,root3,1d0]
 real(kind=8), intent(inout) :: norm(6)

 idx0 = idx
 do i = 1, 6
  idx(i) = idx0(order(i))
  norm(i) = norm(i)*ratio(i)
 end do ! for i
end subroutine fch2mrcc_permute_6d

subroutine fch2mrcc_permute_10f(idx, norm)
 use root_parameter, only: root5, root15
 implicit none
 integer :: i, idx0(10)
 integer, parameter :: order(10) = [3,8,9,2,7,10,4,6,5,1]
 integer, intent(inout) :: idx(10)
 real(kind=8), intent(inout) :: norm(10)
 real(kind=8), parameter :: ratio(10) = [1d0,root5,root5,1d0,root5,root15,&
  root5,root5,root5,1d0]

 idx0 = idx
 do i = 1, 10
  idx(i) = idx0(order(i))
  norm(i) = norm(i)*ratio(i)
 end do ! for i
end subroutine fch2mrcc_permute_10f

subroutine fch2mrcc_permute_15g(norm)
 use root_parameter, only: root7, root35, root105
 implicit none
 integer :: i
 real(kind=8), intent(inout) :: norm(15)
 real(kind=8), parameter :: r105d3 = root105/3d0
 real(kind=8), parameter :: ratio(15) = [1d0,root7,r105d3,root7,1d0,root7,&
  root35,root35,root7,r105d3,root35,r105d3,root7,root7,1d0]

 do i = 1, 15
  norm(i) = norm(i)*ratio(i)
 end do ! for i
end subroutine fch2mrcc_permute_15g

subroutine fch2mrcc_permute_21h(norm)
 use root_parameter, only: root7, root21, root105
 implicit none
 integer :: i
 real(kind=8), intent(inout) :: norm(21)
 real(kind=8), parameter :: p3r7 = 3d0*root7
 real(kind=8), parameter :: ratio(21) = [1d0,3d0,root21,root21,3d0,1d0,3d0,&
  p3r7,root105,p3r7,3d0,root21,root105,root105,root21,root21,p3r7,root21,3d0,&
  3d0,1d0]

 do i = 1, 21
  norm(i) = norm(i)*ratio(i)
 end do ! for i
end subroutine fch2mrcc_permute_21h

subroutine fch2mrcc_permute_cart(n6dmark, n10fmark, n15gmark, n21hmark, k, &
                             d_mark, f_mark, g_mark, h_mark, nbf, idx, norm)
 implicit none
 integer :: i, j
 integer, intent(in) :: n6dmark, n10fmark, n15gmark, n21hmark, k, nbf
 integer, intent(in) :: d_mark(k), f_mark(k), g_mark(k), h_mark(k)
 integer, intent(inout) :: idx(nbf)
 real(kind=8), intent(inout) :: norm(nbf)

 do i = 1, n6dmark, 1
  j = d_mark(i)
  call fch2mrcc_permute_6d(idx(j:j+5), norm(j:j+5))
 end do ! for i

 do i = 1, n10fmark, 1
  j = f_mark(i)
  call fch2mrcc_permute_10f(idx(j:j+9), norm(j:j+9))
 end do ! for i

 do i = 1, n15gmark, 1
  j = g_mark(i)
  call fch2mrcc_permute_15g(norm(j:j+14))
 end do ! for i

 do i = 1, n21hmark, 1
  j = h_mark(i)
  call fch2mrcc_permute_21h(norm(j:j+20))
 end do ! for i
end subroutine fch2mrcc_permute_cart

