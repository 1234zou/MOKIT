! written by jxzou at 20171203: adjust the orders of d,f,g, etc functions in
!  .fch(k) file, to the orders in PySCF

! updated by jxzou at 20180314: add a input parameter 'a' or 'b' to read Alpha or Beta MO in .fchk
! updated by jxzou at 20180323: support the case that D functions preceding L functions
! updated by jxzou at 20180406: code optimization
! updated by jxzou at 20180520: support linear dependence
! updated by jxzou at 20190226: change data to parameter when declare constant arrays
! updated by jxzou at 20190411: pass the Sdiag in
! updated by jxzou at 20200324: renamed from fch2py as fch2py, simplify code
! updated by jxzou at 20210527: remove intent(in) parameter Sdiag, use parameter array
! updated by jxzou at 20220810: support Gau->PySCF complex GHF

! diagonal elements of overlap matrix using Cartesian functions (6D 10F)
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

! read the MOs in .fch(k) file and adjust its d,f,g,h functions order
!  of Gaussian to that of PySCF
subroutine fch2py(fchname, nbf, nif, ab, coeff2)
 implicit none
 integer :: i, j, k, length, nbf0, ncoeff, fchid
 integer, intent(in) :: nbf, nif
!f2py intent(in) :: nbf, nif

 integer :: n6dmark,n10fmark,n15gmark,n21hmark
 integer :: n5dmark,n7fmark, n9gmark, n11hmark
 integer, allocatable :: shell_type(:), shell_to_atom_map(:)
 ! mark the index where d, f, g, h functions begin
 integer, allocatable :: d_mark(:), f_mark(:), g_mark(:), h_mark(:)

 real(kind=8), intent(out) :: coeff2(nbf,nif)
!f2py depend(nbf,nif) :: coeff2
!f2py intent(out) :: coeff2
 real(kind=8), allocatable :: coeff(:)

 character(len=1), intent(in) :: ab
!f2py intent(in) :: ab

 character(len=8) :: key
 character(len=8), parameter :: key1 = 'Alpha MO'
 character(len=7), parameter :: key2 = 'Beta MO'
 character(len=240) :: buffer
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 logical :: alive, ghf

 key = ' '; buffer = ' '
 ncoeff = 0; ghf = .false.

 select case(ab)
 case('a')
  key = key1
 case('b')
  key = key2//' '
 case('r','i') ! r/i for real/imaginary in complex GHF
  key = key1
  ghf = .true.
 case default
  write(6,'(/,A)') 'ERROR in subroutine fch2py: wrong data type of ab.'
  write(6,'(A)') "This argument can only be 'a'/'b'/'r'/'i'. But your input"
  write(6,*) 'ab=', ab
  stop
 end select

 inquire(file=TRIM(fchname),exist=alive)
 if(.not. alive) then
  write(6,'(/,A)') "File '"//TRIM(fchname)//"' does not exist!"
  stop
 end if

 open(newunit=fchid,file=TRIM(fchname),status='old',position='rewind')
 do while(.true.)
  read(fchid,'(A)') buffer
  if(buffer(1:8) == key) exit
 end do
 BACKSPACE(fchid)
 read(fchid,'(A49,2X,I10)') buffer, ncoeff

 if(ghf) then
  if(ncoeff /= 2*nbf*nif) then
   close(fchid)
   write(6,'(/,A)') 'ERROR in subroutine fch2py: ncoeff/=2*nbf*nif! Inconsisten&
                    &t basis sets'
   write(6,'(A)') 'in PySCF script and file '//TRIM(fchname)
   write(6,'(3(A,I0))') 'ncoeff=', ncoeff, ', nbf=', nbf, ', nif=', nif
   stop
  end if
 else ! R(O)HF, UHF
  if(ncoeff /= nbf*nif) then
   close(fchid)
   write(6,'(/,A)') 'ERROR in subroutine fch2py: ncoeff/=nbf*nif! Inconsistent &
                    &basis sets in'
   write(6,'(A)') 'PySCF script and file '//TRIM(fchname)
   write(6,'(3(A,I0))') 'ncoeff=', ncoeff, ', nbf=', nbf, ', nif=', nif
   stop
  end if
 end if

 ! read Alpha MO or Beta MO
 allocate(coeff(ncoeff), source=0d0)
 read(fchid,'(5(1X,ES15.8))') (coeff(i),i=1,ncoeff)

 ! For complex GHF, the MO coefficients of Gaussian are in the order
 !  (C_ar, C_ai, C_br, C_bi)X1 + ()X2 + ... which can be seen in the log file
 !  by specifying the keyword 'pop=full'
 ! While for PySCF, the real part of MO coefficients are in the order
 !  C_ar_1, C_ar_2, ..., C_ar_nbf0, C_br_1, C_br_2, ..., C_br_nbf0
 !  the imaginary part of MO coefficients are in the order
 !  C_ai_1, C_ai_2, ..., C_ai_nbf0, C_bi_1, C_bi_2, ..., C_bi_nbf0

 if(ghf) then ! complex GHF
  select case(ab)
  case('r') ! real part
   k = 0
   do i = 1, nif, 1
    j = k + 2*nbf
    coeff2(1:nbf/2,i) = coeff(k+1:j-3:4)
    coeff2(nbf/2+1:nbf,i) = coeff(k+3:j-1:4)
    k = j ! update k
   end do ! for i
  case('i') ! imarginary part
   k = 0
   do i = 1, nif, 1
    j = k + 2*nbf
    coeff2(1:nbf/2,i) = coeff(k+2:j-2:4)
    coeff2(nbf/2+1:nbf,i) = coeff(k+4:j:4)
    k = j ! update k
   end do ! for i
  case default
   write(6,'(/,A)') 'ERROR in subroutine fch2py: wrong parameter ab.'
   write(6,*) 'ab=', ab
   stop
  end select
 else         ! real R(O)HF, UHF
  coeff2 = RESHAPE(coeff,(/nbf,nif/))
 end if

 deallocate(coeff)

 rewind(fchid)   ! find and read Shell types
 do while(.true.)
  read(fchid,'(A)',iostat=i) buffer
  if(buffer(1:11) == 'Shell types') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine fch2py: missing the 'Shell types'&
                & section in .fch file!"
  write(6,'(A)') TRIM(fchname)
  close(fchid)
  stop
 end if

 BACKSPACE(fchid)
 read(fchid,'(A49,2X,I10)') buffer, k

 allocate(shell_type(2*k), source=0)
 read(fchid,'(6(6X,I6))') (shell_type(i),i=1,k)
 ! read Shell types done

 ! find and read Shell to atom map
 do while(.true.)
  read(fchid,'(A)',iostat=i) buffer
  if(buffer(1:13) == 'Shell to atom') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine fch2py: missing the 'Shell to atom map'&
                & section in .fch file!"
  write(6,'(A)') TRIM(fchname)
  close(fchid)
  stop
 end if

 allocate(shell_to_atom_map(2*k), source=0)
 read(fchid,'(6(6X,I6))') (shell_to_atom_map(i),i=1,k)
 ! read Shell to atom map done

 ! all information in .fchk file read done
 close(fchid)

! first we adjust the basis functions in each MO according to the Shell to atom map
! this is to ensure that D comes after L functions
 ! split the 'L' into 'S' and 'P'
 call split_L_func(k, shell_type, shell_to_atom_map, length)

 ! Sort the shell_type, shell_to_atom_map by ascending order. MOs will be
 ! adjusted accordingly
 if(ghf) then ! GHF
  ! Note: arrays shell_type and shell_to_atom_map will be changed after calling
  !  subroutine sort_shell_and_mo, so we need to backup these two arrays.
  !  Using d_mark and f_mark as intermediate arrays here

  allocate(d_mark(2*k), source=shell_type)
  allocate(f_mark(2*k), source=shell_to_atom_map)
  call sort_shell_and_mo(length, shell_type, shell_to_atom_map, nbf/2, nif, &
                         coeff2(1:nbf/2,:))
  ! restore shell_type and shell_to_atom_map from d_mark and f_mark
  shell_type = d_mark
  shell_to_atom_map = f_mark
  deallocate(d_mark, f_mark)
  call sort_shell_and_mo(length, shell_type, shell_to_atom_map, nbf/2, nif, &
                         coeff2(nbf/2+1:nbf,:))
 else         ! R(O)HF, UHF
  call sort_shell_and_mo(length, shell_type, shell_to_atom_map, nbf, nif, coeff2)
 end if

 deallocate(shell_to_atom_map)
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
 nbf0 = 0

 do i = 1, k, 1
  select case(shell_type(i))
  case( 0)   ! S
   nbf0 = nbf0 + 1
  case( 1)   ! 3P
   nbf0 = nbf0 + 3
  case(-1)   ! SP or L
   nbf0 = nbf0 + 4
  case(-2)   ! 5D
   n5dmark = n5dmark + 1
   d_mark(n5dmark) = nbf0 + 1
   nbf0 = nbf0 + 5
  case( 2)   ! 6D
   n6dmark = n6dmark + 1
   d_mark(n6dmark) = nbf0 + 1
   nbf0 = nbf0 + 6
  case(-3)   ! 7F
   n7fmark = n7fmark + 1
   f_mark(n7fmark) = nbf0 + 1
   nbf0 = nbf0 + 7
  case( 3)   ! 10F
   n10fmark = n10fmark + 1
   f_mark(n10fmark) = nbf0 + 1
   nbf0 = nbf0 + 10
  case(-4)   ! 9G
   n9gmark = n9gmark + 1
   g_mark(n9gmark) = nbf0 + 1
   nbf0 = nbf0 + 9
  case( 4)   ! 15G
   n15gmark = n15gmark + 1
   g_mark(n15gmark) = nbf0 + 1
   nbf0 = nbf0 + 15
  case(-5)   ! 11H
   n11hmark = n11hmark + 1
   h_mark(n11hmark) = nbf0 + 1
   nbf0 = nbf0 + 11
  case( 5)   ! 21H
   n21hmark = n21hmark + 1
   h_mark(n21hmark) = nbf0 + 1
   nbf0 = nbf0 + 21
  case default
   write(6,'(A)') 'ERROR in subroutine fch2py: shell_type(i) out of range.'
   write(6,'(A,3I4)') 'k, i, shell_type(i)=', k, i, shell_type(i)
   stop
  end select
 end do ! for i
 deallocate(shell_type)

 ! adjust the order of d, f, etc. functions
 call fch2py_permute_sph(n5dmark, n7fmark, n9gmark, n11hmark, k, d_mark, &
                         f_mark, g_mark, h_mark, nbf, nif, coeff2)
 call fch2py_permute_cart(n6dmark, n10fmark, n15gmark, n21hmark, k, d_mark, &
                          f_mark, g_mark, h_mark, nbf, nif, coeff2)

 if(ghf) then ! Beta real or Beta imaginary
  d_mark = d_mark + nbf0
  f_mark = f_mark + nbf0
  g_mark = g_mark + nbf0
  h_mark = h_mark + nbf0
  call fch2py_permute_sph(n5dmark, n7fmark, n9gmark, n11hmark, k, d_mark, &
                          f_mark, g_mark, h_mark, nbf, nif, coeff2)
  call fch2py_permute_cart(n6dmark, n10fmark, n15gmark, n21hmark, k, d_mark, &
                           f_mark, g_mark, h_mark, nbf, nif, coeff2)
 end if
! adjustment finished

 deallocate(d_mark, f_mark, g_mark, h_mark)
end subroutine fch2py

subroutine fch2py_permute_sph(n5dmark, n7fmark, n9gmark, n11hmark, k, d_mark, &
                              f_mark, g_mark, h_mark, nbf, nif, coeff2)
 implicit none
 integer :: i
 integer, intent(in) :: n5dmark, n7fmark, n9gmark, n11hmark, k, nbf, nif
 integer, intent(in) :: d_mark(k), f_mark(k), g_mark(k), h_mark(k)
 real(kind=8), intent(inout) :: coeff2(nbf,nif)

 if(n5dmark==0 .and. n7fmark==0 .and. n9gmark==0 .and. n11hmark==0) return

 do i = 1, n5dmark, 1
  call fch2py_permute_5d(nif, coeff2(d_mark(i):d_mark(i)+4,:))
 end do
 do i = 1, n7fmark, 1
  call fch2py_permute_7f(nif, coeff2(f_mark(i):f_mark(i)+6,:))
 end do
 do i = 1, n9gmark, 1
  call fch2py_permute_9g(nif, coeff2(g_mark(i):g_mark(i)+8,:))
 end do
 do i = 1, n11hmark, 1
  call fch2py_permute_11h(nif, coeff2(h_mark(i):h_mark(i)+10,:))
 end do

end subroutine fch2py_permute_sph

subroutine fch2py_permute_cart(n6dmark, n10fmark, n15gmark, n21hmark, k, d_mark, &
                              f_mark, g_mark, h_mark, nbf, nif, coeff2)
 implicit none
 integer :: i
 integer, intent(in) :: n6dmark, n10fmark, n15gmark, n21hmark, k, nbf, nif
 integer, intent(in) :: d_mark(k), f_mark(k), g_mark(k), h_mark(k)
 real(kind=8), intent(inout) :: coeff2(nbf,nif)

 if(n6dmark==0 .and. n10fmark==0 .and. n15gmark==0 .and. n21hmark==0) return

 do i = 1, n6dmark, 1
  call fch2py_permute_6d(nif, coeff2(d_mark(i):d_mark(i)+5,:))
 end do
 do i = 1, n10fmark, 1
  call fch2py_permute_10f(nif, coeff2(f_mark(i):f_mark(i)+9,:))
 end do
 do i = 1, n15gmark, 1
  call fch2py_permute_15g(nif, coeff2(g_mark(i):g_mark(i)+14,:))
 end do
 do i = 1, n21hmark, 1
  call fch2py_permute_21h(nif, coeff2(h_mark(i):h_mark(i)+20,:))
 end do

end subroutine fch2py_permute_cart

subroutine fch2py_permute_5d(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(5) = [5, 3, 1, 2, 4]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(5,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical d functions in Gaussian
! To: the order of spherical d functions in PySCF
! 1    2    3    4    5
! d0 , d+1, d-1, d+2, d-2
! d-2, d-1, d0 , d+1, d+2

 allocate(coeff2(5,nif), source=0d0)
 forall(i = 1:5) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
end subroutine fch2py_permute_5d

subroutine fch2py_permute_6d(nif,coeff)
 use Sdiag_parameter, only: Sdiag_d
 implicit none
 integer :: i
 integer, parameter :: order(6) = [1, 4, 5, 2, 6, 3]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(6,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian d functions in Gaussian
! To: the order of Cartesian d functions in PySCF
! 1  2  3  4  5  6
! XX,YY,ZZ,XY,XZ,YZ
! XX,XY,XZ,YY,YZ,ZZ

 allocate(coeff2(6,nif), source=coeff)
 forall(i = 1:6) coeff(i,:) = coeff2(order(i),:)/Sdiag_d(i)
 deallocate(coeff2)
end subroutine fch2py_permute_6d

subroutine fch2py_permute_7f(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(7) = [7, 5, 3, 1, 2, 4, 6]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(7,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical f functions in Gaussian
! To: the order of spherical f functions in PySCF
! 1    2    3    4    5    6    7
! f0 , f+1, f-1, f+2, f-2, f+3, f-3
! f-3, f-2, f-1, f0 , f+1, f+2, f+3

 allocate(coeff2(7,nif), source=0d0)
 forall(i = 1:7) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
end subroutine fch2py_permute_7f

subroutine fch2py_permute_10f(nif,coeff)
 use Sdiag_parameter, only: Sdiag_f
 implicit none
 integer :: i
 integer, parameter :: order(10) = [1, 5, 6, 4, 10, 7, 2, 9, 8, 3]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(10,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian f functions in Gaussian
! To: the order of Cartesian f functions in PySCF
! 1   2   3   4   5   6   7   8   9   10
! XXX,YYY,ZZZ,XYY,XXY,XXZ,XZZ,YZZ,YYZ,XYZ
! XXX,XXY,XXZ,XYY,XYZ,XZZ,YYY,YYZ,YZZ,ZZZ

 allocate(coeff2(10,nif), source=coeff)
 forall(i = 1:10) coeff(i,:) = coeff2(order(i),:)/Sdiag_f(i)
 deallocate(coeff2)
end subroutine fch2py_permute_10f

subroutine fch2py_permute_9g(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(9) = [9, 7, 5, 3, 1, 2, 4, 6, 8]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(9,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical g functions in Gaussian
! To: the order of spherical g functions in PySCF
! 1    2    3    4    5    6    7    8    9
! g0 , g+1, g-1, g+2, g-2, g+3, g-3, g+4, g-4
! g-4, g-3, g-2, g-1, g0 , g+1, g+2, g+3, g+4

 allocate(coeff2(9,nif), source=0d0)
 forall(i = 1:9) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
end subroutine fch2py_permute_9g

subroutine fch2py_permute_15g(nif,coeff)
 use Sdiag_parameter, only: Sdiag_g
 implicit none
 integer :: i
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(15,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian g functions in Gaussian
! To: the order of Cartesian g functions in PySCF
! 1    2    3    4    5    6    7    8    9    10   11   12   13   14   15
! ZZZZ,YZZZ,YYZZ,YYYZ,YYYY,XZZZ,XYZZ,XYYZ,XYYY,XXZZ,XXYZ,XXYY,XXXZ,XXXY,XXXX
! xxxx,xxxy,xxxz,xxyy,xxyz,xxzz,xyyy,xyyz,xyzz,xzzz,yyyy,yyyz,yyzz,yzzz,zzzz

 allocate(coeff2(15,nif), source=coeff)
 forall(i = 1:15) coeff(i,:) = coeff2(16-i,:)/Sdiag_g(i)
 deallocate(coeff2)
end subroutine fch2py_permute_15g

subroutine fch2py_permute_11h(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(11) = [11, 9, 7, 5, 3, 1, 2, 4, 6, 8, 10]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(11,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical h functions in Gaussian
! To: the order of spherical h functions in PySCF
! 1    2    3    4    5    6    7    8    9    10   11
! h0 , h+1, h-1, h+2, h-2, h+3, h-3, h+4, h-4, h+5, h-5
! h-5, h-4, h-3, h-2, h-1, h0 , h+1, h+2, h+3, h+4, h+5

 allocate(coeff2(11,nif), source=0d0)
 forall(i = 1:11) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
end subroutine fch2py_permute_11h

subroutine fch2py_permute_21h(nif,coeff)
 use Sdiag_parameter, only: Sdiag_h
 implicit none
 integer :: i
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(21,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian h functions in Gaussian
! To: the order of Cartesian h functions in PySCF
! 1     2     3     4     5     6     7     8     9     10    11    12    13    14    15    16    17    18    19    20    21
! ZZZZZ,YZZZZ,YYZZZ,YYYZZ,YYYYZ,YYYYY,XZZZZ,XYZZZ,XYYZZ,XYYYZ,XYYYY,XXZZZ,XXYZZ,XXYYZ,XXYYY,XXXZZ,XXXYZ,XXXYY,XXXXZ,XXXXY,XXXXX
! xxxxx,xxxxy,xxxxz,xxxyy,xxxyz,xxxzz,xxyyy,xxyyz,xxyzz,xxzzz,xyyyy,xyyyz,xyyzz,xyzzz,xzzzz,yyyyy,yyyyz,yyyzz,yyzzz,yzzzz,zzzzz

 allocate(coeff2(21,nif), source=coeff)
 forall(i = 1:21) coeff(i,:) = coeff2(22-i,:)/Sdiag_h(i)
 deallocate(coeff2)
end subroutine fch2py_permute_21h

! complex GHF wrapper of fch2py
subroutine fch2py_cghf(fchname, nbf, nif, coeff)
 implicit none
 integer :: i, j
 integer, intent(in) :: nbf, nif
!f2py intent(in) :: nbf, nif

 complex(kind=8), intent(out) :: coeff(nbf,nif)
!f2py intent(out) :: coeff
!f2py depend(nbf,nif) :: coeff
 real(kind=8), allocatable :: real_c(:,:), imag_c(:,:)

 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 allocate(real_c(nbf,nif), source=0d0)
 allocate(imag_c(nbf,nif), source=0d0)

 call fch2py(fchname, nbf, nif, 'r', real_c) ! real part
 call fch2py(fchname, nbf, nif, 'i', imag_c) ! imaginary part
 forall(i=1:nbf, j=1:nif) coeff(i,j) = CMPLX(real_c(i,j), imag_c(i,j))

 deallocate(real_c, imag_c)
end subroutine fch2py_cghf

