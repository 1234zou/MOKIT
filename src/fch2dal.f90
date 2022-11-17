! written by jxzou at 20210613: adjust the orders of d,f,g, etc functions in
!  .fch(k) file, to the orders in Dalton

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
 use util_wrapper, only: formchk, fch2inp_wrap
 implicit none
 integer :: i, system
 character(len=240) :: fchname, inpname
 logical :: sph

 i = iargc()
 if(i /= 1) then
  write(6,'(/,A)') ' ERROR in subroutine fch2dal: wrong command line argument!'
  write(6,'(A,/)') ' Example (R(O)HF, CAS): fch2dal a.fch'
  stop
 end if

 fchname = ' '
 sph = .true.

 call getarg(1, fchname)
 call require_file_exist(fchname)

 ! if .chk file provided, convert into .fch file automatically
 i = LEN_TRIM(fchname)
 if(fchname(i-3:i) == '.chk') then
  call formchk(fchname)
  fchname = fchname(1:i-3)//'fch'
 end if

 i = index(fchname, '.fch', back=.true.)
 inpname = fchname(1:i-1)//'.inp'

 call fch2inp_wrap(fchname, .false., 0, 0) ! generate GAMESS .inp file

 call check_sph(fchname, sph)
 if(sph) then
  i = system('bas_gms2dal '//TRIM(inpname)//' -sph')
 else
  i = system('bas_gms2dal '//TRIM(inpname))
 end if

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine fch2dal: failed to call utility&
                  & bas_gms2dal. Two possible reasons:'
  write(6,'(A)')   '(1) The file '//TRIM(fchname)//' may be incomplete.'
  write(6,'(A)')   '(2) You forgot to compile the utility bas_gms2dal.'
  write(6,'(A,/)') '(3) This is a bug of the utility bas_gms2dal.'
  stop
 end if

 call delete_file(inpname)
 call fch2dal(fchname)
 stop
end program main

! read the MOs in .fch(k) file and adjust its d,f,g,h functions order of Gaussian
!  to that of Dalton
subroutine fch2dal(fchname)
 use fch_content, only: check_uhf_in_fch
 implicit none
 integer :: i, k, length, fid, fid1, nbf, nif, RENAME
 integer :: n6dmark,n10fmark,n15gmark,n21hmark
 integer :: n5dmark,n7fmark, n9gmark, n11hmark
 integer, allocatable :: shell_type(:), shell2atom_map(:)
 ! mark the index where d, f, g, h functions begin
 integer, allocatable :: d_mark(:), f_mark(:), g_mark(:), h_mark(:)
 real(kind=8), allocatable :: coeff(:,:)
 character(len=24) :: data_string
 character(len=240) :: buf, dalfile, dalfile1
 character(len=240), intent(in) :: fchname
 logical :: uhf, sph

 buf = ' '
 call check_uhf_in_fch(fchname, uhf)
 if(uhf) then
  write(6,'(/,A)') 'ERROR in subroutine fch2dal: UHF wave function not supported.'
  write(6,'(A)') 'Because there is no UHF method in Dalton. You can compute&
                   & ROHF instead.'
  stop
 end if

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(coeff(nbf,nif))
 call read_mo_from_fch(fchname, nbf, nif, 'a', coeff)

 call read_ncontr_from_fch(fchname, k)
 allocate(shell_type(2*k), source=0)
 allocate(shell2atom_map(2*k), source=0)
 call read_shltyp_and_shl2atm_from_fch(fchname, k, shell_type, shell2atom_map)

 if(ANY(shell_type<-1) .and. ANY(shell_type>1)) then
  write(6,'(A)') 'ERROR in subroutine fch2dal: mixed spherical harmonic/&
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
! this is to ensure that D comes after L functions
 ! split the 'L' into 'S' and 'P'
 call split_L_func(k, shell_type, shell2atom_map, length)

 ! sort the shell_type, shell_to_atom_map by ascending order
 ! MOs will be adjusted accordingly
 call sort_shell_and_mo(length, shell_type, shell2atom_map, nbf, nif, coeff)
 deallocate(shell2atom_map)
! adjust done

! then we adjust the basis functions in each MO according to the type of basis functions
 k = length  ! update k
 n6dmark = 0; n10fmark = 0; n15gmark = 0; n21hmark = 0
 n5dmark = 0; n7fmark = 0; n9gmark = 0; n11hmark = 0
 allocate(d_mark(k), source=0)
 allocate(f_mark(k), source=0)
 allocate(g_mark(k), source=0)
 allocate(h_mark(k), source=0)
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
 deallocate(shell_type)

 ! adjust the order of d, f, etc. functions
 do i = 1, n5dmark, 1
  call fch2dal_permute_5d(nif,coeff(d_mark(i):d_mark(i)+4,:))
 end do
 do i = 1, n6dmark, 1
  call fch2dal_permute_6d(nif,coeff(d_mark(i):d_mark(i)+5,:))
 end do
 do i = 1, n7fmark, 1
  call fch2dal_permute_7f(nif,coeff(f_mark(i):f_mark(i)+6,:))
 end do
 do i = 1, n10fmark, 1
  call fch2dal_permute_10f(nif,coeff(f_mark(i):f_mark(i)+9,:))
 end do
 do i = 1, n9gmark, 1
  call fch2dal_permute_9g(nif,coeff(g_mark(i):g_mark(i)+8,:))
 end do
 do i = 1, n15gmark, 1
  call fch2dal_permute_15g(nif,coeff(g_mark(i):g_mark(i)+14,:))
 end do
 do i = 1, n11hmark, 1
  call fch2dal_permute_11h(nif,coeff(h_mark(i):h_mark(i)+10,:))
 end do
 do i = 1, n21hmark, 1
  call fch2dal_permute_21h(nif,coeff(h_mark(i):h_mark(i)+20,:))
 end do
! adjustment finished

 deallocate(d_mark, f_mark, g_mark, h_mark)

 i = index(fchname, '.fch', back=.true.)
 dalfile = fchname(1:i-1)//'.dal'
 dalfile1 = fchname(1:i-1)//'.dal1'
 open(newunit=fid,file=TRIM(dalfile),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(dalfile1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:9) == '.PUNCHOUT') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine fch2dal: no '' found&
                   & in file "//TRIM(dalfile)
  close(fid)
  close(fid1,status='delete')
  stop
 end if
 close(fid,status='delete')

 write(fid1,'(A)') '.MOSTART'
 write(fid1,'(A)') 'FORM18'
 write(fid1,'(A)') '.PUNCHOUTPUTORBITALS'
 data_string = ' '
 call fdate(data_string)
 write(fid1,'(A)') '**MOLORB (punched by fch2dal of MOKIT '//TRIM(data_string)//')'
 do i = 1, nif, 1
  write(fid1,'(4F18.14)') (coeff(k,i),k=1,nbf)
 end do ! for i
 deallocate(coeff)

 write(fid1,'(A)') '**END OF INPUT'
 close(fid1)
 i = RENAME(TRIM(dalfile1), TRIM(dalfile))
 return
end subroutine fch2dal

subroutine fch2dal_permute_5d(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(5) = [5, 3, 1, 2, 4]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(5,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical d functions in Gaussian
! To: the order of spherical d functions in Dalton
! 1    2    3    4    5
! d0 , d+1, d-1, d+2, d-2
! d-2, d-1, d0 , d+1, d+2

 allocate(coeff2(5,nif), source=0d0)
 forall(i = 1:5) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
 return
end subroutine fch2dal_permute_5d

subroutine fch2dal_permute_6d(nif,coeff)
 use Sdiag_dalton, only: Sdiag_d
 implicit none
 integer :: i
 integer, parameter :: order(6) = [1, 4, 5, 2, 6, 3]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(6,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian d functions in Gaussian
! To: the order of Cartesian d functions in Dalton
! 1  2  3  4  5  6
! XX,YY,ZZ,XY,XZ,YZ
! XX,XY,XZ,YY,YZ,ZZ

 allocate(coeff2(6,nif), source=coeff)
 forall(i = 1:6) coeff(i,:) = coeff2(order(i),:)/Sdiag_d(i)
 deallocate(coeff2)
 return
end subroutine fch2dal_permute_6d

subroutine fch2dal_permute_7f(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(7) = [7, 5, 3, 1, 2, 4, 6]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(7,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical f functions in Gaussian
! To: the order of spherical f functions in Dalton
! 1    2    3    4    5    6    7
! f0 , f+1, f-1, f+2, f-2, f+3, f-3
! f-3, f-2, f-1, f0 , f+1, f+2, f+3

 allocate(coeff2(7,nif), source=0d0)
 forall(i = 1:7) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
 return
end subroutine fch2dal_permute_7f

subroutine fch2dal_permute_10f(nif,coeff)
 use Sdiag_dalton, only: Sdiag_f
 implicit none
 integer :: i
 integer, parameter :: order(10) = [1, 5, 6, 4, 10, 7, 2, 9, 8, 3]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(10,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian f functions in Gaussian
! To: the order of Cartesian f functions in Dalton
! 1   2   3   4   5   6   7   8   9   10
! XXX,YYY,ZZZ,XYY,XXY,XXZ,XZZ,YZZ,YYZ,XYZ
! XXX,XXY,XXZ,XYY,XYZ,XZZ,YYY,YYZ,YZZ,ZZZ

 allocate(coeff2(10,nif), source=coeff)
 forall(i = 1:10) coeff(i,:) = coeff2(order(i),:)/Sdiag_f(i)
 deallocate(coeff2)
 return
end subroutine fch2dal_permute_10f

subroutine fch2dal_permute_9g(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(9) = [9, 7, 5, 3, 1, 2, 4, 6, 8]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(9,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical g functions in Gaussian
! To: the order of spherical g functions in Dalton
! 1    2    3    4    5    6    7    8    9
! g0 , g+1, g-1, g+2, g-2, g+3, g-3, g+4, g-4
! g-4, g-3, g-2, g-1, g0 , g+1, g+2, g+3, g+4

 allocate(coeff2(9,nif), source=0d0)
 forall(i = 1:9) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
 return
end subroutine fch2dal_permute_9g

subroutine fch2dal_permute_15g(nif,coeff)
 use Sdiag_dalton, only: Sdiag_g
 implicit none
 integer :: i
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(15,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian g functions in Gaussian
! To: the order of Cartesian g functions in Dalton
! 1    2    3    4    5    6    7    8    9    10   11   12   13   14   15
! ZZZZ,YZZZ,YYZZ,YYYZ,YYYY,XZZZ,XYZZ,XYYZ,XYYY,XXZZ,XXYZ,XXYY,XXXZ,XXXY,XXXX
! xxxx,xxxy,xxxz,xxyy,xxyz,xxzz,xyyy,xyyz,xyzz,xzzz,yyyy,yyyz,yyzz,yzzz,zzzz

 allocate(coeff2(15,nif), source=coeff)
 forall(i = 1:15) coeff(i,:) = coeff2(16-i,:)/Sdiag_g(i)
 deallocate(coeff2)
 return
end subroutine fch2dal_permute_15g

subroutine fch2dal_permute_11h(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(11) = [11, 9, 7, 5, 3, 1, 2, 4, 6, 8, 10]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(11,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical h functions in Gaussian
! To: the order of spherical h functions in Dalton
! 1    2    3    4    5    6    7    8    9    10   11
! h0 , h+1, h-1, h+2, h-2, h+3, h-3, h+4, h-4, h+5, h-5
! h-5, h-4, h-3, h-2, h-1, h0 , h+1, h+2, h+3, h+4, h+5

 allocate(coeff2(11,nif), source=0d0)
 forall(i = 1:11) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
 return
end subroutine fch2dal_permute_11h

subroutine fch2dal_permute_21h(nif,coeff)
 use Sdiag_dalton, only: Sdiag_h
 implicit none
 integer :: i
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(21,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian h functions in Gaussian
! To: the order of Cartesian h functions in Dalton
! 1     2     3     4     5     6     7     8     9     10    11    12    13    14    15    16    17    18    19    20    21
! ZZZZZ,YZZZZ,YYZZZ,YYYZZ,YYYYZ,YYYYY,XZZZZ,XYZZZ,XYYZZ,XYYYZ,XYYYY,XXZZZ,XXYZZ,XXYYZ,XXYYY,XXXZZ,XXXYZ,XXXYY,XXXXZ,XXXXY,XXXXX
! xxxxx,xxxxy,xxxxz,xxxyy,xxxyz,xxxzz,xxyyy,xxyyz,xxyzz,xxzzz,xyyyy,xyyyz,xyyzz,xyzzz,xzzzz,yyyyy,yyyyz,yyyzz,yyzzz,yzzzz,zzzzz

 allocate(coeff2(21,nif), source=coeff)
 forall(i = 1:21) coeff(i,:) = coeff2(22-i,:)/Sdiag_h(i)
 deallocate(coeff2)
 return
end subroutine fch2dal_permute_21h

