! written by jxzou at 20210210: transfer MOs from Gaussian to Psi4
! updated by jxzou at 20210215: supported ECP/PP

! This utility/subroutine will generate PSI4 input file from Gaussian .fch(k)
! file (with coordinates, basis sets and MOs written in)

module root_param_fch2psi
 implicit none
 real(kind=8), parameter :: root3 = DSQRT(3d0)
 real(kind=8), parameter :: root5 = DSQRT(5d0)
 real(kind=8), parameter :: root7 = DSQRT(7d0)
 real(kind=8), parameter :: root15 = DSQRT(15d0)
 real(kind=8), parameter :: root21 = DSQRT(21d0)
 real(kind=8), parameter :: root35 = DSQRT(35d0)
 real(kind=8), parameter :: root35_3 = DSQRT(35d0/3d0)
 real(kind=8), parameter :: root63 = 3d0*DSQRT(7d0)
 real(kind=8), parameter :: root105 = DSQRT(105d0)
end module root_param_fch2psi

program main
 implicit none
 integer :: i
 integer, parameter :: iout = 6
 character(len=240) :: fchname = ' '

 i = iargc()
 if(i /= 1) then
  write(iout,'(/,A)') ' ERROR in subroutine fch2psi: wrong command line arguments!'
  write(iout,'(A,/)') ' Example (R(O)HF, UHF, CAS): fch2psi a.fch'
  stop
 end if

 call getarg(1, fchname)
 call require_file_exist(fchname)
 call fch2psi(fchname)
 stop
end program main

! transfer MOs from Gaussian to Psi4
subroutine fch2psi(fchname)
 use util_wrapper, only: fch2inp_wrap
 use fch_content, only: check_uhf_in_fch
 implicit none
 integer :: i, k, fid, system
 integer :: nbf0, nbf, nif, ncontr
 integer :: n3pmark, n6dmark, n10fmark, n15gmark, n21hmark
 integer, parameter :: iout = 6
 integer, allocatable :: shell_type(:), shl2atm(:)
 integer, allocatable :: p_mark(:), d_mark(:), f_mark(:), g_mark(:), h_mark(:)
 real(kind=8), allocatable :: coeff(:,:)
 character(len=240) :: inpname, fileA, fileB
 character(len=240), intent(in) :: fchname
 logical :: sph, uhf

 call check_uhf_in_fch(fchname, uhf)
 call fch2inp_wrap(fchname, .false., 0, 0)
 call check_sph(fchname, sph)

 i = index(fchname,'.fch', back=.true.)
 inpname = fchname(1:i-1)//'.inp'
 fileA = fchname(1:i-1)//'.A'
 fileB = fchname(1:i-1)//'.B'

 if(sph) then
  i = system('bas_gms2psi '//TRIM(inpname)//' -sph')
 else
  i = system('bas_gms2psi '//TRIM(inpname))
 end if

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine fch2psi: call utility bas_gms2psi failed.'
  write(iout,'(A)') 'Did you forget to compile bas_gms2psi? Or the file '//&
                     TRIM(fchname)//' may be incomplete.'
  stop
 end if
 call delete_file(inpname)

 i = index(fchname,'.fch', back=.true.)
 inpname = fchname(1:i-1)//'_psi.inp'

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 nbf0 = nbf

 ! read MO Coefficients
 if(uhf) then
  allocate(coeff(nbf,2*nif))
  call read_mo_from_fch(fchname, nbf, nif, 'a', coeff(:,1:nif))
  call read_mo_from_fch(fchname, nbf, nif, 'b', coeff(:,nif+1:2*nif))
  nif = 2*nif   ! double the size
 else
  allocate(coeff(nbf,nif))
  call read_mo_from_fch(fchname, nbf, nif, 'a', coeff)
 end if

 call read_ncontr_from_fch(fchname, ncontr)
 allocate(shell_type(ncontr), source=0)
 allocate(shl2atm(ncontr), source=0)
 call read_shltyp_and_shl2atm_from_fch(fchname, ncontr, shell_type, shl2atm)
 deallocate(shl2atm)

 n3pmark = 0; n6dmark = 0; n10fmark = 0; n15gmark = 0; n21hmark = 0
 allocate(p_mark(ncontr), source=0)
 allocate(d_mark(ncontr), source=0)
 allocate(f_mark(ncontr), source=0)
 allocate(g_mark(ncontr), source=0)
 allocate(h_mark(ncontr), source=0)
 nbf = 0
 do i = 1, ncontr, 1
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
   nbf = nbf + 5
  case( 2)   ! 6D
   n6dmark = n6dmark + 1
   d_mark(n6dmark) = nbf + 1
   nbf = nbf + 6
  case(-3)   ! 7F
   nbf = nbf + 7
  case( 3)   ! 10F
   n10fmark = n10fmark + 1
   f_mark(n10fmark) = nbf + 1
   nbf = nbf + 10
  case(-4)   ! 9G
   nbf = nbf + 9
  case( 4)   ! 15G
   n15gmark = n15gmark + 1
   g_mark(n15gmark) = nbf + 1
   nbf = nbf + 15
  case(-5)   ! 11H
   nbf = nbf + 11
  case( 5)   ! 21H
   n21hmark = n21hmark + 1
   h_mark(n21hmark) = nbf + 1
   nbf = nbf + 21
  case default
   write(iout,'(A)') 'ERROR in subroutine fch2psi: angular momentum too high!'
   write(iout,'(3(A,I0))') 'ncontr=', ncontr, ', i=', i, ', nbf=', nbf
   stop
  end select
 end do ! for i
 deallocate(shell_type)

 if(nbf0 /= nbf) then
  write(iout,'(A)') 'ERROR in subroutine fch2psi: inconsistent nbf.'
  write(iout,'(2(A,I0))') 'nbf0=', nbf0, ', nbf=', nbf
  stop
 end if

 ! adjust MO coefficients according to the order of AOs (p,)
 if(sph) then
  do i = 1, n3pmark, 1
   call fch2psi_permute_3p(nif, coeff(p_mark(i):p_mark(i)+2,:))
  end do ! for i
 end if
 do i = 1, n6dmark, 1
  call fch2psi_permute_6d(nif, coeff(d_mark(i):d_mark(i)+5,:))
 end do ! for i
 do i = 1, n10fmark, 1
  call fch2psi_permute_10f(nif, coeff(f_mark(i):f_mark(i)+9,:))
 end do ! for i
 do i = 1, n15gmark, 1
  call fch2psi_permute_15g(nif, coeff(g_mark(i):g_mark(i)+14,:))
 end do ! for i
 do i = 1, n21hmark, 1
  call fch2psi_permute_21h(nif, coeff(h_mark(i):h_mark(i)+20,:))
 end do ! for i
 deallocate(p_mark, d_mark, f_mark, g_mark, h_mark)

 if(uhf) nif = nif/2
 call write_mo_into_psi_mat(fileA, nbf, nif, coeff(:,1:nif))
 if(uhf) call write_mo_into_psi_mat(fileB, nbf, nif, coeff(:,nif+1:2*nif))

 deallocate(coeff)
 return
end subroutine fch2psi

! this subroutine is used only when spherical harmonic functions are used
subroutine fch2psi_permute_3p(nif, coeff)
 implicit none
 integer :: i
 integer, parameter :: order(3) = [3,1,2]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(3,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical p functions in Gaussian
! To: the order of spherical p functions in PSI4
! 1    2    3
! +1, -1,  0
! 0 , +1, -1

 allocate(coeff2(3,nif), source=0d0)
 forall(i = 1:3) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
 return
end subroutine fch2psi_permute_3p

! The basis order of 6D, 10F, 15G, 21H of PSI4 is the same as that of (OpenMolcas).
! But their normalization factors are not the same.
subroutine fch2psi_permute_6d(nif, coeff)
 use root_param_fch2psi, only: root3
 implicit none
 integer :: i, j
 integer, parameter :: order(6) = [1, 4, 5, 2, 6, 3]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(6,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian d functions in Gaussian
! To: the order of Cartesian d functions in PSI4
! 1  2  3  4  5  6
! XX,YY,ZZ,XY,XZ,YZ
! XX,XY,XZ,YY,YZ,ZZ

 forall(j=4:6, i=1:nif) coeff(j,i) = coeff(j,i)*root3

 allocate(coeff2(6,nif), source=0d0)
 forall(i = 1:6) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
 return
end subroutine fch2psi_permute_6d

subroutine fch2psi_permute_10f(nif, coeff)
 use root_param_fch2psi, only: root5, root15
 implicit none
 integer :: i, j
 integer, parameter :: order(10) = [1, 5, 6, 4, 10, 7, 2, 9, 8, 3]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(10,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian f functions in Gaussian
! To: the order of Cartesian f functions in PSI4
! 1   2   3   4   5   6   7   8   9   10
! XXX,YYY,ZZZ,XYY,XXY,XXZ,XZZ,YZZ,YYZ,XYZ
! XXX,XXY,XXZ,XYY,XYZ,XZZ,YYY,YYZ,YZZ,ZZZ

 forall(j=4:9, i=1:nif) coeff(j,i) = coeff(j,i)*root5
 coeff(10,:) = coeff(10,:)*root15

 allocate(coeff2(10,nif), source=0d0)
 forall(i = 1:10) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
 return
end subroutine fch2psi_permute_10f

subroutine fch2psi_permute_15g(nif, coeff)
 use root_param_fch2psi, only: root7, root35, root35_3
 implicit none
 integer :: i, j
 integer, intent(in) :: nif
 real(kind=8), parameter :: ratio(15) = [1d0,root7,root35_3,root7,1d0,root7,&
  root35,root35,root7,root35_3,root35,root35_3,root7,root7,1d0]
 real(kind=8), intent(inout) :: coeff(15,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian g functions in Gaussian
! To: the order of Cartesian g functions in PSI4
! 1    2    3    4    5    6    7    8    9    10   11   12   13   14   15
! ZZZZ,YZZZ,YYZZ,YYYZ,YYYY,XZZZ,XYZZ,XYYZ,XYYY,XXZZ,XXYZ,XXYY,XXXZ,XXXY,XXXX
! xxxx,xxxy,xxxz,xxyy,xxyz,xxzz,xyyy,xyyz,xyzz,xzzz,yyyy,yyyz,yyzz,yzzz,zzzz

 forall(j=1:15, i=1:nif) coeff(j,i) = coeff(j,i)*ratio(j)

 allocate(coeff2(15,nif), source=0d0)
 forall(i = 1:15) coeff2(i,:) = coeff(16-i,:)
 coeff = coeff2
 deallocate(coeff2)
 return
end subroutine fch2psi_permute_15g

subroutine fch2psi_permute_21h(nif, coeff)
 use root_param_fch2psi, only: root21, root63, root105
 implicit none
 integer :: i, j
 integer, intent(in) :: nif
 real(kind=8), parameter :: ratio(21) = [1d0,3d0,root21,root21,3d0,1d0,3d0,&
  root63,root105,root63,3d0,root21,root105,root105,root21,root21,root63,root21,&
  3d0,3d0,1d0]
 real(kind=8), intent(inout) :: coeff(21,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian h functions in Gaussian
! To: the order of Cartesian h functions in PSI4
! 1     2     3     4     5     6     7     8     9     10    11    12    13    14    15    16    17    18    19    20    21
! ZZZZZ,YZZZZ,YYZZZ,YYYZZ,YYYYZ,YYYYY,XZZZZ,XYZZZ,XYYZZ,XYYYZ,XYYYY,XXZZZ,XXYZZ,XXYYZ,XXYYY,XXXZZ,XXXYZ,XXXYY,XXXXZ,XXXXY,XXXXX
! xxxxx,xxxxy,xxxxz,xxxyy,xxxyz,xxxzz,xxyyy,xxyyz,xxyzz,xxzzz,xyyyy,xyyyz,xyyzz,xyzzz,xzzzz,yyyyy,yyyyz,yyyzz,yyzzz,yzzzz,zzzzz

 forall(j=1:21, i=1:nif) coeff(j,i) = coeff(j,i)*ratio(j)

 allocate(coeff2(21,nif), source=0d0)
 forall(i = 1:21) coeff2(i,:) = coeff(22-i,:)
 coeff = coeff2
 deallocate(coeff2)
 return
end subroutine fch2psi_permute_21h

