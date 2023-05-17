! written by jxzou at 20230516: generate the wave function of the mirror of a
! molecule

! .chk -> .chk, a wrapper of {formchk, mirror_wfn, and unfchk}
subroutine mirror_c2c(chkname)
 use util_wrapper, only: formchk, unfchk
 implicit none
 integer :: i
 character(len=240) :: fchname1, fchname2
 character(len=240), intent(in) :: chkname
!f2py intent(in) :: chkname

 call formchk(chkname)
 i = index(chkname, '.chk', back=.true.)
 fchname1 = chkname(1:i-1)//'.fch'
 fchname2 = chkname(1:i-1)//'_m.fch'
 call mirror_wfn(fchname1)
 call unfchk(fchname2)
end subroutine mirror_c2c

! Get the wave function of the mirror of a molecule.
! The z-component of Cartesian coordinates are multiplied by -1.
! Some of the MO coefficients on basis functions are multiplied by -1
!  5D: D+1, D-1
!  7F: F0, F+2, F-2
!  9G: G+1, G-1, G+3, G-3
! 11H: H+2, H-2, H+4, H-4
subroutine mirror_wfn(fchname)
 use fch_content, only: nbf, nif, ncontr, alpha_coeff, beta_coeff, shell_type, &
  is_uhf, coor, tot_dm, spin_dm, check_uhf_in_fch, read_fch
 implicit none
 integer :: i, j, k, nif1
 integer, parameter :: a1(3) = [1,4,5]
 integer, parameter :: a2(4) = [3,6,9,10]
 integer, parameter :: a3(4) =[2,3,6,7]
 integer, parameter :: a4(6) = [2,4,6,8,11,13]
 integer, parameter :: a5(4) = [4,5,8,9]
 integer, parameter :: a6(9) = [1,3,5,8,10,12,14,17,19]
 character(len=240) :: new_fch
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname
 real(kind=8), allocatable :: coeff(:,:)
 logical, allocatable :: pos(:)
 logical :: uhf

 call require_file_exist(fchname)
 i = index(fchname, '.fch', back=.true.)
 if(i == 0) then
  write(6,'(/,A)') "ERROR in subroutine mirror_wfn: '.fch' suffix not found in&
                   & file "//TRIM(fchname)
  stop
 end if
 new_fch = fchname(1:i-1)//'_m.fch'

 call check_uhf_in_fch(fchname, uhf)
 is_uhf = uhf ! be careful of this variable

 call read_fch(fchname, uhf)
 coor(3,:) = -coor(3,:) ! z*(-1)

 nif1 = nif
 if(uhf) nif1 = 2*nif
 allocate(coeff(nbf,nif1))
 if(uhf) then
  coeff(:,1:nif) = alpha_coeff
  coeff(:,nif+1:nif1) = beta_coeff
 else
  coeff = alpha_coeff
 end if

 allocate(pos(nbf)) ! record the positions to be positive or to be negative
 pos = .true.
 k = 0

 do i = 1, ncontr, 1
  select case(shell_type(i))
  case( 0)   ! S
   k = k + 1
  case( 1)   ! 3P
   k = k + 3
   pos(k) = .false.
  case(-1)   ! SP or L
   k = k + 4
   pos(k) = .false.
  case(-2)   ! 5D
   pos(k+2:k+3) = .false.
   k = k + 5
  case( 2)   ! 6D
   k = k + 6
   pos(k-1:k) = .false.
  case(-3)   ! 7F
   forall(j = 1:3) pos(k+a1(j)) = .false.
   k = k + 7
  case( 3)   ! 10F
   forall(j = 1:4) pos(k+a2(j)) = .false.
   k = k + 10
  case(-4)   ! 9G
   forall(j = 1:4) pos(k+a3(j)) = .false.
   k = k + 9
  case( 4)   ! 15G
   forall(j = 1:6) pos(k+a4(j)) = .false.
   k = k + 15
  case(-5)   ! 11H
   forall(j = 1:4) pos(k+a5(j)) = .false.
   k = k + 11
  case( 5)   ! 21H
   forall(j = 1:9) pos(k+a6(j)) = .false.
   k = k + 21
  case default
   write(6,'(/,A,I0)') 'ERROR in subroutine mirror_wfn: invalid shell_type(i)=',&
                        shell_type(i)
   stop
  end select
 end do ! for i

 if(k /= nbf) then
  write(6,'(/,A)') 'ERROR in subroutine mirror_wfn: internal inconsistency.'
  write(6,'(A)') 'fchname='//TRIM(fchname)
  write(6,'(2(A,I0))') 'nbf=', nbf, ', k=', k
  stop
 end if

 ! adjust MO coefficients
 forall(i=1:nbf, (.not.pos(i))) coeff(i,:) = -coeff(i,:)
 if(uhf) then
  alpha_coeff = coeff(:,1:nif)
  beta_coeff = coeff(:,nif+1:nif1)
 else
  alpha_coeff = coeff
 end if
 deallocate(coeff)

 ! adjust Total SCF Density
 forall(i=1:nbf, j=1:nbf, (pos(i).neqv.pos(j))) tot_dm(j,i) = -tot_dm(j,i)

 ! adjust Spin SCF Density, if any
 if(allocated(spin_dm)) then
  forall(i=1:nbf, j=1:nbf, (pos(i).neqv.pos(j)))
   spin_dm(j,i) = -spin_dm(j,i)
  end forall
 end if
 deallocate(pos)

 ! generate a new .fch file
 call write_fch(new_fch)
end subroutine mirror_wfn

