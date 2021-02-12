! written by jxzou at 20210210: transfer MOs from Gaussian to Psi4

! This utility/subroutine will generate PSI4 input file from Gaussian .fch(k)
! file (with coordinates, basis sets and MOs written in)

program main
 implicit none
 integer :: i
 integer, parameter :: iout = 6
 character(len=4) :: str
 character(len=240) :: fchname = ' '
 logical :: uhf

 i = iargc()
 if(i<1 .or. i>2) then
  write(iout,'(/,A)') ' ERROR in subroutine fch2psi: wrong command line arguments!'
  write(iout,'(A)')   ' Example 1 (R(O)HF, CAS): fch2psi a.fch'
  write(iout,'(A,/)') ' Example 2 (UHF)        : fch2psi a.fch -uhf'
  stop
 end if

 call getarg(1, fchname)
 call require_file_exist(fchname)

 str = ' '
 uhf = .false.

 if(i == 2) then
  call getarg(2, str)
  if(str == '-uhf') then
   uhf = .true.
  else
   write(iout,'(A)') 'ERROR in subroutine fch2psi: wrong command line arguments.'
   write(iout,'(A)') "The 2nd argument can only be '-uhf'. But got '"//str//"'"
   stop
  end if
 end if

 call fch2psi(fchname, uhf)
 stop
end program main

! transfer MOs from Gaussian to Psi4
subroutine fch2psi(fchname, uhf)
 use util_wrapper, only: fch2inp_wrap
 implicit none
 integer :: i, k, fid, system
 integer :: nbf0, nbf, nif, ncontr, n3pmark
 integer, parameter :: iout = 6
 integer, allocatable :: shell_type(:), shl2atm(:), p_mark(:)
 real(kind=8), allocatable :: coeff(:,:)
 character(len=240) :: inpname, fileA, fileB
 character(len=240), intent(in) :: fchname
 logical, intent(in) :: uhf
 logical :: sph

 call fch2inp_wrap(fchname, uhf, .false., 0, 0)
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

 n3pmark = 0
 allocate(p_mark(ncontr), source=0)
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
  case(-6)   ! 13I
   nbf = nbf + 13
  case( 6)   ! 28I
   nbf = nbf + 28
  end select
 end do ! for i
 deallocate(shell_type)

 if(nbf0 /= nbf) then
  write(iout,'(A)') 'ERROR in subroutine fch2psi: inconsistent nbf.'
  write(iout,'(2(A,I0))') 'nbf0=', nbf0, ', nbf=', nbf
  stop
 end if

 ! adjust MO coefficients according to the order of AOs (p,)
 do i = 1, n3pmark, 1
  call fch2psi_permute_3p(nif, coeff(p_mark(i):p_mark(i)+2,:))
 end do ! for i
 deallocate(p_mark)

 if(uhf) nif = nif/2
 call write_mo_into_psi_mat(fileA, nbf, nif, coeff(:,1:nif))
 if(uhf) call write_mo_into_psi_mat(fileB, nbf, nif, coeff(:,nif+1:2*nif))

 deallocate(coeff)
 return
end subroutine fch2psi

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
! Px,  Py,  Pz
! Pz,  Px,  Py

 allocate(coeff2(3,nif), source=0d0)
 forall(i = 1:3) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
 return
end subroutine fch2psi_permute_3p

