! written by jxzou at 20180208
! updated by jxzou at 20200402: renamed from fchk_mo_copy to fch_mo_copy
! updated by jxzou at 20200622: simplify code, add support input idx1, idx2, and combined with rwwfn.f90

! This is a program/subroutine used to copy MOs from one .fchk file into another .fchk file

program main
 implicit none
 integer :: i, idx1, idx2
 character(len=3) :: ab
 character(len=8) :: str
 character(len=240) :: fchname1, fchname2

 idx1 = 0; idx2 = 0
 fchname1 = ' '; fchname2 = ' '
 ab = '-aa'

 i = iargc()
 if(i<2 .or. i>4) then
  write(6,'(/,A)') ' ERROR in subroutine fch_mo_copy: wrong command line arguments!'
  write(6,'(A)')   ' Example 1: fch_mo_copy a.fch b.fch      (Alpha MOs in a.fch -> Alpha MOs in b.fch)'
  write(6,'(A)')   ' Example 2: fch_mo_copy a.fch b.fch -aa  (eqv. to Example 1)'
  write(6,'(A)')   ' Example 3: fch_mo_copy a.fch b.fch -ab  (Alpha MOs in a.fch -> Beta  MOs in b.fch)'
  write(6,'(A)')   ' Example 4: fch_mo_copy a.fch b.fch -ba  (Beta  MOs in a.fch -> Alpha MOs in b.fch)'
  write(6,'(A)')   ' Example 5: fch_mo_copy a.fch b.fch -bb  (Beta  MOs in a.fch -> Beta  MOs in b.fch)'
  write(6,'(A,/)') ' Example 6: fch_mo_copy a.fch b.fch 5 10 (copy 5-10 Alpha MOs in a.fch -> b.fch)'
  stop
 end if

 call getarg(1, fchname1)
 call getarg(2, fchname2)

 if(i == 3) then
  call getarg(3,ab)
  if(ab(1:1) /= '-') then
   write(6,'(A)') "ERROR in subroutine fch_mo_copy: the 3rd argument must be '-aa'/'-ab'/'-ba'/'-bb'."
   stop
  end if
 else if(i == 4) then
  call getarg(3,str)
  read(str,*) idx1
  call getarg(4,str)
  read(str,*) idx2
  if(idx1<1 .or. idx2<1 .or. idx1>idx2) then
   write(6,'(A)') 'ERROR in subroutine fch_mo_copy: invalid indices idx1 and/or idx2.'
   write(6,'(2(A,I0))') 'idx1=', idx1, ', idx2=', idx2
   stop
  end if
 end if

 call fch_mo_copy(fchname1, fchname2, ab, idx1, idx2)
end program main

! copy Alpha/Beta MOs from a .fch file into Alpha/Beta MOs of another .fch file
subroutine fch_mo_copy(fchname1, fchname2, ab, idx1, idx2)
 implicit none
 integer :: nbf, nif, nbf2, nif2
 integer, intent(in) :: idx1, idx2
 real(kind=8), allocatable :: coeff1(:,:), coeff2(:,:)
 character(len=3), intent(in) :: ab
 character(len=240), intent(in) :: fchname1, fchname2

 call read_nbf_and_nif_from_fch(fchname1, nbf, nif)
 call read_nbf_and_nif_from_fch(fchname2, nbf2, nif2)

 if(nbf/=nbf2 .or. nif/=nif2) then
  write(6,'(/,A)') 'ERROR in subroutine fch_mo_copy: nbf and/or nif are not equ&
                   &al in two files:'
  write(6,'(A)') TRIM(fchname1)
  write(6,'(A)') TRIM(fchname2)
  stop
 end if

 allocate(coeff1(nbf,nif), coeff2(nbf,nif))
 call read_mo_from_fch(fchname1, nbf, nif, ab(2:2), coeff1)
 call read_mo_from_fch(fchname2, nbf, nif, ab(3:3), coeff2)

 if(idx1==0 .and. idx2==0) then
  coeff2 = coeff1
 else ! 0 < idx1 <= idx2
  coeff2(:,idx1:idx2) = coeff1(:,idx1:idx2)
 end if
 deallocate(coeff1)

 call write_mo_into_fch(fchname2, nbf, nif, ab(3:3), coeff2)
 deallocate(coeff2)
end subroutine fch_mo_copy

