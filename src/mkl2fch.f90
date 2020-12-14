! written by jxzou at 20200627: transfer MOs from .mkl to .fch(k)
! updated by jxzou at 20201213: read NOONs from .mkl and save to .fch

! The 'Shell types' array in Gaussian .fch file:
!
!   Spherical     |     Cartesian
! -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5
!  H  G  F  D  L  S  P  D  F  G  H
!
! 'L' is 'SP' in Pople-type basis sets

program main
 implicit none
 integer :: i
 integer, parameter :: iout = 6
 character(len=4) :: str = ' '
 character(len=240) :: mklname, fchname
 logical :: uhf, read_on

 i = iargc()
 if(i<1 .or. i>3) then
  write(iout,'(/,A)')  ' ERROR in subroutine mkl2fch: wrong command line arguments.'
  write(iout,'(A)')  ' Example 1 (R(O)HF, CAS): mkl2fch a.mkl a.fch'
  write(iout,'(A)')  ' Example 2 (CAS NO)     : mkl2fch a.mkl a.fch -no'
  write(iout,'(A,/)')  ' Example 3 (UHF)        : mkl2fch a.mkl a.fch -uhf'
  stop
 end if

 mklname = ' '; fchname = ' '
 call getarg(1, mklname)
 call getarg(2, fchname)

 uhf = .false.; read_on = .false.

 if(i == 3) then
  call getarg(3, str)
  str = ADJUSTL(str)

  select case(TRIM(str))
  case('-uhf')
   uhf = .true.
  case('-no')
   read_on = .true.
  case default
   write(iout,'(A)') 'ERROR in subroutine mkl2fch: wrong command line arguments.'
   write(iout,'(A)') "The 3rd input parameter can only be '-uhf' or '-no'. But&
                    & you specify '"//TRIM(str)//"'"
   stop
  end select
 end if

 call mkl2fch(mklname, fchname, uhf, read_on)
 stop
end program main

! convert .fch(k) file (Gaussian) to .mkl file (Molekel, ORCA)
subroutine mkl2fch(mklname, fchname, uhf, read_on)
 use fch_content
 implicit none
 integer :: i, j, k, nf3mark, ng3mark, nh3mark
 integer, allocatable :: f3_mark(:), g3_mark(:), h3_mark(:)
 real(kind=8), allocatable :: noon(:)
 character(len=240), intent(in) :: mklname, fchname
 logical, intent(in) :: uhf, read_on

 i = INDEX(fchname,'.fch',back=.true.)
 if(i == 0) then
  write(iout,'(A)') "ERROR in subroutine mkl2fch: input filename does not&
                   & contain '.fch' suffix!"
  write(iout,'(A)') 'fchname='//TRIM(fchname)
  stop
 end if

 call read_fch(fchname, uhf) ! read content in .fch(k) file

 ! check if any Cartesian functions
 if( ANY(shell_type > 1) ) then
  write(iout,'(A)') 'ERROR in subroutine mkl2fch: Cartesian functions detected&
                   & in file '//TRIM(fchname)//'.'
  write(iout,'(A)') "ORCA supports only spherical functions. You need to add&
                  & '5D 7F' keywords in Gaussian."
  stop
 else if( ANY(shell_type < -5) ) then
  write(iout,'(A)') 'ERROR in subroutine mkl2fch: angular momentum too high! not supported.'
  stop
 end if
 ! check done

 ! find F+3, G+3 and H+3 functions, multiply them by -1
 nf3mark = 0; ng3mark = 0; nh3mark = 0
 allocate(f3_mark(nbf), source=0)
 allocate(g3_mark(nbf), source=0)
 allocate(h3_mark(nbf), source=0)
 k = 0
 do i = 1, ncontr, 1
  select case(shell_type(i))
  case(0) ! S
   k = k + 1
  case(1) ! P
   k = k + 3
  case(-1) ! L
   k = k + 4
  case(-2) ! D
   k = k + 5
  case(-3) ! F
   k = k + 7
   nf3mark = nf3mark + 1
   f3_mark(nf3mark) = k - 1
  case(-4) ! G
   k = k + 9
   ng3mark = ng3mark + 1
   g3_mark(ng3mark) = k - 3
  case(-5) ! H
   k = k + 11
   nh3mark = nh3mark + 1
   h3_mark(nh3mark) = k - 5
  end select
 end do ! for i

 if(k /= nbf) then
  write(iout,'(A)') 'ERROR in subroutine mkl2fch: k /= nbf!'
  write(iout,'(2(A,I0))') 'k=', k, ', nbf=', nbf
  stop
 end if

 ! read Alpha and/or Beta MOs from .fch(k) file
 call read_mo_from_mkl(mklname, nbf, nif, 'a', alpha_coeff)
 if(uhf) call read_mo_from_mkl(mklname, nbf, nif, 'b', beta_coeff)

 do i = 1, nf3mark, 1
  alpha_coeff(f3_mark(i),:) = -alpha_coeff(f3_mark(i),:)
  alpha_coeff(f3_mark(i)+1,:) = -alpha_coeff(f3_mark(i)+1,:)
 end do ! for i
 do i = 1, ng3mark, 1
  alpha_coeff(g3_mark(i):g3_mark(i)+3,:) = -alpha_coeff(g3_mark(i):g3_mark(i)+3,:)
 end do ! for i
 do i = 1, nh3mark, 1
  alpha_coeff(h3_mark(i):h3_mark(i)+3,:) = -alpha_coeff(h3_mark(i):h3_mark(i)+3,:)
 end do ! for i

 if(uhf) then
  do i = 1, nf3mark, 1
   beta_coeff(f3_mark(i),:) = -beta_coeff(f3_mark(i),:)
   beta_coeff(f3_mark(i)+1,:) = -beta_coeff(f3_mark(i)+1,:)
  end do ! for i
  do i = 1, ng3mark, 1
   beta_coeff(g3_mark(i):g3_mark(i)+3,:) = -beta_coeff(g3_mark(i):g3_mark(i)+3,:)
  end do ! for i
  do i = 1, nh3mark, 1
   beta_coeff(h3_mark(i):h3_mark(i)+3,:) = -beta_coeff(h3_mark(i):h3_mark(i)+3,:)
  end do ! for i
 end if
 deallocate(f3_mark, g3_mark, h3_mark)

 call write_mo_into_fch(fchname, nbf, nif, 'a', alpha_coeff)
 deallocate(alpha_coeff)

 if(uhf) then
  call write_mo_into_fch(fchname, nbf, nif, 'b', beta_coeff)
  deallocate(beta_coeff)
 end if

 if(read_on) then
  allocate(noon(nif))
  call read_on_from_mkl(mklname, nif, noon)
  call write_eigenvalues_to_fch(fchname, nif, 'a', noon, .true.)
  deallocate(noon)
 end if
 return
end subroutine mkl2fch

