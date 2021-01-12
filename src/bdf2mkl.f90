! written by jxzou at 20210111: combined two utilities - bdf2fch and fch2mkl

program main
 implicit none
 integer :: i
 integer, parameter :: iout = 6
 character(len=4) :: ab
 character(len=240) :: orbname, fchname
 logical :: alive

 i = iargc()
 if(.not. (i==2 .or. i==3)) then
  write(iout,'(/,1X,A)') 'ERROR in subroutine bdf2mkl: wrong command line arguments!'
  write(iout,'(1X,A)') 'Example 1 (for R(O)HF): bdf2mkl a.scforb a.fch'
  write(iout,'(1X,A)') 'Example 2 (for CAS)   : bdf2mkl a.casorb a.fch'
  write(iout,'(1X,A)') 'Example 3 (for UHF)   : bdf2mkl a.scforb a.fch -uhf'
  write(iout,'(1X,A,/)') 'Example 4 (for CAS NO): bdf2mkl a.casorb a.fch -no'
  stop
 end if

 ab = ' '; fchname = ' '
 call getarg(1,orbname)
 inquire(file=TRIM(orbname),exist=alive)
 if(.not. alive) then
  write(iout,'(A)') 'ERROR in subroutine bdf2mkl: file '//TRIM(orbname)//' does not exist.'
  stop
 end if

 call getarg(2,fchname)
 inquire(file=TRIM(fchname),exist=alive)
 if(.not. alive) then
  write(iout,'(A)') 'ERROR in subroutine bdf2mkl: file '//TRIM(fchname)//' does not exist.'
  stop
 end if

 if(i == 3) then
  call getarg(3, ab)
  ab = ADJUSTL(ab)
  if(ab/='-uhf' .and. ab/='-no') then
   write(iout,'(A)') "ERROR in subroutine bdf2mkl: the 3rd argument is&
                    & wrong! Only '-uhf' or '-no' is accepted."
   stop
  end if
 end if

 call bdf2mkl(orbname, fchname, ab)
 stop
end program main

! Step 1: call utility bdf2fch to transfer orbitals back into .fch(k) file
! Step 2: call utility fch2mkl to generate ORCA .inp and .mkl file
subroutine bdf2mkl(orbname, fchname, ab)
 implicit none
 integer :: i, system
 integer, parameter :: iout = 6
 character(len=4), intent(in) :: ab
 character(len=240), intent(in) :: orbname, fchname

 i = system('bdf2fch '//TRIM(orbname)//' '//TRIM(fchname)//' '//ab)
 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine bdf2mkl: failed to call utility bdf2fch.'
  write(iout,'(A)') 'Did you forget to compile utility bdf2fch?'
  stop
 end if

 if(TRIM(ab)=='-no' .or. LEN_TRIM(ab)==0) then
  i = system('fch2mkl '//TRIM(fchname))
 else
  i = system('fch2mkl '//TRIM(fchname)//' '//ab)
 end if

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine bdf2mkl: failed to call utility fch2mkl.'
  write(iout,'(A)') 'Did you forget to compile utility fch2mkl?'
  stop
 end if
 return
end subroutine bdf2mkl

