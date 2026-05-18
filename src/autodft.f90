! first written by jxzou at 20260325: modern and advanced DFT calculations

! the hybrid functional part of DSD-PBEP86-D3(BJ):
! #p PBEP86 IOP(3/76=1000006900,3/77=0310003100,3/78=0440004400)

! the hybrid functional part of revDSD-PBEP86-D3(BJ):
! #p PBEP86 IOP(3/76=1000006900,3/77=0310003100,3/78=0429604296)

program main
 use mokit_version_info, only: version, date
 implicit none
 integer :: i
 character(len=29), parameter :: error_warn = 'ERROR in subroutine autodft: '
 character(len=240) :: fname

 i = iargc()
 if(i /= 1) then
  write(6,'(/,1X,A)') error_warn//'wrong command line arguments!'
  write(6,'(A)')   " Example: autodft h2o.gjf >h2o.out 2>&1 &"
  write(6,'(A,/)') ' See help: autodft -h'
  stop
 end if

 fname = ' '
 call getarg(1, fname)

 select case(TRIM(fname))
 case('-v', '-V', '--version')
  write(6,'(A)') 'AutoDFT '//TRIM(version)//' :: MOKIT, release date: '//TRIM(date)
  stop
 case('-h','-help','--help')
  write(6,'(/,A)') 'Usage: autodft [gjfname] > [outname]'
  write(6,'(A)')   "  Example: autodft h2o.gjf >h2o.out 2>&1 &"
  write(6,'(/,A)') 'Options:'
  write(6,'(A)')   '  -h, -help, --help: Print this message and exit.'
  write(6,'(A)')   '  -v, -V, --version: Print the version number of autodft and exit.'
  write(6,'(A)')   '  -t, --testprog: Print the path of programs detected by autodft and exit.'
  stop
 case('-t','--testprog')
  call check_mokit_root()
  call read_program_path()
  stop
 end select

 call require_gjf_suffix(fname)
 call require_file_exist(fname)
 call autodft(fname)
end program main

subroutine autodft(fname)
 implicit none
 character(len=240), intent(in) :: fname

end subroutine autodft

