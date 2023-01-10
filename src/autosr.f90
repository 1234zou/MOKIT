! written by jxzou at 20221230: automatic single reference calculations

program main
 use mr_keyword, only: read_program_path
 implicit none
 integer :: i, j
 character(len=240) :: fname = ' '

 i = iargc()
 if(i /= 1) then
  write(6,'(/,A)') ' ERROR in subroutine autosr: wrong command line argument!'
  write(6,'(A)')   " Example 1 (in bash): autosr a.gjf >& a.out &"
  write(6,'(A)')   " Example 2 (in dash): autosr a.gjf >a.out 2>&1 &"
  write(6,'(A,/)') ' See help: autosr -h'
  stop
 end if

 call getarg(1, fname)

 select case(TRIM(fname))
 case('-v', '-V', '--version')
  write(6,'(A)') 'AutoSR 1.2.5 :: MOKIT, release date: 2023-Jan-10'
  stop
 case('-h','-help','--help')
  write(6,'(/,A)') "Usage: autosr [gjfname] >& [outname]"
  write(6,'(A)')   "  Example 1 (in bash): autosr a.gjf >& a.out &"
  write(6,'(A)')   "  Example 2 (in dash): autosr a.gjf >a.out 2>&1 &"
  write(6,'(/,A)') 'Options:'
  write(6,'(A)')   '  -h, -help, --help: Print this message and exit.'
  write(6,'(A)')   '  -v, -V, --version: Print the version number of autosr and exit.'
  write(6,'(A)')   '  -t, --testprog: Print the path of programs detected by autosr and exit.'
  write(6,'(/,A)') 'Methods(#p ...):'
  write(6,'(A)')   '  MP2, RI-MP2, CCD, CCSD, CCSD(T), CCSD-F12, CCSD(T)-F12, &
                   &DLPNO-MP2, DLPNO-CCSD,'
  write(6,'(A)')   '  DLPNO-CCSD(T), DLPNO-CCSD(T1), EOM-CCSD, EOM-SF-CCSD, EOM&
                   &-IP-CCSD, EOM-DIP-CCSD,'
  write(6,'(A)')   '  EOM-EA-CCSD, CC2, CC3, ADC(2), ADC(3), CIS(D), ROCIS, XCI&
                   &S, SF-XCIS, SA-SF-CIS'
  write(6,'(/,A)') 'Frequently used keywords in MOKIT{}:'
  write(6,'(A)')   '     HF_prog=Gaussian/PySCF/ORCA/PSI4'
  write(6,'(A)')   '    MP2_prog=Molpro/ORCA/Gaussian/GAMESS/PySCF/OpenMolcas/Dalton/QChem'
  write(6,'(A)')   '     CC_prog=Molpro/ORCA/Gaussian/GAMESS/PySCF/OpenMolcas/Dalton/QChem'
  write(6,'(A)')   '    CIS_prog=Molpro/ORCA/Gaussian/QChem'
  write(6,'(A)')   '    ADC_prog=Molpro/ORCA/PySCF/QChem'
  write(6,'(A,/)') '  EOMCC_prog=Molpro/ORCA/Gaussian/GAMESS/PySCF/QChem'
  stop
 case('-t','--testprog')
  call read_program_path()
  stop
 end select

 i = index(fname, '.gjf', back=.true.)
 j = index(fname, '.fch', back=.true.)
 if(i/=0 .and. j/=0) then
  write(6,'(/,A)') "ERROR in subroutine autosr: both '.gjf' and '.fch' keys&
                   & detected in filename "//TRIM(fname)//'.'
  write(6,'(A)') "Better to use a filename only with suffix '.gjf'."
  stop
 else if(i == 0) then
  write(6,'(/,A)') "ERROR in subroutine autosr: '.gjf' key not found in&
                   & filename "//TRIM(fname)
  stop
 end if

 call require_file_exist(fname)
 call autosr(fname)
end program main

! automatically do single reference calculations
subroutine autosr(fname)
 use mr_keyword, only: gjfname, read_program_path, parse_keyword, check_kywd_compatible
 implicit none
 character(len=24) :: data_string
 character(len=240), intent(in) :: fname

 gjfname = fname
 ! read paths of various programs from environment variables
 call read_program_path()
 call parse_keyword()
 call check_kywd_compatible()

 call do_hf()  ! RHF, ROHF, UHF
 call do_mp2() ! (RI-)MP2
 call do_cc()  ! (DLPNO-)CCD, CCSD, CCSD(T), etc
 call do_cis() ! CIS(D), ROCISD, XCIS, etc
 call do_adc() ! ADC(2), ADC(3)
 call do_eomcc() ! EOM-CCSD, EOM-SF-CCSD

 call fdate(data_string)
 write(6,'(/,A)') 'Normal termination of AutoSR at '//TRIM(data_string)
end subroutine autosr

subroutine do_mp2()
 implicit none
end subroutine do_mp2

subroutine do_cc()
 implicit none
end subroutine do_cc

subroutine do_cis()
 implicit none
end subroutine do_cis

subroutine do_adc()
 implicit none
end subroutine do_adc

subroutine do_eomcc()
 implicit none
end subroutine do_eomcc

