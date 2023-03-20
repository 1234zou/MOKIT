! written by jxzou at 20221230: automatic single reference calculations

! keywords information of single reference calculations (default values are set)
module sr_keyword
 use mr_keyword, only: gjfname, mem, nproc, method, basis, bgchg, cart, DKH2, &
  X2C, RI, F12, DLPNO, RIJK_bas, RIC_bas, F12_cabs, localm, readrhf, readuhf, &
  mo_rhf, hf_prog, hfonly, hf_fch, chgname, orca_path
 implicit none
 real(kind=8) :: ref_e = 0d0    ! reference wfn energy
 real(kind=8) :: corr_e = 0d0   ! total correlation energy
 real(kind=8) :: mp2_e = 0d0    ! MP2 total energy
 real(kind=8) :: ccd_e = 0d0    ! CCD total energy
 real(kind=8) :: ccsd_e = 0d0   ! CCSD total energy
 real(kind=8) :: ccsd_t_e = 0d0 ! CCSD(T) total energy
 real(kind=8) :: t1diag = 0d0   ! T1 diagnostic
 character(len=10) :: mp2_prog = 'orca'
 character(len=10) :: cc_prog = 'orca'
 character(len=10) :: adc_prog = 'pyscf'
 logical :: uhf_based = .true. ! UHF based MP2 or CC
 logical :: noRI = .false. ! turn on RI by default
 logical :: mp2 = .false.
 logical :: ccd = .false.
 logical :: ccsd = .false.
 logical :: ccsd_t = .false.
 logical :: cis = .false.
 logical :: adc = .false.
 logical :: eom = .false.
 logical :: iterative_t = .false. ! default DLPNO-CCSD(T0)

contains

subroutine read_sr_program_path()
 use mr_keyword, only: mokit_root, gau_path, molpro_path
 implicit none
 integer :: i
 integer(kind=4) :: hostnm
 character(len=8) :: hostname
 character(len=24) :: data_string
 character(len=240), external :: get_mokit_root 

 write(6,'(A)') '------ Output of AutoSR of MOKIT(Molecular Orbital Kit) ------'
 write(6,'(A)') '       GitLab page: https://gitlab.com/jxzou/mokit'
 write(6,'(A)') '     Documentation: https://jeanwsr.gitlab.io/mokit-doc-mdbook'
 write(6,'(A)') '           Version: 1.2.5rc19 (2023-Mar-20)'
 write(6,'(A)') '       How to cite: see README.md or $MOKIT_ROOT/doc/'

 hostname = ' '
 data_string = ' '
 i = hostnm(hostname)
 call fdate(data_string)
 write(6,'(/,A)') 'HOST '//TRIM(hostname)//', '//TRIM(data_string)

 write(6,'(/,A)') 'Read program paths from environment variables:'
 mokit_root = get_mokit_root()
 write(6,'(A)') 'MOKIT_ROOT  = '//TRIM(mokit_root)

 call get_gau_path(gau_path)
 call get_orca_path(orca_path)
 call get_molpro_path(molpro_path)

 write(6,'(A)') 'gau_path    = '//TRIM(gau_path)
 write(6,'(A)') 'orca_path   = '//TRIM(orca_path)
 write(6,'(A)') 'molpro_path = '//TRIM(molpro_path)
end subroutine read_sr_program_path

! Parse keywords of single reference calculations
subroutine parse_sr_keyword()
 use mr_keyword, only: basname, check_readfch
 implicit none
 integer :: i, j, k, ifail, fid
 character(len=24) :: method0 = ' '
 character(len=240) :: buf = ' '
 character(len=1000) :: longbuf = ' '

 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
 call read_mem_nproc_route(fid, mem, nproc, buf)

 i = index(buf,'/')
 if(i == 0) then
  write(6,'(A)') "ERROR in subroutine parse_sr_keyword: no '/' symbol detected&
                 &in keyword line."
  write(6,'(A)') "The method and basis set must be specified via '/' symbol,&
                 & e.g. CCSD(T)/cc-pVTZ."
  close(fid)
  stop
 end if

 j = index(buf(1:i-1),' ', back=.true.)
 if(j == 0) then
  write(6,'(A)') 'ERROR in subroutine parse_sr_keyword: syntax error detected in'
  write(6,'(A)') "the current line '"//TRIM(buf)//"'"
  stop
 end if
 method0 = buf(j+1:i-1)

 if(i /= 0) then
  method = method0(1:i-1)

  select case(TRIM(method))
  case('ccd','ccsd','ccsd(t)','ccsd(t)-f12','dlpno-ccsd','dlpno-ccsd(t)',&
       'dlpno-ccsd(t0)','dlpno-ccsd(t1)')
  case default
   write(6,'(/,A)') 'ERROR in subroutine parse_sr_keyword: unsupported method '&
                   //TRIM(method)
   close(fid)
   stop
  end select
 else ! i = 0
  method = TRIM(method0)
 end if

 RI = .true. ! turn on RI by default
 select case(TRIM(method))
 case('mp2')
  mp2 = .true.
 case('ccd')
  ccd = .true.
 case('ccsd')
  ccsd = .true.
 case('ccsd(t)')
  ccsd_t = .true.
 case('ccsd(t)-f12')
  ccsd_t = .true.; F12 = .true.
 case('dlpno-mp2')
  mp2 = .true.; DLPNO = .true.
 case('dlpno-ccsd')
  ccsd = .true.; DLPNO = .true.
 case('dlpno-ccsd(t)','dlpno-ccsd(t0)')
  ccsd_t = .true.; DLPNO = .true.
 case('dlpno-ccsd(t1)')
  ccsd_t = .true.; DLPNO = .true.; iterative_t = .true.
 case default
  write(6,'(/,A)') 'ERROR in subroutine parse_sr_keyword: unrecognized method='&
                  //TRIM(method)
  stop
 end select

 i = index(buf,'/')
 j = i - 1 + index(buf(i+1:),' ')
 if(j == 0) j = LEN_TRIM(buf)
 basis = buf(i+1:j)
 write(6,'(/,2(A,I4))',advance='no') 'memory =', mem, 'GB, nproc =', nproc
 write(6,'(A)') ', method/basis = '//TRIM(method)//'/'//TRIM(basis)

 if(basis(1:5) == 'def2-') then
  write(6,'(A)') "ERROR in subroutine parse_sr_keyword: 'def2-' prefix detec&
                 &ted in given basis set."
  write(6,'(A)') 'Basis set in Gaussian syntax should be like def2TZVP, not &
                 &def2-TZVP.'
  stop
 end if

 if(basis(1:3) == 'gen') then
  close(fid)
  if(.not. check_readfch(gjfname)) then
   call record_gen_basis_in_gjf(gjfname, basname, .true.)
  end if

  open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
  do while(.true.)
   read(fid,'(A)') buf
   if(buf(1:1) == '#') exit
  end do ! for while
 end if

 read(fid,'(A)') buf ! skip a blank line
 read(fid,'(A)') buf ! Title Card, the 1st line of keywords
 call lower(buf)

 if(buf(1:6) /= 'mokit{') then
  write(6,'(A)') "ERROR in subroutine parse_sr_keyword: 'mokit{' not detected&
                 & in file "//TRIM(gjfname)
  write(6,'(A)') "Syntax error. You must put 'mokit{' in leading position of &
                 &the Title Card line."
  stop
 end if

 j = index(buf,'}')
 if(j==7 .or. (j>7 .and. LEN_TRIM(buf(7:j-1))==0)) then ! mokit{}
  close(fid)
  return
 else if(j == 0) then ! keywords written in more than 1 line
  j = LEN_TRIM(buf) + 1
 else ! j > 0
  j = LEN_TRIM(buf)
 end if

 longbuf(1:j-7) = buf(7:j-1) ! some keywords specified
 k = j - 6 ! the beginning index for next keyword in longbuf

 if(index(buf,'}') == 0) then ! keywords are written in at least two lines
  do while(.true.)
   read(fid,'(A)',iostat=ifail) buf
   if(ifail /= 0) exit

   call lower(buf)
   i = LEN_TRIM(buf)
   if(index(buf,'}') > 0) i = i - 1
   longbuf(k:k+i-1) = buf(1:i)
   k = k + i

   if(index(buf,'}') > 0) exit
  end do ! for while

  if(ifail /= 0) then
   write(6,'(A)') 'ERROR in subroutine parse_sr_keyword: end-of-file detected.'
   write(6,'(A)') 'The provided .gjf file may be incomplete.'
   close(fid)
   stop
  end if
 end if

 close(fid)
 ! now all keywords are stored in longbuf

 write(6,'(/,A)') 'Keywords in MOKIT{} are merged and shown as follows:'
 write(6,'(A)') TRIM(longbuf)

 do while(.true.)
  i = index(longbuf,',')
  if(i == 0) i = LEN_TRIM(longbuf) + 1

  j = index(longbuf(1:i-1),'=')
  if(j == 0) j = i   ! in case this keyword has no value assigned

  select case(longbuf(1:j-1))
  case('readrhf')
   mo_rhf = .true.
   readrhf = .true.
   read(longbuf(j+1:i-1),*) hf_fch
  case('readuhf')
   readuhf = .true.
   read(longbuf(j+1:i-1),*) hf_fch
  case('cart')      ! use Cartesian functions
   cart = .true.
  case('dkh2')
   DKH2 = .true.
  case('x2c')
   X2C = .true.
  case('localm')    ! localization method
   read(longbuf(j+1:i-1),*) localm
   localm = ADJUSTL(localm)
  case('hf_prog')
   read(longbuf(j+1:i-1),*) hf_prog
  case('mp2_prog')
   read(longbuf(j+1:i-1),*) mp2_prog
  case('cc_prog')
   read(longbuf(j+1:i-1),*) cc_prog
  case('adc_prog')
   read(longbuf(j+1:i-1),*) adc_prog
  case('charge')
   bgchg = .true.
  case('nori')
   noRI = .true.; RI = .false.
  case('rijk_bas')
   read(longbuf(j+1:i-1),*) RIJK_bas
  case('ric_bas')
   read(longbuf(j+1:i-1),*) RIC_bas
  case('f12')
   F12 = .true.; RI = .true.
  case('f12_cabs')
   read(longbuf(j+1:i-1),*) F12_cabs
  case('dlpno')
   DLPNO = .true.; RI = .true.
  case('hfonly')
   HFonly = .true.
  case default
   write(6,'(/,A)') "ERROR in subroutine parse_sr_keyword: keyword '"//longbuf(1:j-1)&
                   //"' not recognized in {}."
   stop
  end select

  ! delete the keyword which has been specified
  longbuf(1:i) = ' '
  longbuf = ADJUSTL(longbuf)
  if(LEN_TRIM(longbuf) == 0) exit
 end do ! for while
end subroutine parse_sr_keyword

! check the compatiblity of keywords in single reference calculations
subroutine check_sr_kywd_compatible()
 implicit none

 if(ccd .and. TRIM(cc_prog)=='molpro') then
  write(6,'(/,A)') 'ERROR in subroutine check_sr_kywd_compatible: CC_prog=Molpr&
                   &o is incompatible with'
  write(6,'(A)') 'the CCD method, you may want to try CC_prog=ORCA.'
  stop
 end if

 if(RI .eqv. noRI) then
  write(6,'(A)') 'ERROR in subroutine check_sr_kywd_compatible: RI=noRI. Confus&
                 &ed keywords.'
  stop
 end if
 if(RI) call determine_auxbas(basis, RIJK_bas, .true., RIC_bas, F12, F12_cabs)
 call prt_sr_strategy()
end subroutine check_sr_kywd_compatible

! print the internal variables of single reference calculations
subroutine prt_sr_strategy()
 implicit none

 if(DLPNO .and. ccsd_t .and. (.not.iterative_t)) then
  write(6,'(/,A)') 'Remark: DLPNO-CCSD(T) is DLPNO-CCSD(T0) by default in ORCA.&
                   & If you want higher'
  write(6,'(A)') 'accuracy, it is recommended to use the DLPNO-CCSD(T1) method.'
 end if

 if(.not. RI) then
  write(6,'(/,A)') 'Remark: RI status is off currently. It is strongly recomme&
                   &nded to turn on'
  write(6,'(A)') 'the RI technique unless you want to do some tests with no RI.'
 end if

 write(6,'(/,A)') 'Internal variables:'
 write(6,'(5(A,L1,3X))') 'MP2     = ', mp2, 'CCD     = ', ccd, 'CCSD    = ', &
                          ccsd, 'CCSD(T) = ', ccsd_t, 'ADC     = ', adc
 write(6,'(5(A,L1,3X))') 'EOM     = ', eom, 'RI      = ',  RI, 'F12     = ', &
                          F12,  'DLPNO   = ',  DLPNO, '(T1)    = ', iterative_t
 write(6,'(A)') 'RIJK_bas='//TRIM(RIJK_bas)//'  RIC_bas='//TRIM(RIC_bas)//&
                ' F12_cabs='//TRIM(F12_cabs)
end subroutine prt_sr_strategy

end module sr_keyword

program main
 use sr_keyword, only: read_sr_program_path
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
  write(6,'(A)') 'AutoSR 1.2.5rc19 :: MOKIT, release date: 2023-Mar-20'
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
  write(6,'(A)')   '  DLPNO-CCSD(T), DLPNO-CCSD(T0), DLPNO-CCSD(T1), EOM-CCSD,&
                   & EOM-SF-CCSD, EOM-IP-CCSD,'
  write(6,'(A)')   '  EOM-DIP-CCSD, EOM-EA-CCSD, CC2, CC3, ADC(2), ADC(3), CIS&
                   &(D), ROCIS, XCIS, SF-XCIS,'
  write(6,'(A)')   '  SA-SF-CIS'
  write(6,'(/,A)') 'Frequently used keywords in MOKIT{}:'
  write(6,'(A)')   '     HF_prog=Gaussian/PySCF/ORCA/PSI4'
  write(6,'(A)')   '    MP2_prog=Molpro/ORCA/Gaussian/GAMESS/PySCF/Dalton/QChem'
  write(6,'(A)')   '     CC_prog=Molpro/CFOUR/ORCA/PSI4/Gaussian/GAMESS/PySCF/Dalton/QChem'
  write(6,'(A)')   '    CIS_prog=Molpro/ORCA/Gaussian/QChem'
  write(6,'(A)')   '    ADC_prog=Molpro/ORCA/PySCF/QChem'
  write(6,'(A,/)') '  EOMCC_prog=Molpro/CFOUR/ORCA/Gaussian/GAMESS/PySCF/QChem'
  stop
 case('-t','--testprog')
  call read_sr_program_path()
  stop
 end select

 i = index(fname, '.gjf', back=.true.)
 j = index(fname, '.fch', back=.true.)
 if(i>0 .and. j>0) then
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
 use sr_keyword, only: gjfname, read_sr_program_path, parse_sr_keyword, &
  check_sr_kywd_compatible
 implicit none
 character(len=24) :: data_string
 character(len=240), intent(in) :: fname

 gjfname = fname
 ! read paths of various programs from environment variables
 call read_sr_program_path()
 call parse_sr_keyword()
 call check_sr_kywd_compatible()

 call do_hf(.false.) ! RHF, ROHF, UHF
 call do_mp2()   ! (DLPNO-, RI-)MP2
 call do_cc()    ! (DLPNO-)CCD, CCSD, CCSD(T), etc
 call do_cis()   ! CIS(D), ROCISD, XCIS, etc
 call do_adc()   ! ADC(2), ADC(3)
 call do_eomcc() ! EOM-CCSD, EOM-SF-CCSD

 call fdate(data_string)
 write(6,'(/,A)') 'Normal termination of AutoSR at '//TRIM(data_string)
end subroutine autosr

subroutine do_mp2()
 use sr_keyword, only: mp2
 implicit none
 character(len=24) :: data_string = ' '

 if(.not. mp2) return
 write(6,'(//,A)') 'Enter subroutine do_mp2...'

 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine do_mp2 at '//TRIM(data_string)
end subroutine do_mp2

subroutine do_cc()
 use sr_keyword, only: mem, nproc, ccd, ccsd, ccsd_t, bgchg, hf_fch, chgname, &
  cc_prog, orca_path, method, ccsd_e, ccsd_t_e, ref_e, corr_e, t1diag
 use util_wrapper, only: fch2mkl_wrap, mkl2gbw, fch2com_wrap
 implicit none
 integer :: i, RENAME, SYSTEM
 real(kind=8) :: e = 0d0
 character(len=15) :: method0 = ' '
 character(len=24) :: data_string = ' '
 character(len=240) :: old_inp, inpname, mklname, outname

 if(.not. (ccd .or. ccsd .or. ccsd_t)) return
 write(6,'(//,A)') 'Enter subroutine do_cc...'
 write(6,'(A)') 'CC using program '//TRIM(cc_prog)

 method0 = method
 call upper(method0)
 i = index(hf_fch, '.fch', back=.true.)
 outname = hf_fch(1:i-1)//'_CC.out'

 select case(TRIM(cc_prog))
 case('orca')
  inpname = hf_fch(1:i-1)//'_CC.inp'
  old_inp = hf_fch(1:i-1)//'_o.inp'
  mklname = hf_fch(1:i-1)//'_CC.mkl'
  call fch2mkl_wrap(hf_fch, mklname)
  i = RENAME(TRIM(old_inp), TRIM(inpname))
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  ! if bgchg = .True., .inp and .mkl file will be updated
  call mkl2gbw(mklname)
  call delete_file(mklname)
  call prt_cc_orca_inp(inpname)
  call submit_orca_job(orca_path, inpname)
  call read_cc_e_from_orca_out(outname, (.not.ccd), t1diag, ref_e, e)
 case('molpro')
  inpname = hf_fch(1:i-1)//'_CC.com'
  call fch2com_wrap(hf_fch, inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call prt_cc_molpro_inp(inpname)
  call submit_molpro_job(inpname, mem, nproc)
  call read_cc_e_from_molpro_out(outname, (.not.ccd), t1diag, ref_e, e)
 case default
  write(6,'(A)') 'ERROR in subroutine do_cc: invalid CC_prog='//TRIM(cc_prog)
  stop
 end select

 if(ccsd) then
  ccsd_e = e
 else
  ccsd_t_e = e
 end if
 corr_e = e - ref_e
 write(6,'(/,A,F18.8,A)')'E(ref)  = ', ref_e, ' a.u.'
 write(6,'(A,F18.8,A)')  'E(corr) = ', corr_e, ' a.u.'
 write(6,'(A,F18.8,A)')  'E('//TRIM(method0)//') = ', e, ' a.u.'
 if(.not. ccd) write(6,'(A,F10.5,A)')  'T1diag = ', t1diag

 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine do_cc at '//TRIM(data_string)
end subroutine do_cc

subroutine do_cis()
 use sr_keyword, only: cis
 implicit none
 character(len=24) :: data_string = ' '

 if(.not. cis) return
 write(6,'(//,A)') 'Enter subroutine do_cis...'

 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine do_cis at '//TRIM(data_string)
end subroutine do_cis

subroutine do_adc()
 use sr_keyword, only: adc
 implicit none
 character(len=24) :: data_string = ' '

 if(.not. adc) return
 write(6,'(//,A)') 'Enter subroutine do_adc...'

 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine do_adc at '//TRIM(data_string)
end subroutine do_adc

subroutine do_eomcc()
 use sr_keyword, only: eom
 implicit none
 character(len=24) :: data_string = ' '

 if(.not. eom) return
 write(6,'(//,A)') 'Enter subroutine do_eomcc...'

 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine do_eomcc at '//TRIM(data_string)
end subroutine do_eomcc

! add CC keywords into a ORCA input file
subroutine prt_cc_orca_inp(inpname)
 use sr_keyword, only: mem, nproc, DLPNO, F12, RIJK_bas, RIC_bas, F12_cabs, &
  mo_rhf, ccd, ccsd, ccsd_t, iterative_t
 use mol, only: mult
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 i = index(inpname, '.inp', back=.true.)
 inpname1 = inpname(1:i-1)//'.t'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:6) == '%coord') exit
 end do ! for while

 write(fid1,'(A,I0,A)') '%pal nprocs ',nproc,' end'
 i = FLOOR(0.8d0*DBLE(mem*1000)/DBLE(nproc))
 write(fid1,'(A,I0)') '%maxcore ', i ! MB

 if(mo_rhf) then
  if(mult > 1) then
   write(fid1,'(A)',advance='no') '! RO'
  else
   write(fid1,'(A)',advance='no') '! R'
  end if
 else
  write(fid1,'(A)',advance='no') '! U'
 end if
 write(fid1,'(A)',advance='no') 'HF TightSCF RIJK '//TRIM(RIJK_bas)//' '//&
                                TRIM(RIC_bas)
 if(F12) write(fid1,'(A)',advance='no') ' '//TRIM(F12_cabs)

 if(DLPNO) then
  write(fid1,'(A)',advance='no') ' TightPNO DLPNO-CCSD'
  if(iterative_t) then
   write(fid1,'(A)',advance='no') '(T1)'
  else if(ccsd_t) then
   write(fid1,'(A)',advance='no') '(T)'
  end if
 else
  if(ccsd_t) then
   write(fid1,'(A)',advance='no') ' CCSD(T)'
  else if(ccsd) then
   write(fid1,'(A)',advance='no') ' CCSD'
  else if(ccd) then
   write(fid1,'(A)',advance='no') ' CCD'
  end if
 end if
 if(F12) then
  if(DLPNO) then
   write(fid1,'(A)') '-F12'
  else
   write(fid1,'(A)') '-F12/RI'
  end if
 else
  write(fid1,'(/)',advance='no')
 end if

 write(fid1,'(A)') TRIM(buf)
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine prt_cc_orca_inp

subroutine prt_cc_molpro_inp(inpname)
 use sr_keyword, only: ccsd, ccsd_t, RI, F12, mo_rhf, uhf_based, basis, RIJK_bas,&
  RIC_bas
 use mol, only: mult
 implicit none
 integer :: i, fid
 character(len=21) :: MP2FIT_bas, RIJK_bas1
 character(len=240), intent(in) :: inpname

 call auxbas_convert(RIJK_bas, RIJK_bas1, 2)
 i = index(RIJK_bas, '/')
 RIC_bas = RIJK_bas(1:i)//'optri'
 MP2FIT_bas = RIJK_bas(1:i)//'mp2fit'
 open(newunit=fid,file=TRIM(inpname),status='old',position='append')
 BACKSPACE(fid)

 if(RI) then
  ! do 1 cycle of R(O)HF, let molpro recognize auxiliary basis sets
  write(fid,'(A,/,A)') 'basis='//TRIM(basis), 'DF-HF'
  write(fid,'(A)') 'basis={'
  write(fid,'(A,/,A)') 'set,df', 'default,'//TRIM(MP2FIT_bas)
  write(fid,'(A,/,A)') 'set,jk', 'default,'//TRIM(RIJK_bas1)
  if(F12) then
   write(fid,'(A,/,A)') 'set,ri', 'default,'//TRIM(RIC_bas)
  end if
  write(fid,'(A)') '}'
  write(fid,'(A)',advance='no') 'DF-'
  if(mult /= 1) then
   if(mo_rhf) then
    if(uhf_based) then  ! ROHF-UCCSD
     write(fid,'(A)',advance='no') 'U'
    else                ! ROHF-ROCCSD
     write(fid,'(A)',advance='no') 'R'
    end if
   else
    write(6,'(A)') 'ERROR in subroutine prt_cc_molpro_inp: UHF-based CC is not&
                   & supported in Molpro.'
    close(fid)
    stop
   end if
  end if
  write(fid,'(A)',advance='no') 'CCSD'
  if(ccsd_t) write(fid,'(A)',advance='no') '(T)'
  if(F12) write(fid,'(A)',advance='no') '-F12,ri_basis=ri'
  write(fid,'(A)') ',df_basis=df,df_basis_exch=jk'

 else ! RI is turned off
  if(ccsd_t) then
   write(fid,'(A)',advance='no') 'CCSD(T)'
  else if(ccsd) then
   write(fid,'(A)',advance='no') 'CCSD'
  end if
  if(F12) then
   write(fid,'(A)') '-F12'
  else
   write(fid,'(/)',advance='no')
  end if
 end if

 close(fid)
end subroutine prt_cc_molpro_inp

subroutine read_cc_e_from_orca_out(outname, read_t1diag, t1diag, ref_e, tot_e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: t1diag, ref_e, tot_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: read_t1diag

 t1diag = 0d0; ref_e = 0d0; tot_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:10) == 'Total Ener') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_orca_out: 'Total Ener' n&
                   &ot found in file "//TRIM(outname)
  close(fid)
  stop
 end if

 i = index(buf, ':')
 read(buf(i+1:),*) ref_e

 if(read_t1diag) then
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:10) == 'T1 diagnos') exit
  end do ! for while
  read(buf(47:),*) t1diag
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:10) == 'FINAL SING') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_orca_out: 'FINAL SING' n&
                   &ot found in file "//TRIM(outname)
  stop
 end if

 read(buf(26:),*) tot_e
end subroutine read_cc_e_from_orca_out

subroutine read_cc_e_from_molpro_out(outname, read_t1diag, t1diag, ref_e, tot_e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: t1diag, ref_e, tot_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: read_t1diag

 t1diag = 0d0; ref_e = 0d0; tot_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(2:9) == 'Checking') then
   i = -1
   exit
  end if
  if(index(buf,'HF-SCF') > 0) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_orca_out: 'Total Ener' n&
                   &ot found in file "//TRIM(outname)
  close(fid)
  stop
 end if

 read(fid,*) tot_e, ref_e

 if(read_t1diag) then
  do while(.true.)
   BACKSPACE(fid,iostat=i)
   if(i /= 0) exit
   BACKSPACE(fid,iostat=i)
   if(i /= 0) exit
   read(fid,'(A)') buf
   if(buf(2:9) == 'Checking') then
    i = -1
    exit
   end if
   if(index(buf,'T1 diag') > 0) exit
  end do ! for while
  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_orca_out: 'T1 diag' not&
                   & found in file "//TRIM(outname)
   close(fid)
   stop
  end if
  i = index(buf, 'T1 diagnostic:')
  read(buf(i+14:),*) t1diag
 end if

 close(fid)
end subroutine read_cc_e_from_molpro_out

