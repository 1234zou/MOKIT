! written by jxzou at 20221230: automatic single reference calculations

! keywords information of single reference calculations (default values are set)
module sr_keyword
 use mr_keyword, only: gjfname, mem, nproc, method, basis, bgchg, cart, force, &
  DKH2, X2C, dkh2_or_x2c, RI, F12, DLPNO, RIJK_bas, RIC_bas, F12_cabs, localm, &
  readrhf, readuhf, mo_rhf, hf_prog, hfonly, hf_fch, chgname, gau_path, &
  orca_path, psi4_path, gms_path, gms_scr_path, check_gms_path
 use mol, only: chem_core, ecp_core
 implicit none
 integer :: core_wish = 0 ! the number of frozen core orbitals the user wishes
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
 logical :: uhf_based = .false. ! UHF based MP2 or CC
 logical :: noRI = .false. ! turn on RI by default
 logical :: mp2 = .false.
 logical :: ccd = .false.
 logical :: ccsd = .false.
 logical :: ccsd_t = .false.
 logical :: cc_enabled = .false. ! (ccd .or. ccsd .or. ccsd_t)
 logical :: cis = .false.
 logical :: adc = .false.
 logical :: eom = .false.
 logical :: iterative_t = .false. ! default DLPNO-CCSD(T0)
 logical :: mp2_force = .false.
 logical :: ccd_force = .false.
 logical :: ccsd_force = .false.
 logical :: ccsd_t_force = .false.
 logical :: customized_core = .false.
 ! whether the user has specified the number of frozen core orbitals

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
 write(6,'(A)') '           Version: 1.2.6rc7 (2023-Jun-22)'
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
 call get_psi4_path(psi4_path)
 call getenv('GMS', gms_path)
 if(LEN_TRIM(gms_path) == 0) gms_path = 'NOT FOUND'

 write(6,'(A)') 'gau_path    = '//TRIM(gau_path)
 write(6,'(A)') 'gms_path    = '//TRIM(gms_path)
 write(6,'(A)') 'orca_path   = '//TRIM(orca_path)
 write(6,'(A)') 'molpro_path = '//TRIM(molpro_path)
 write(6,'(A)') 'psi4_path   = '//TRIM(psi4_path)
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
  case('mp2','ri-mp2','ccd','ccsd','ccsd(t)','ccsd(t)-f12','dlpno-ccsd',&
       'dlpno-ccsd(t)','dlpno-ccsd(t0)','dlpno-ccsd(t1)')
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
 case('mp2','ri-mp2')
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
  case('force')
   force = .true.
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
  case('fc')
   read(longbuf(j+1:i-1),*) core_wish
   customized_core = .true.
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

 cc_enabled = (ccd .or. ccsd .or. ccsd_t)
end subroutine parse_sr_keyword

! check the compatiblity of keywords in single reference calculations
subroutine check_sr_kywd_compatible()
 implicit none

 if(customized_core) then
  if(core_wish < 0) then
   write(6,'(/,A,I0)') 'ERROR in subroutine check_sr_kywd_compatible: invalid c&
                       &ore_wish=', core_wish
   write(6,'(A)') 'The number of frozen core orbitals must be non-negative.'
   stop
  end if
 end if

 if(ccd .and. TRIM(cc_prog)=='molpro') then
  write(6,'(/,A)') 'ERROR in subroutine check_sr_kywd_compatible: CC_prog=Molpr&
                   &o is incompatible with'
  write(6,'(A)') 'the CCD method, you may want to try CC_prog=ORCA.'
  stop
 end if

 if(RI) then
  if(TRIM(mp2_prog)=='gaussian' .or. TRIM(cc_prog)=='gaussian') then
   write(6,'(/,A)') 'ERROR in subroutine check_sr_kywd_compatible: RI is not su&
                    &pported in Gaussian.'
   write(6,'(A)') "If you still want to use Gaussian, you need to write 'noRI' &
                  &in mokit{}."
   stop
  end if
 end if

 if(RI .eqv. noRI) then
  write(6,'(A)') 'ERROR in subroutine check_sr_kywd_compatible: RI=noRI. Confus&
                 &ed keywords.'
  stop
 end if

 if(force) then
  if(mp2) mp2_force = .true.
  if(ccd) ccd_force = .true.
  if(ccsd) ccsd_force = .true.
  if(ccsd_t) ccsd_t_force = .true.
  if(RI .and. TRIM(mp2_prog)=='gamess') then
   write(6,'(/,A)') 'ERROR in subroutine check_sr_kywd_compatible: RI-MP2 analy&
                    &tical gradients'
   write(6,'(A)') 'are not supported in GAMESS.'
   stop
  end if
  if(cc_enabled .and. TRIM(cc_prog)=='orca') then
   write(6,'(/,A)') 'ERROR in subroutine check_sr_kywd_compatible: CC analytica&
                    &l gradients are'
   write(6,'(A)') 'not supported in ORCA. You can choose CC_prog=Molpro/PySCF/P&
                  &SI4, etc.'
   stop
  end if
  if(ccsd_t .and. TRIM(cc_prog)=='gaussian') then
   write(6,'(/,A)') 'ERROR in subroutine check_sr_kywd_compatible: CCSD(T) anal&
                    &ytical gradients'
   write(6,'(A)') 'are not supported in Gaussian. You can choose CC_prog=Molpro&
                  &/PSI4/PySCF.'
   stop
  end if
  if(TRIM(cc_prog) == 'gamess') then
   write(6,'(/,A)') 'ERROR in subroutine check_sr_kywd_compatible: CC analytica&
                    &l gradient is not'
   write(6,'(A)') 'supported in GAMESS.'
   stop
  end if
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
  write(6,'(A)') 'AutoSR 1.2.6rc7 :: MOKIT, release date: 2023-Jun-22'
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
  check_sr_kywd_compatible, hf_fch, chem_core, ecp_core
 implicit none
 character(len=24) :: data_string
 character(len=240), intent(in) :: fname

 gjfname = fname
 ! read paths of various programs from environment variables
 call read_sr_program_path()
 call parse_sr_keyword()
 call check_sr_kywd_compatible()

 call do_hf(.false.) ! RHF, ROHF, UHF
 call calc_ncore(hf_fch, chem_core, ecp_core) ! get number of core orbitals
 call do_mp2()   ! (DLPNO-, RI-)MP2
 call do_cc()    ! (DLPNO-)CCD, CCSD, CCSD(T), etc
 call do_cis()   ! CIS(D), ROCISD, XCIS, etc
 call do_adc()   ! ADC(2), ADC(3)
 call do_eomcc() ! EOM-CCSD, EOM-SF-CCSD

 call delete_file('fort.7')
 call fdate(data_string)
 write(6,'(/,A)') 'Normal termination of AutoSR at '//TRIM(data_string)
end subroutine autosr

! print frozen core information
subroutine prt_fc_info(chem_core, ecp_core, customized_core, core_wish)
 implicit none
 integer, intent(inout) :: chem_core
 integer, intent(in) :: ecp_core, core_wish
 logical, intent(in) :: customized_core
 logical :: fc

 if(customized_core) then
  if(chem_core /= core_wish) then
   write(6,'(/,A)') REPEAT('-',79)
   write(6,'(A)') 'Warning: you are changing the default number of frozen core &
                  &orbtials.'
   write(6,'(2(A,I0))') 'Recommended FC=',chem_core,', your request=',core_wish
   write(6,'(A)') 'Setting FC as your request. You should know what you are cal&
                  &culating.'
   write(6,'(A)') REPEAT('-',79)
   chem_core = core_wish
  end if
 end if

 fc = .true.
 if(chem_core == 0) fc = .false.
 write(6,'(A,L1,2(A,I0))') 'Frozen_core = ',fc,', chem_core=', chem_core, &
                           ', ecp_core=', ecp_core
end subroutine prt_fc_info

subroutine do_mp2()
 use sr_keyword, only: mem, nproc, hf_fch, bgchg, chgname, ref_e, corr_e, mp2_e,&
  mp2, mp2_prog, chem_core, ecp_core, customized_core, core_wish, gau_path, &
  psi4_path, gms_path, gms_scr_path, check_gms_path, orca_path, force
 use util_wrapper, only: unfchk, fch2mkl_wrap, mkl2gbw, bas_fch2py_wrap, &
  fch2com_wrap, fch2psi_wrap, fch2inp_wrap
 use mol, only: natom, grad
 implicit none
 integer :: i, RENAME, SYSTEM
 real(kind=8) :: rtmp
 character(len=24) :: data_string = ' '
 character(len=240) :: old_inp, inpname, chkname, mklname, outname, datname

 if(.not. mp2) return
 write(6,'(//,A)') 'Enter subroutine do_mp2...'
 write(6,'(A)') 'MP2 using program '//TRIM(mp2_prog)
 call prt_fc_info(chem_core, ecp_core, customized_core, core_wish)

 if(chem_core>0 .and. TRIM(mp2_prog)=='pyscf') then
  write(6,'(/,A)') 'ERROR in subroutine do_mp2: MP2 analytical gradients in PyS&
                   &CF cannot be run'
  write(6,'(A)') 'using frozen core. If you still want to use PySCF, you need t&
                 &o write FC=0 in mokit{}.'
  stop
 end if

 i = index(hf_fch, '.fch', back=.true.)
 outname = hf_fch(1:i-1)//'_MP2.out'

 select case(TRIM(mp2_prog))
 case('gaussian')
  inpname = hf_fch(1:i-1)//'_MP2.gjf'
  chkname = hf_fch(1:i-1)//'_MP2.chk'
  outname = hf_fch(1:i-1)//'_MP2.log'
  call unfchk(hf_fch, chkname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call prt_posthf_gau_inp(inpname)
  call submit_gau_job(gau_path, inpname)
  call read_mp2_e_from_gau_out(outname, ref_e, mp2_e)
 case('orca')
  inpname = hf_fch(1:i-1)//'_MP2.inp'
  old_inp = hf_fch(1:i-1)//'_o.inp'
  chkname = hf_fch(1:i-1)//'_MP2.engrad'
  mklname = hf_fch(1:i-1)//'_MP2.mkl'
  call fch2mkl_wrap(hf_fch, mklname)
  i = RENAME(TRIM(old_inp), TRIM(inpname))
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  ! if bgchg = .True., .inp and .mkl file will be updated
  call mkl2gbw(mklname)
  call delete_file(mklname)
  call prt_posthf_orca_inp(inpname)
  call submit_orca_job(orca_path, inpname)
  call read_posthf_e_from_orca_out(outname, .false., rtmp, ref_e, mp2_e)
 case('molpro')
  inpname = hf_fch(1:i-1)//'_MP2.com'
  call fch2com_wrap(hf_fch, inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call prt_posthf_molpro_inp(inpname)
  call submit_molpro_job(inpname, mem, nproc)
  call read_mp2_e_from_molpro_out(outname, ref_e, mp2_e)
 case('pyscf')
  inpname = hf_fch(1:i-1)//'_MP2.py'
  call bas_fch2py_wrap(hf_fch, .false., inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call prt_posthf_pyscf_inp(inpname)
  call submit_pyscf_job(inpname)
  call read_mp2_e_from_pyscf_out(outname, ref_e, mp2_e)
 case('psi4')
  inpname = hf_fch(1:i-1)//'_MP2.inp'
  call fch2psi_wrap(hf_fch, inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call prt_posthf_psi4_inp(inpname)
  call submit_psi4_job(psi4_path, inpname, nproc)
  !call read_mp2_e_from_psi4_out(outname, ref_e, mp2_e)
  stop
 case('gamess')
  call check_gms_path()
  datname = hf_fch(1:i-1)//'_MP2.dat'
  inpname = hf_fch(1:i-1)//'_MP2.inp'
  outname = hf_fch(1:i-1)//'_MP2.gms'
  call fch2inp_wrap(hf_fch, .false., 0, 0)
  old_inp = hf_fch(1:i-1)//'.inp'
  i = RENAME(TRIM(old_inp), TRIM(inpname))
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call prt_posthf_gms_inp(inpname)
  call submit_gms_job(gms_path, gms_scr_path, inpname, nproc)
  call read_mp2_e_from_gms_out(outname, ref_e, mp2_e)
 case default
  write(6,'(/,A)') 'ERROR in subroutine do_mp2: invalid MP2_prog='//TRIM(mp2_prog)
  stop
 end select

 corr_e = mp2_e - ref_e
 write(6,'(/,A,F18.8,A)')'E(ref)  = ', ref_e, ' a.u.'
 write(6,'(A,F18.8,A)')  'E(corr) = ', corr_e,' a.u.'
 write(6,'(A,F18.8,A)')  'E(MP2)  = ', mp2_e, ' a.u.'

 if(force) then
  allocate(grad(3*natom))

  select case(mp2_prog)
  case('gaussian')
   call read_grad_from_gau_log(outname, natom, grad)
  case('orca')
   call read_grad_from_engrad(chkname, natom, grad)
  case('molpro')
   call read_grad_from_molpro_out(outname, natom, grad)
  case('pyscf')
   call read_grad_from_pyscf_out(outname, natom, grad)
  case('gamess')
   call read_grad_from_gms_dat(datname, natom, grad)
  case default
   write(6,'(A)') 'ERROR in subroutine do_mp2: program cannot be identified.'
   write(6,'(A)') 'MP2_prog='//TRIM(mp2_prog)
   stop
  end select

  write(6,'(/,A)') 'Cartesian gradients (HARTREE/BOHR):'
  write(6,'(5(1X,ES15.8))') (grad(i),i=1,3*natom)
 end if

 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine do_mp2 at '//TRIM(data_string)
end subroutine do_mp2

subroutine do_cc()
 use sr_keyword, only: mem, nproc, ccd, ccsd, ccsd_t, cc_enabled, bgchg, hf_fch,&
  chgname, cc_prog, method, ccsd_e, ccsd_t_e, ref_e, corr_e, t1diag, chem_core,&
  ecp_core, customized_core, core_wish, RI, gau_path, orca_path, psi4_path, &
  gms_path, gms_scr_path, check_gms_path, force
 use util_wrapper, only: unfchk, fch2mkl_wrap, mkl2gbw, fch2com_wrap, &
  bas_fch2py_wrap, fch2psi_wrap, fch2inp_wrap
 use mol, only: natom, grad
 implicit none
 integer :: i, RENAME, SYSTEM
 real(kind=8) :: e = 0d0
 character(len=15) :: method0 = ' '
 character(len=24) :: data_string = ' '
 character(len=240) :: chkname, old_inp, inpname, mklname, outname

 if(.not. cc_enabled) return
 write(6,'(//,A)') 'Enter subroutine do_cc...'
 write(6,'(A)') 'CC using program '//TRIM(cc_prog)
 call prt_fc_info(chem_core, ecp_core, customized_core, core_wish)

 if(chem_core > 0) then
  select case(TRIM(cc_prog))
  case('pyscf')
   write(6,'(/,A)') 'ERROR in subroutine do_cc: CC analytical gradients in PySCF&
                    & cannot be run'
   write(6,'(A)') 'using frozen core. If you still want to use PySCF, you need t&
                  &o write FC=0 in mokit{}.'
   stop
  case('psi4')
   write(6,'(/,A)') 'ERROR in subroutine do_cc: CC analytical gradients in PSI4 &
                    &cannot be run'
   write(6,'(A)') 'using frozen core. If you still want to use PSI4, you need to&
                  & write FC=0 in mokit{}.'
   stop
  end select
 end if

 method0 = method
 call upper(method0)
 i = index(hf_fch, '.fch', back=.true.)
 outname = hf_fch(1:i-1)//'_CC.out'

 select case(TRIM(cc_prog))
 case('gaussian')
  inpname = hf_fch(1:i-1)//'_CC.gjf'
  chkname = hf_fch(1:i-1)//'_CC.chk'
  outname = hf_fch(1:i-1)//'_CC.log'
  call unfchk(hf_fch, chkname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call prt_posthf_gau_inp(inpname)
  call submit_gau_job(gau_path, inpname)
  call read_cc_e_from_gau_out(outname, t1diag, ref_e, e)
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
  call prt_posthf_orca_inp(inpname)
  call submit_orca_job(orca_path, inpname)
  call read_posthf_e_from_orca_out(outname, (.not.ccd), t1diag, ref_e, e)
 case('molpro')
  inpname = hf_fch(1:i-1)//'_CC.com'
  call fch2com_wrap(hf_fch, inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call prt_posthf_molpro_inp(inpname)
  call submit_molpro_job(inpname, mem, nproc)
  call read_cc_e_from_molpro_out(outname, (.not.ccd), t1diag, ref_e, e)
 case('pyscf')
  inpname = hf_fch(1:i-1)//'_CC.py'
  call bas_fch2py_wrap(hf_fch, .false., inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call prt_posthf_pyscf_inp(inpname)
  call submit_pyscf_job(inpname)
  call read_cc_e_from_pyscf_out(outname, t1diag, ref_e, e)
 case('psi4')
  inpname = hf_fch(1:i-1)//'_CC.inp'
  call fch2psi_wrap(hf_fch, inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call prt_posthf_psi4_inp(inpname)
  call submit_psi4_job(psi4_path, inpname, nproc)
  call read_cc_e_from_psi4_out(outname, RI, ccd, ccsd, ccsd_t, t1diag, ref_e, e)
 case('gamess')
  call check_gms_path()
  inpname = hf_fch(1:i-1)//'_CC.inp'
  outname = hf_fch(1:i-1)//'_CC.gms'
  call fch2inp_wrap(hf_fch, .false., 0, 0)
  old_inp = hf_fch(1:i-1)//'.inp'
  i = RENAME(TRIM(old_inp), TRIM(inpname))
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call prt_posthf_gms_inp(inpname)
  i = nproc
  if(ccd) i = 1 ! CCD in GAMESS is serial
  call submit_gms_job(gms_path, gms_scr_path, inpname, i)
  call read_cc_e_from_gms_out(outname, ccd, ccsd, ccsd_t, t1diag, ref_e, e)
 case default
  write(6,'(/,A)') 'ERROR in subroutine do_cc: invalid CC_prog='//TRIM(cc_prog)
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
 if(.not. ccd) write(6,'(A,F10.5,A)')  'T1_diag = ', t1diag

 if(force) then
  allocate(grad(3*natom))

  select case(cc_prog)
  case('pyscf')
   call read_grad_from_pyscf_out(outname, natom, grad)
  case('gaussian')
   call read_grad_from_gau_log(outname, natom, grad)
  case('openmolcas')
   call read_grad_from_molcas_out(outname, natom, grad)
  case('molpro')
   call read_grad_from_molpro_out(outname, natom, grad)
  case('psi4')
   call read_grad_from_psi4_out(outname, natom, grad)
  case default
   write(6,'(A)') 'ERROR in subroutine do_cc: program cannot be identified.'
   write(6,'(A)') 'CC_prog='//TRIM(cc_prog)
   stop
  end select

  write(6,'(/,A)') 'Cartesian gradients (HARTREE/BOHR):'
  write(6,'(5(1X,ES15.8))') (grad(i),i=1,3*natom)
 end if

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

! print a Gaussian MP2/CCSD/CCSD/CCSD(T) gjf file
subroutine prt_posthf_gau_inp(gjfname)
 use sr_keyword, only: mem, nproc, DKH2, mp2, ccd, ccsd, ccsd_t, chem_core, &
  force
 implicit none
 integer :: i, fid
 character(len=240), intent(in) :: gjfname

 i = index(gjfname, '.gjf', back=.true.)
 open(newunit=fid,file=TRIM(gjfname),status='replace')

 write(fid,'(A)') '%chk='//gjfname(1:i-1)//'.chk'
 write(fid,'(A,I0,A)') '%mem=', mem, 'GB'
 write(fid,'(A,I0)') '%nprocshared=',nproc
 write(fid,'(A)',advance='no') '#p '

 if(mp2) then
  write(fid,'(A)',advance='no') 'MP2'
 else if(ccd) then
  write(fid,'(A)',advance='no') 'CCD'
 else if(ccsd) then
  write(fid,'(A)',advance='no') 'CCSD(T1diag)'
 else if(ccsd_t) then
  write(fid,'(A)',advance='no') 'CCSD(T,T1diag)'
 end if

 write(fid,'(A)',advance='no') ' chkbasis nosymm guess=read geom=allcheck'

 if(DKH2) then
  write(fid,'(A)',advance='no') ' int(nobasistransform,DKH2) iop(3/93=1)'
 else
  write(fid,'(A)',advance='no') ' int=nobasistransform'
 end if

 write(fid,'(A,I0,A)',advance='no') ' window=(',chem_core+1,',0)'
 if(force) write(fid,'(A)',advance='no') ' force'

 write(fid,'(/)',advance='no')
 close(fid)
end subroutine prt_posthf_gau_inp

! add CC keywords into a ORCA input file
subroutine prt_posthf_orca_inp(inpname)
 use sr_keyword, only: mem, nproc, RI, DLPNO, F12, RIJK_bas, RIC_bas, F12_cabs,&
  mo_rhf, mp2, ccd, ccsd, ccsd_t, iterative_t, chem_core, force
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
 write(fid1,'(A)',advance='no') 'HF VeryTightSCF'
 if(RI) then
  write(fid1,'(A)',advance='no') ' RIJK '//TRIM(RIJK_bas)//' '//TRIM(RIC_bas)
 else
  write(fid1,'(A)',advance='no') ' noRI'
 end if
 if(F12) write(fid1,'(A)',advance='no') ' '//TRIM(F12_cabs)
 if(chem_core == 0) write(fid1,'(A)',advance='no') ' NoFrozenCore'
 if(force) write(fid1,'(A)',advance='no') ' EnGrad'

 if(DLPNO) then
  write(fid1,'(A)',advance='no') ' TightPNO DLPNO-'
  if(iterative_t) then
   write(fid1,'(A)',advance='no') 'CCSD(T1)'
  else if(ccsd_t) then
   write(fid1,'(A)',advance='no') 'CCSD(T)'
  else if(mp2) then
   write(fid1,'(A)',advance='no') 'MP2'
  end if
 else
  if(ccsd_t) then
   write(fid1,'(A)',advance='no') ' CCSD(T)'
  else if(ccsd) then
   write(fid1,'(A)',advance='no') ' CCSD'
  else if(ccd) then
   write(fid1,'(A)',advance='no') ' CCD'
  else if(mp2) then
   write(fid1,'(A)',advance='no') ' MP2'
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
end subroutine prt_posthf_orca_inp

subroutine prt_posthf_molpro_inp(inpname)
 use sr_keyword, only: mp2, ccsd_t, RI, F12, mo_rhf, uhf_based, basis, RIJK_bas,&
  RIC_bas, chem_core, force
 use mol, only: mult
 implicit none
 integer :: i, fid
 character(len=21) :: MP2FIT_bas, RIJK_bas1
 character(len=240), intent(in) :: inpname

 if(RI) then
  call auxbas_convert(RIJK_bas, RIJK_bas1, 2)
  i = index(RIJK_bas, '/')
  RIC_bas = RIJK_bas(1:i)//'optri'
  MP2FIT_bas = RIJK_bas(1:i)//'mp2fit'
 end if

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
  write(fid,'(A)',advance='no') '{DF-'
  if(mult /= 1) then
   if(mo_rhf) then
    if(uhf_based) then  ! ROHF-UCCSD
     write(fid,'(A)',advance='no') 'U'
    else                ! ROHF-ROCCSD
     write(fid,'(A)',advance='no') 'R'
    end if
   else
    write(6,'(A)') 'ERROR in subroutine prt_posthf_molpro_inp: UHF-based CC is &
                   &not supported in Molpro.'
    close(fid)
    stop
   end if
  end if
  if(mp2) then
   write(fid,'(A)',advance='no') 'MP2'
  else
   write(fid,'(A)',advance='no') 'CCSD'
  end if
  if(ccsd_t) write(fid,'(A)',advance='no') '(T)'
  if(F12) write(fid,'(A)',advance='no') '-F12,ri_basis=ri'
  write(fid,'(A)',advance='no') ',df_basis=df,df_basis_exch=jk'

 else ! RI is turned off
  if(mp2) then
   write(fid,'(A)',advance='no') '{MP2'
  else
   write(fid,'(A)',advance='no') '{CCSD'
  end if
  if(ccsd_t) write(fid,'(A)',advance='no') '(T)'
  if(F12) write(fid,'(A)',advance='no') '-F12'
 end if

 write(fid,'(A,I0,A)') ';core,',chem_core,'}'
 if(force) write(fid,'(A)') 'Force'
 close(fid)
end subroutine prt_posthf_molpro_inp

subroutine prt_posthf_pyscf_inp(pyname)
 use sr_keyword, only: mem, nproc, mp2, ccd, ccsd_t, chem_core, force, RI
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=240) :: buf, pyname1
 character(len=240), intent(in) :: pyname

 i = index(pyname, '.py', back=.true.)
 pyname1 = pyname(1:i-1)//'.t'
 open(newunit=fid,file=TRIM(pyname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(pyname1),status='replace')

 write(fid1,'(A)',advance='no') 'from pyscf import lib, '
 if(mp2) then
  write(fid1,'(A)') 'mp'
 else
  write(fid1,'(A)') 'cc'
 end if
 if(ccsd_t) write(fid1,'(A)') 'from pyscf.cc import ccsd_t'

 do while(.true.)
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while
 write(fid1,'(/,A,I0,A)') 'lib.num_threads(',nproc,')'

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:3) == '#dm') exit
  if(buf(1:9) == 'mf.kernel') then
   write(fid1,'(A,I0,A)') 'mf.max_memory = ',mem*1000,' #MB'
  end if
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 write(fid1,'(A)') TRIM(buf(2:))
 read(fid,'(A)') buf
 write(fid1,'(A)') TRIM(buf(2:))
 read(fid,'(A)') buf
 write(fid1,'(A)') TRIM(buf(2:))
 close(fid,status='delete')

 if(mp2) then
  if(RI) then
   write(fid1,'(/,A)') 'mc = mp.dfmp2.DFMP2(mf)'
  else
   write(fid1,'(/,A)') 'mc = mp.MP2(mf)'
  end if
 else
  write(fid1,'(/,A)') 'mc = cc.CCSD(mf)'
 end if
 write(fid1,'(A)') 'mc.verbose = 4'
 if(ccd) then
  write(fid1,'(A)') 'old_update_amps = mc.update_amps'
  write(fid1,'(A)') 'def update_amps(t1, t2, eris):'
  write(fid1,'(A)') '  t1, t2 = old_update_amps(t1, t2, eris)'
  write(fid1,'(A)') '  return t1*0, t2'
  write(fid1,'(A)') 'mc.update_amps = update_amps'
 end if
 write(fid1,'(A,I0)') 'mc.frozen = ', chem_core
 write(fid1,'(A)') 'mc.kernel()'
 if((.not.mp2) .and. (.not.ccd)) then
  write(fid1,'(A)') 'T1diag = mc.get_t1_diagnostic()'
  write(fid1,'(A)') "print('T1_diag =', T1diag)"
 end if
 if(ccsd_t) write(fid1,'(A)') 'mc.ccsd_t()'

 close(fid1)
 i = RENAME(TRIM(pyname1), TRIM(pyname))
 if(force) call add_force_key2py_script(mem, pyname)
end subroutine prt_posthf_pyscf_inp

subroutine prt_posthf_psi4_inp(inpname)
 use sr_keyword, only: mem, mp2, ccd, ccsd, ccsd_t, chem_core, RI, force
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 inpname1 = TRIM(inpname)//'.t'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 read(fid,'(A)') buf
 write(fid1,'(A)') TRIM(buf)
 read(fid,'(A)') buf
 write(fid1,'(A,I0,A)') 'memory ',mem,' GB'

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:9) == 'scf_type') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(RI) then
  write(fid1,'(A)') ' scf_type df'
 else
  write(fid1,'(A)') ' scf_type pk'
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:11) == 'scfenergy =') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid, status='delete')
 write(fid1,'(A,I0)') 'set freeze_core ', chem_core
 if(RI) then
  write(fid1,'(A)') 'set cc_type df'
 else
  write(fid1,'(A)') 'set cc_type conv'
  if(ccd) then
   write(6,'(/,A)') 'ERROR in subroutine prt_posthf_psi4_inp: conventional CCD &
                    &is not supported'
   write(6,'(A)') 'in PSI4. Please turn on the RI approximation.'
   close(fid1)
   stop
  end if
 end if

 if(force) then
  write(fid1,'(A)',advance='no') "gradient('"
 else
  write(fid1,'(A)',advance='no') "energy('"
 end if
 if(mp2) then
  write(fid1,'(A)') "mp2')"
 else if(ccd) then
  write(fid1,'(A)') "ccd')"
 else if(ccsd) then
  write(fid1,'(A)') "ccsd')"
 else if(ccsd_t) then
  write(fid1,'(A)') "ccsd(t)')"
 end if
 if((.not.mp2) .and. (.not.ccd)) write(fid1,'(A)') "print_variables()"

 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine prt_posthf_psi4_inp

subroutine prt_posthf_gms_inp(inpname)
 use sr_keyword, only: mem, nproc, mp2, ccd, ccsd, ccsd_t, cc_enabled, force, &
  chem_core, RI, uhf_based
 implicit none
 integer :: i, j, fid, fid1, RENAME
 character(len=13) :: method
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 if(mp2) then
  method = 'MPLEVL=2'
 else if(ccd) then
  method = 'CCTYP=CCD'
 else if(ccsd) then
  method = 'CCTYP=CCSD'
 else if(ccsd_t) then
  method = 'CCTYP=CCSD(T)'
 end if

 inpname1 = TRIM(inpname)//'.t'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 read(fid,'(A)') buf
 if(force) then
  i = index(buf, 'ENERGY')
  buf = buf(1:i-1)//'GRADIENT'//TRIM(buf(i+6:))
 end if
 write(fid1,'(A)') TRIM(buf)

 read(fid,'(A)') buf
 i = index(buf, '$END')
 write(fid1,'(A)') buf(1:i-1)//TRIM(method)//' $END'

 read(fid,'(A)') buf ! assuming this is $SYSTEM
 i = FLOOR(DBLE(mem)*475d0/(4d0*DBLE(nproc)))
 j = FLOOR(DBLE(mem)*25d0/4d0)
 write(fid1,'(2(A,I0),A)') " $SYSTEM MWORDS=",i,' MEMDDI=',j," $END"

 if(mp2) then
  write(fid1,'(A,I0)',advance='no') " $MP2 NACORE=", chem_core
  if(uhf_based) write(fid1,'(A,I0)',advance='no') " NBCORE=", chem_core
 else if(cc_enabled) then
  write(fid1,'(A,I0)',advance='no') " $CCINP NCORE=", chem_core
 end if
 if(RI) then
  if(mp2) then
   write(fid1,'(A)',advance='no') ' CODE=RIMP2'
  else if(cc_enabled) then
   write(fid1,'(A)',advance='no') ' CCERI=RI'
  end if
 end if
 write(fid1,'(A)') " $END"

 if(RI) then
  if(mp2) then
   write(fid1,'(A)') " $AUXBAS CABNAM=TZVPP $END"
  else if(cc_enabled) then
   write(fid1,'(A)') " $RICC CABNAM=TZVPP $END"
  end if
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine prt_posthf_gms_inp

subroutine read_mp2_e_from_gau_out(outname, ref_e, tot_e)
 implicit none
 integer :: i, fid
 real(kind=8) :: ss
 real(kind=8), intent(out) :: ref_e, tot_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 ref_e = 0d0; tot_e = 0d0
 call read_hf_e_and_ss_from_gau_out(outname, ref_e, ss)

 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(28:34) == 'EUMP2 =') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_mp2_e_from_gau_out: 'EUMP2 =' not &
                   &found in"
  write(6,'(A)') 'file '//TRIM(outname)
  close(fid)
  stop
 end if

 read(buf(35:),*) tot_e ! MP2 total energy
 close(fid)
end subroutine read_mp2_e_from_gau_out

subroutine read_mp2_e_from_molpro_out(outname, ref_e, tot_e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: ref_e, tot_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 ref_e = 0d0; tot_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(index(buf,'MP2') > 0) exit 
 end do ! for while

 read(fid,*) tot_e, ref_e
 close(fid)
end subroutine read_mp2_e_from_molpro_out

subroutine read_mp2_e_from_pyscf_out(outname, ref_e, tot_e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: ref_e, tot_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 ref_e = 0d0; tot_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(1:8)=='E(MP2) =' .or. buf(1:10)=='E(DFMP2) =') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_mp2_e_from_pyscf_out: 'E(MP2) =' n&
                   &ot found in"
  write(6,'(A)') 'file '//TRIM(outname)
  stop
 end if

 i = index(buf, '=')
 read(buf(i+1:),*) tot_e
 i = index(buf, '=', back=.true.)
 read(buf(i+1:),*) ref_e ! here is MP2 correlation energy
 ref_e = tot_e - ref_e
end subroutine read_mp2_e_from_pyscf_out

subroutine read_scf_e_from_gms_out(outname, scf_e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: scf_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 scf_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:7) == ' FINAL ') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_scf_e_from_gms_out: ' FINAL ' not &
                   &found in file "//TRIM(outname)
  close(fid)
  stop
 end if

 i = index(buf, 'IS')
 buf = ADJUSTL(buf(i+2:))
 read(buf,*) scf_e
end subroutine read_scf_e_from_gms_out

subroutine read_mp2_e_from_gms_out(outname, ref_e, tot_e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: ref_e, tot_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 ref_e = 0d0; tot_e = 0d0
 call read_scf_e_from_gms_out(outname, ref_e)

 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(index(buf,'E(MP2)=') > 0) exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_mp2_e_from_gms_out: no 'E(MP2)=' f&
                   &ound in file "//TRIM(outname)
  stop
 end if

 i = index(buf,'E(MP2)=')
 read(buf(i+7:),*) tot_e
end subroutine read_mp2_e_from_gms_out

subroutine read_cc_e_from_gau_out(outname, t1diag, ref_e, tot_e)
 implicit none
 integer :: i, fid
 real(kind=8) :: ss
 real(kind=8), intent(out) :: t1diag, ref_e, tot_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 t1diag = 0d0; ref_e = 0d0; tot_e = 0d0
 call read_hf_e_and_ss_from_gau_out(outname, ref_e, ss)

 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:17) == 'Wavefunction amp') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_orca_out: 'Wavefunction &
                   &amp' not found"
  write(6,'(A)') 'in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = index(buf, '=')
 read(buf(i+1:),*) tot_e ! CCD/CCSD total energy

 read(fid,'(A)') buf
 if(buf(2:8) == 'T1 Diag') then
  i = index(buf, '=')
  read(buf(i+1:),*) t1diag
 end if

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:14) == 'Discarding MO') exit
  if(buf(2:8) == 'CCSD(T)') then
   i = index(buf, '=')
   read(buf(i+1:),*) tot_e ! CCSD(T) total energy
   exit
  end if
 end do ! for while

 close(fid)
end subroutine read_cc_e_from_gau_out

subroutine find_hf_type_in_orca_out(outname, hf_type)
 implicit none
 integer :: i, fid
 integer, intent(out) :: hf_type
 character(len=4) :: str
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:18) == 'Hartree-Fock type') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine find_hf_type_in_orca_out: no 'Hartree-F&
                   &ock type' found "
  write(6,'(A)') 'in file '//TRIM(outname)
  stop
 end if

 i = index(buf, '.', back=.true.)
 read(buf(i+1:),*) str

 select case(TRIM(str))
 case('RHF')
  hf_type = 1
 case('ROHF')
  hf_type = 2
 case('UHF')
  hf_type = 3
 case default
  write(6,'(/,A)') 'ERROR in subroutine find_hf_type_in_orca_out: invalid HF_ty&
                   &pe='//TRIM(str)
  stop
 end select
end subroutine find_hf_type_in_orca_out

subroutine read_posthf_e_from_orca_out(outname, read_t1diag, t1diag, ref_e, tot_e)
 implicit none
 integer :: i, hf_type, fid
 real(kind=8) :: ss
 real(kind=8), intent(out) :: t1diag, ref_e, tot_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: read_t1diag

 t1diag = 0d0; ref_e = 0d0; tot_e = 0d0
 call find_hf_type_in_orca_out(outname, hf_type)
 call read_hf_e_and_ss_from_orca_out(outname, hf_type, ref_e, ss)

 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

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
  write(6,'(/,A)') "ERROR in subroutine read_posthf_e_from_orca_out: 'FINAL SIN&
                   &G' not found in"
  write(6,'(A)') 'file '//TRIM(outname)
  stop
 end if

 read(buf(26:),*) tot_e
end subroutine read_posthf_e_from_orca_out

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

subroutine read_cc_e_from_pyscf_out(outname, t1diag, ref_e, tot_e)
 implicit none
 integer :: i, k, fid
 real(kind=8) :: ccsd_t_corr = 0d0
 real(kind=8), intent(out) :: t1diag, ref_e, tot_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 t1diag = 0d0; ref_e = 0d0; tot_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')
 write(fid,'(/)',advance='no')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf

  select case(buf(1:9))
  case('CCSD(T) c')
   k = index(buf, '=')
   read(buf(k+1:),*) ccsd_t_corr
  case('E(CCSD) =')
   read(buf(10:),*) tot_e
  case('T1_diag =')
   read(buf(10:),*) t1diag
  case('converged')
   exit
  end select
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine read_cc_e_from_pyscf_out: incomplete fil&
                 &e '//TRIM(outname)
  stop
 end if

 k = index(buf, '=')
 read(buf(k+1:),*) ref_e
 tot_e = tot_e + ccsd_t_corr
end subroutine read_cc_e_from_pyscf_out

subroutine read_cc_e_from_psi4_out(outname, RI, ccd, ccsd, ccsd_t, t1diag, ref_e,&
                                   tot_e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: t1diag, ref_e, tot_e
 character(len=14), parameter :: key0(5) = ['DF-CCD Total E','* CCSD total e',&
                           'DF-CCSD Total ','* CCSD(T) tota','DF-CCSD(T) Tot']
 character(len=14) :: key
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: RI, ccd, ccsd, ccsd_t

 t1diag = 0d0; ref_e = 0d0; tot_e = 0d0
 if(ccd) then
  key = key0(1)
 else if(ccsd) then
  if(RI) then
   key = key0(3)
  else
   key = key0(2)
  end if
 else if(ccsd_t) then
  if(RI) then
   key = key0(5)
  else
   key = key0(4)
  end if
 end if

 open(newunit=fid,file=TRIM(outname),status='old',position='append')
 if(.not. ccd) then
  do while(.true.)
   BACKSPACE(fid, iostat=i)
   if(i /= 0) exit
   BACKSPACE(fid, iostat=i)
   if(i /= 0) exit
   read(fid,'(A)') buf
   if(buf(1:6) == 'memory') then
    i = -1; exit
   end if
   if(index(buf,'T1 DIAG')>0 .or. index(buf,'T1 diag')>0) exit
  end do ! for while
  i = index(buf, '=>')
  if(i == 0) i = index(buf, ':')
  read(buf(i+2:),*) t1diag
 end if

 rewind(fid)
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(index(buf, key) > 0) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_cc_e_from_psi4_out: key '"//key//"'"
  write(6,'(A)') 'not found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = index(buf, '=')
 if(i == 0) i = index(buf, ':')
 read(buf(i+1:),*) tot_e

 do while(.true.)
  BACKSPACE(fid, iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid, iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(1:6) == 'memory') then
   i = -1; exit
  end if
  if(index(buf,'SCF E')>0 .or. index(buf,'SCF e')>0) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_cc_e_from_psi4_out: 'SCF E' or 'SCF&
                 & e' not found in"
  write(6,'(A)') 'file '//TRIM(outname)
  close(fid)
  stop
 end if

 close(fid)
 i = index(buf, '=')
 if(i == 0) i = index(buf, ':')
 read(buf(i+1:),*) ref_e
end subroutine read_cc_e_from_psi4_out

subroutine read_cc_e_from_gms_out(outname, ccd, ccsd, ccsd_t, t1diag, ref_e, &
                                  tot_e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: t1diag, ref_e, tot_e
 character(len=10), parameter :: key0(3) = ['CCD ENERGY','CCSD    EN','CCSD(T) EN']
 character(len=10) :: key
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: ccd, ccsd, ccsd_t

 t1diag = 0d0; ref_e = 0d0; tot_e = 0d0
 if(ccd) then
  key = key0(1)
 else if(ccsd) then
  key = key0(2)
 else if(ccsd_t) then
  key = key0(3)
 end if

 call read_scf_e_from_gms_out(outname, ref_e)
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 if(.not. ccd) then
  do while(.true.)
   read(fid,'(A)') buf
   if(buf(2:8) == 'T1 DIAG') exit
  end do ! for while
  i = index(buf, '=')
  read(buf(i+1:),*) t1diag
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(index(buf, key) > 0) exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_cc_e_from_gms_out: key '"//key//"'"
  write(6,'(A)') 'not found in file '//TRIM(outname)
  stop
 end if

 i = index(buf, ':')
 read(buf(i+1:),*) tot_e
end subroutine read_cc_e_from_gms_out

