! written by jxzou at 20221230: automatic single reference calculations

! keywords information of single reference calculations (default values are set)
module sr_keyword
 use mr_keyword, only: gjfname, mem, nproc, method, basis, bgchg, cart, force, &
  DKH2, X2C, dkh2_or_x2c, RI, F12, DLPNO, RIJK_bas, RIC_bas, F12_cabs, localm, &
  nstate, readrhf, readuhf, mo_rhf, hf_prog, hfonly, hf_fch, skiphf, chgname, &
  gau_path, molcas_path, orca_path, psi4_path, dalton_path, gms_path, gms_scr_path,&
  nmr, check_gms_path, molcas_omp, dalton_mpi
 use mol, only: chem_core, ecp_core, ptchg_e, nuc_pt_e, lin_dep
 implicit none
 integer :: core_wish = 0 ! the number of frozen core orbitals the user wishes
 real(kind=8) :: ref_e = 0d0    ! reference wfn energy
 real(kind=8) :: corr_e = 0d0   ! total correlation energy
 real(kind=8) :: mp2_e = 0d0    ! MP2 total energy
 real(kind=8) :: ccd_e = 0d0    ! CCD total energy
 real(kind=8) :: ccsd_e = 0d0   ! CCSD total energy
 real(kind=8) :: ccsd_t_e = 0d0 ! CCSD(T) total energy
 real(kind=8) :: t1diag = 0d0   ! T1 diagnostic
 real(kind=8), allocatable :: ex_elec_e(:)
 character(len=10) :: mp2_prog = 'orca'
 character(len=10) :: cc_prog = 'pyscf' ! very fast CCSD and (T)
 character(len=10) :: adc_prog = 'pyscf'
 character(len=10) :: eom_prog = 'orca'
 character(len=240) :: no_fch = ' ' ! filename which includes MP2/CC NOs
 logical :: noRI = .false.          ! turn on RI by default
 logical :: mp2 = .false.
 logical :: ccd = .false.
 logical :: ccsd = .false.
 logical :: ccsd_t = .false.
 logical :: cc_enabled = .false. ! (ccd .or. ccsd .or. ccsd_t)
 logical :: cis = .false.
 logical :: adc = .false.
 logical :: eom = .false.
 logical :: ip = .false.
 logical :: ea = .false.
 logical :: iterative_t = .false. ! default DLPNO-CCSD(T0)
 logical :: gen_no = .false.      ! whether to generate NOs
 logical :: relaxed_dm = .false.  ! relaxed/unrelaxed density matrix
 logical :: customized_core = .false.
 ! whether the user has specified the number of frozen core orbitals

contains

subroutine read_sr_program_path()
 use mr_keyword, only: mokit_root, gau_path, molpro_path, get_molcas_path
 implicit none
 integer :: i
 integer(kind=4) :: hostnm
 character(len=8) :: hostname
 character(len=24) :: data_string
 character(len=240), external :: get_mokit_root 

 write(6,'(A)') '------ Output of AutoSR of MOKIT(Molecular Orbital Kit) ------'
 write(6,'(A)') '       GitLab page: https://gitlab.com/jxzou/mokit'
 write(6,'(A)') '     Documentation: https://jeanwsr.gitlab.io/mokit-doc-mdbook'
 write(6,'(A)') '           Version: 1.2.6rc39 (2024-Sep-14)'
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
 call get_molcas_path()
 call check_molcas_is_omp(molcas_omp)
 call get_orca_path(orca_path)
 call get_molpro_path(molpro_path)
 call get_psi4_path(psi4_path)
 call get_dalton_path(dalton_path)
 if(TRIM(dalton_path) /= 'NOT FOUND') call check_dalton_is_mpi(dalton_mpi)
 call getenv('GMS', gms_path)
 if(LEN_TRIM(gms_path) == 0) gms_path = 'NOT FOUND'

 write(6,'(A)') 'gau_path    = '//TRIM(gau_path)
 write(6,'(A)') 'gms_path    = '//TRIM(gms_path)
 write(6,'(A)') 'orca_path   = '//TRIM(orca_path)
 write(6,'(A)') 'molpro_path = '//TRIM(molpro_path)
 write(6,'(A)') 'molcas_path = '//TRIM(molcas_path)
 write(6,'(A)') 'psi4_path   = '//TRIM(psi4_path)
 write(6,'(A)') 'dalton_path = '//TRIM(dalton_path)
end subroutine read_sr_program_path

! Parse keywords of single reference calculations
subroutine parse_sr_keyword()
 use mr_keyword, only: basname, check_readfch, lower
 implicit none
 integer :: i, j, k, ifail, fid
 character(len=24) :: method0 = ' '
 character(len=240) :: buf = ' '
 character(len=1000) :: longbuf = ' '

 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
 call read_mem_nproc_route(fid, mem, nproc, buf)

 i = INDEX(buf,'/')
 if(i == 0) then
  write(6,'(/,A)') "ERROR in subroutine parse_sr_keyword: no '/' symbol found i&
                   &n the keyword line."
  write(6,'(A)') "The method and basis set must be specified via '/' symbol, e.&
                 &g. CCSD(T)/cc-pVTZ."
  close(fid)
  stop
 end if
 if(mem < nproc) then
  write(6,'(/,A)') 'ERROR in subroutine parse_sr_keyword: please specify larger&
                   & memory. Post-'
  write(6,'(A)') 'HF calculations usually requires large memory.'
  stop
 end if

 j = INDEX(buf(1:i-1),' ', back=.true.)
 if(j == 0) then
  write(6,'(/,A)') 'ERROR in subroutine parse_sr_keyword: syntax error detected'
  write(6,'(A)') "in the current line '"//TRIM(buf)//"'"
  stop
 end if
 method0 = buf(j+1:i-1)

 if(i /= 0) then
  method = method0(1:i-1)

  select case(TRIM(method))
  case('mp2','ri-mp2','ccd','ccsd','ccsd(t)','ccsd(t)-f12','ccsd(t)-f12a', &
       'ccsd(t)-f12b','dlpno-ccsd','dlpno-ccsd(t)','dlpno-ccsd(t0)', &
       'dlpno-ccsd(t1)','eomccsd','eom-ccsd','eom-ip-ccsd','ip-eom-ccsd', &
       'eom-ea-ccsd','ea-eom-ccsd')
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
 case('ccsd(t)-f12','ccsd(t)-f12a','ccsd(t)-f12b')
  ccsd_t = .true.; F12 = .true.
 case('dlpno-mp2')
  mp2 = .true.; DLPNO = .true.
 case('dlpno-ccsd')
  ccsd = .true.; DLPNO = .true.
 case('dlpno-ccsd(t)','dlpno-ccsd(t0)')
  ccsd_t = .true.; DLPNO = .true.
 case('dlpno-ccsd(t1)')
  ccsd_t = .true.; DLPNO = .true.; iterative_t = .true.
 case('eom-ip-ccsd','ip-eom-ccsd')
  ccsd = .true.; eom = .true.; ip = .true.
 case('eom-ea-ccsd','ea-eom-ccsd')
  ccsd = .true.; eom = .true.; ea = .true.
 case('eomccsd','eom-ccsd')
  ccsd = .true.; eom = .true.
 case default
  write(6,'(/,A)') 'ERROR in subroutine parse_sr_keyword: unrecognized method='&
                  //TRIM(method)
  stop
 end select

 i = INDEX(buf,'/')
 j = i - 1 + INDEX(buf(i+1:),' ')
 if(j == 0) j = LEN_TRIM(buf)
 basis = buf(i+1:j)
 write(6,'(/,2(A,I4))',advance='no') 'memory =', mem, 'GB, nproc =', nproc
 write(6,'(A)') ', method/basis = '//TRIM(method)//'/'//TRIM(basis)

 if(basis(1:5) == 'def2-') then
  write(6,'(/,A)') "ERROR in subroutine parse_sr_keyword: 'def2-' prefix detect&
                   &ed in given basis set."
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
 ! This lower() is not from string_manipulate.f90, but from module mr_keyword

 if(buf(1:6) /= 'mokit{') then
  write(6,'(/,A)') "ERROR in subroutine parse_sr_keyword: 'mokit{' not detected&
                  & in file "//TRIM(gjfname)
  write(6,'(A)') "Syntax error. You must put 'mokit{' in leading position of th&
                 &e Title Card line."
  stop
 end if

 j = INDEX(buf,'}')
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

 if(INDEX(buf,'}') == 0) then ! keywords are written in at least two lines
  do while(.true.)
   read(fid,'(A)',iostat=ifail) buf
   if(ifail /= 0) exit

   call lower(buf)
   i = LEN_TRIM(buf)
   if(INDEX(buf,'}') > 0) i = i - 1
   longbuf(k:k+i-1) = buf(1:i)
   k = k + i

   if(INDEX(buf,'}') > 0) exit
  end do ! for while

  if(ifail /= 0) then
   write(6,'(/,A)') 'ERROR in subroutine parse_sr_keyword: end-of-file detected.'
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
  i = INDEX(longbuf,',')
  if(i == 0) i = LEN_TRIM(longbuf) + 1

  j = INDEX(longbuf(1:i-1),'=')
  if(j == 0) j = i   ! in case this keyword has no value assigned

  select case(longbuf(1:j-1))
  case('readrhf')
   mo_rhf = .true.
   readrhf = .true.
   read(longbuf(j+1:i-1),*) hf_fch
  case('readuhf')
   mo_rhf = .false.
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
  case('eom_prog')
   read(longbuf(j+1:i-1),*) eom_prog
  case('fc')
   read(longbuf(j+1:i-1),*) core_wish
   customized_core = .true.
  case('nstates')
   read(longbuf(j+1:i-1),*) nstate
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
  case('no')
   gen_no = .true.
  case('relaxed_dm')
   relaxed_dm = .true.
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

 if(readrhf .or. readuhf) then
  skiphf = .true.
  call check_sanity_of_provided_fch(DKH2, X2C, hf_fch)
 end if
end subroutine parse_sr_keyword

! check the compatiblity of keywords in single reference calculations
subroutine check_sr_kywd_compatible()
 implicit none
 character(len=46), parameter :: error_warn = 'ERROR in subroutine check_sr_kyw&
                                              &d_compatible: '
 cc_enabled = (ccd .or. ccsd .or. ccsd_t)

 if(cc_enabled .and. DLPNO .and. TRIM(cc_prog)/='orca') then
  write(6,'(/,A)') REPEAT('-',79)
  write(6,'(A)') 'Warning from subroutine check_sr_kywd_compatible: a DLPNO-CC &
                 &job is requested.'
  write(6,'(A)') 'CC_prog is automatically switched to ORCA.'
  write(6,'(A)') REPEAT('-',79)
  cc_prog = 'orca'
 end if

 if(customized_core) then
  if(core_wish < 0) then
   write(6,'(/,A,I0)') error_warn//'invalid core_wish=', core_wish
   write(6,'(A)') 'The number of frozen core orbitals must be non-negative.'
   stop
  end if
 end if

 if(ccd .and. TRIM(cc_prog)=='openmolcas') then
  write(6,'(/,A)') error_warn//'CCD not supported in OpenMolcas.'
  write(6,'(A)') 'You can specify another CC_prog.'
  stop
 end if

 if(gen_no) then
  if(force .and. (.not.relaxed_dm)) then
   write(6,'(/,A)') error_warn//'force v.s. unrelaxed density is'
   write(6,'(A)') 'incompatible, since analytic force means a calculation of re&
                  &laxed density. You'
   write(6,'(A)') 'are supposed to choose only one of force/unrelaxed density.'
   stop
  end if
  if(cc_enabled .and. relaxed_dm .and. TRIM(cc_prog)=='orca') then
   write(6,'(/,A)') error_warn//'CC relaxed density are'
   write(6,'(A)') 'not supported when CC_prog=ORCA.'
   stop
  end if
  if(ccsd_t) then
   select case(TRIM(cc_prog))
   case('orca','gaussian','qchem')
    write(6,'(/,A)') error_warn//'CCSD(T) density is not'
    write(6,'(A)') 'supported by Gaussian/ORCA/Q-Chem. You may try CC_prog=Molp&
                   &ro/PSI4.'
    stop
   end select
  end if
  if(cc_enabled .and. (.not. relaxed_dm)) then
   select case(TRIM(cc_prog))
   case('gaussian','openmolcas')
    write(6,'(/,A)') error_warn//'CC unrelaxed density is'
    write(6,'(A)') 'not supported when CC_prog=Gaussian/OpenMolcas.'
    stop
   end select
  end if
  if(cc_enabled .and. relaxed_dm) then
   select case(TRIM(cc_prog))
   case('pyscf','openmolcas')
    write(6,'(/,A)') error_warn//'CC relaxed density is'
    write(6,'(A)') 'not supported when CC_prog=PySCF/OpenMolcas.'
    stop
   end select
  end if
  if(mp2) then
   if(noRI .and. TRIM(mp2_prog)=='psi4') then
    write(6,'(/,A)') error_warn//'MP2 density is not supported'
    write(6,'(A)') 'for MP2_prog=PSI4 when RI is turned off. Please turn RI on.'
    stop
   end if
   if((.not.relaxed_dm) .and. TRIM(mp2_prog)=='gamess') then
    write(6,'(/,A)') error_warn//'MP2 unrelaxed density is not'
    write(6,'(A)') 'supported MP2_prog=GAMESS. You can use another MP2_prog. If&
                   & you insist on'
    write(6,'(A)') 'using GAMESS, you can calculate relaxed density instead.'
    stop
   end if
   if(noRI .and. TRIM(mp2_prog)=='pyscf' .and. relaxed_dm) then
    write(6,'(/,A)') error_warn//'MP2 relaxed density is not'
    write(6,'(A)') 'supported when MP2_prog=PySCF. You can use another MP2_prog&
                   &. If you insist'
    write(6,'(A)') 'on using PySCF, you need to turn on the RI approximation.'
    stop
   end if
   if((.not.relaxed_dm) .and. TRIM(mp2_prog)=='qchem') then
    write(6,'(/,A)') error_warn//'MP2 unrelaxed density is'
    write(6,'(A)') 'not supported when MP2_prog=QChem. You can use another MP2_&
                   &prog. If you'
    write(6,'(A)') 'insist on using QChem, you can calculate the relaxed densit&
                   &y instead.'
    stop
   end if
  end if
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
  write(6,'(/,A)') 'ERROR in subroutine check_sr_kywd_compatible: RI=noRI. Conf&
                   &used keywords.'
  stop
 end if

 if(force) then
  if(RI .and. TRIM(mp2_prog)=='gamess') then
   write(6,'(/,A)') 'ERROR in subroutine check_sr_kywd_compatible: RI-MP2 analy&
                    &tical gradients'
   write(6,'(A)') 'are not supported in GAMESS.'
   stop
  end if
  if(cc_enabled) then
   select case(TRIM(cc_prog))
   case('orca','gamess','openmolcas')
    write(6,'(/,A)') 'ERROR in subroutine check_sr_kywd_compatible: CC analytic&
                     &al gradients are'
    write(6,'(A)') 'not supported for the current CC_prog. You can use CC_prog=&
                   &Molpro/PySCF/PSI4, etc.'
    stop
   end select
  end if
  if((ccd.or.ccsd_t) .and. TRIM(cc_prog)=='qchem') then
   write(6,'(/,A)') 'ERROR in subroutine check_sr_kywd_compatible: CCD or CCSD(&
                    &T) analytical grad-'
   write(6,'(A)') 'ients are not supported in QChem.'
   stop
  end if
  if(ccsd_t .and. TRIM(cc_prog)=='gaussian') then
   write(6,'(/,A)') 'ERROR in subroutine check_sr_kywd_compatible: CCSD(T) anal&
                    &ytical gradients'
   write(6,'(A)') 'are not supported in Gaussian. You can choose CC_prog=Molpro&
                  &/CFOUR/PSI4.'
   stop
  end if
  if(ccd .and. (.not.RI) .and. TRIM(cc_prog)=='psi4') then
   write(6,'(/,A)') 'ERROR in subroutine check_sr_kywd_compatible: CCD analytic&
                    &al gradients'
   write(6,'(A)') 'in PSI4 should be used with RI.'
   stop
  end if
 end if

 if(ip .or. ea) then
  select case(eom_prog)
  case('gaussian', 'psi4')
   write(6,'(/,A)') 'ERROR in subroutine check_sr_kywd_compatible: IP-/EA-EOM-C&
                    &CSD is not supported'
   write(6,'(A)') 'in Gaussian/PSI4. Please use another EOM_prog.'
   stop
  end select
 end if

 if(RI) call determine_auxbas(basis, RIJK_bas, .true., RIC_bas, F12, F12_cabs)
 call prt_sr_strategy()
end subroutine check_sr_kywd_compatible

! print the internal variables of single reference calculations
subroutine prt_sr_strategy()
 implicit none

 if(DLPNO .and. ccsd_t .and. (.not.iterative_t)) then
  write(6,'(/,A)') REPEAT('-',79)
  write(6,'(A)') 'Remark: DLPNO-CCSD(T) is DLPNO-CCSD(T0) by default in ORCA. I&
                 &f you want higher'
  write(6,'(A)') 'accuracy, it is recommended to use the DLPNO-CCSD(T1) method.'
  write(6,'(A)') REPEAT('-',79)
 end if

 if(.not. RI) then
  write(6,'(/,A)') REPEAT('-',79)
  write(6,'(A)') 'Remark: RI is turned off currently. It is strongly recommende&
                 &d to turn on'
  write(6,'(A)') 'RI if you are running an MP2 job.'
  write(6,'(A)') REPEAT('-',79)
 end if

 write(6,'(/,A)') 'Internal variables:'
 write(6,'(5(A,L1,3X))') 'MP2     = ', mp2, 'CCD     = ', ccd, 'CCSD    = ', &
                          ccsd, 'CCSD(T) = ', ccsd_t, 'ADC     = ', adc
 write(6,'(5(A,L1,3X))') 'EOM     = ', eom, 'RI      = ',  RI, 'F12     = ', &
                          F12,  'DLPNO   = ',  DLPNO, '(T1)    = ', iterative_t
 write(6,'(2(A,L1,3X))') 'IP      = ',  ip, 'EA      = ',  ea
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
  write(6,'(A)')   " Example: autosr h2o.gjf >& h2o.out &"
  write(6,'(A,/)') ' See help: autosr -h'
  stop
 end if

 call getarg(1, fname)

 select case(TRIM(fname))
 case('-v', '-V', '--version')
  write(6,'(A)') 'AutoSR 1.2.6rc39 :: MOKIT, release date: 2024-Sep-14'
  stop
 case('-h','-help','--help')
  write(6,'(/,A)') "Usage: autosr [gjfname] > [outname]"
  write(6,'(A)')   "  Example: autosr h2o.gjf >h2o.out 2>&1 &"
  write(6,'(/,A)') 'Options:'
  write(6,'(A)')   '  -h, -help, --help: Print this message and exit.'
  write(6,'(A)')   '  -v, -V, --version: Print the version number of autosr and exit.'
  write(6,'(A)')   '  -t, --testprog: Print the path of programs detected by autosr and exit.'
  write(6,'(/,A)') 'Methods(#p ...):'
  write(6,'(A)')   '  MP2, RI-MP2, CCD, CCSD, CCSD(T), CCSD-F12, CCSD(T)-F12, C&
                   &CSD(T)-F12a, CCSD(T)-F12b,'
  write(6,'(A)')   '  DLPNO-MP2, DLPNO-CCSD, DLPNO-CCSD(T), DLPNO-CCSD(T0), DLP&
                   &NO-CCSD(T1), EOM-CCSD,'
  write(6,'(A)')   '  EOM-SF-CCSD, EOM-IP-CCSD, EOM-DIP-CCSD, EOM-EA-CCSD, CC2,&
                   & CC3, ADC(2), ADC(3),'
  write(6,'(A)')   '  CIS(D), ROCIS, XCIS, SF-XCIS, SA-SF-CIS'
  write(6,'(/,A)') 'Frequently used keywords in MOKIT{}:'
  write(6,'(A)')   '  HF_prog=Gaussian/PySCF/ORCA/PSI4'
  write(6,'(A)')   ' MP2_prog=Molpro/ORCA/Gaussian/GAMESS/PySCF/Dalton/QChem/CFOUR'
  write(6,'(A)')   '  CC_prog=PySCF/Molpro/CFOUR/ORCA/PSI4/Gaussian/GAMESS/Dalton/QChem'
  write(6,'(A)')   ' CIS_prog=Molpro/ORCA/Gaussian/QChem'
  write(6,'(A)')   ' ADC_prog=Molpro/ORCA/PySCF/QChem'
  write(6,'(A,/)') ' EOM_prog=Molpro/CFOUR/ORCA/Gaussian/GAMESS/PySCF/QChem'
  stop
 case('-t','--testprog')
  call read_sr_program_path()
  stop
 end select

 i = INDEX(fname, '.gjf', back=.true.)
 j = INDEX(fname, '.fch', back=.true.)
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
   write(6,'(A)') REPEAT('-',79)
   write(6,'(A)') 'Warning: you are changing the default number of frozen core &
                  &orbtials.'
   write(6,'(2(A,I0),A)') 'Recommended FC=', chem_core, ', your request=', &
                          core_wish, '.'
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
 use sr_keyword, only: mem, nproc, hf_fch, mo_rhf, bgchg, chgname, ref_e, &
  corr_e, ptchg_e, mp2_e, mp2, mp2_prog, chem_core, ecp_core, customized_core, &
  core_wish, gau_path, psi4_path, gms_path, gms_scr_path, check_gms_path, orca_path,&
  dalton_mpi, gen_no, no_fch, relaxed_dm, force, molcas_omp, RI
 use util_wrapper, only: formchk, unfchk, fch2mkl_wrap, mkl2gbw, bas_fch2py_wrap,&
  fch2com_wrap, fch2psi_wrap, fch2inp_wrap, fch_u2r_wrap, fch2qchem_wrap, &
  fch2cfour_wrap, fch2dal_wrap, fch2inporb_wrap
 use mol, only: natom, grad
 implicit none
 integer :: i, RENAME, SYSTEM
 real(kind=8) :: rtmp
 character(len=24) :: data_string = ' '
 character(len=240) :: old_inp, inpname, chkname, mklname, outname, datname, &
  molname, no_chk, qcscratch

 if(.not. mp2) return
 write(6,'(//,A)') 'Enter subroutine do_mp2...'
 write(6,'(A)') 'MP2 using program '//TRIM(mp2_prog)
 call prt_fc_info(chem_core, ecp_core, customized_core, core_wish)

 if(force .and. chem_core>0) then
  if(TRIM(mp2_prog) == 'pyscf') then
   write(6,'(/,A)') 'ERROR in subroutine do_mp2: MP2 analytical gradients in Py&
                    &SCF cannot be run'
   write(6,'(A)') 'using frozen core. If you still want to use PySCF, you need &
                  &to write FC=0 in mokit{}.'
   stop
  end if
  if((.not.RI) .and. TRIM(mp2_prog)=='psi4') then
   write(6,'(/,A)') 'ERROR in subroutine do_mp2: conventional MP2 analytical gr&
                    &adients in PSI4'
   write(6,'(A)') 'cannot be run using frozen core. If you still want to use PS&
                  &I4, you can either'
   write(6,'(A)') 'write FC=0 in mokit{}, or switch to RI-MP2.'
   stop
  end if
 end if

 if(gen_no .and. (.not.mo_rhf) .and. (.not.relaxed_dm) .and. TRIM(MP2_prog)=='psi4') then
  write(6,'(/,A)') 'ERROR in subroutine do_mp2: UMP2 unrelaxed density is not s&
                   &upported in PSI4'
  write(6,'(A)') 'currently. You can either calculate the relaxed density, or c&
                 &hange another MP2_prog.'
  write(6,'(A)') 'For example, MP2_prog=ORCA.'
  stop
 end if

 i = INDEX(hf_fch, '.fch', back=.true.)
 outname = hf_fch(1:i-1)//'_MP2.out'
 no_fch = hf_fch(1:i-1)//'_MP2_NO.fch'

 select case(TRIM(mp2_prog))
 case('pyscf')
  inpname = hf_fch(1:i-1)//'_MP2.py'
  call bas_fch2py_wrap(hf_fch, .false., inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call prt_posthf_pyscf_inp(inpname, .false.)
  call submit_pyscf_job(inpname, .true.)
  call read_mp2_e_from_pyscf_out(outname, ref_e, mp2_e)
  ref_e = ref_e + ptchg_e
  mp2_e = mp2_e + ptchg_e
 case('gaussian')
  inpname = hf_fch(1:i-1)//'_MP2.gjf'
  chkname = hf_fch(1:i-1)//'_MP2.chk'
  outname = hf_fch(1:i-1)//'_MP2.log'
  call unfchk(hf_fch, chkname)
  call prt_posthf_gau_inp(inpname, .false.)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call submit_gau_job(gau_path, inpname, .true.)
  call read_mp2_e_from_gau_log(outname, ref_e, mp2_e)
  if(gen_no) then
   call formchk(chkname, no_fch)
   if(relaxed_dm) then ! relaxed density
    call copy_dm_and_gen_no(no_fch, hf_fch, 5)
   else                ! unrelaxed density
    call copy_dm_and_gen_no(no_fch, hf_fch,11)
   end if
  end if
  call delete_file(TRIM(chkname))
 case('gamess')
  call check_gms_path()
  datname = hf_fch(1:i-1)//'_MP2.dat'
  inpname = hf_fch(1:i-1)//'_MP2.inp'
  outname = hf_fch(1:i-1)//'_MP2.gms'
  call fch2inp_wrap(hf_fch, .false., 0, 0)
  old_inp = hf_fch(1:i-1)//'.inp'
  i = RENAME(TRIM(old_inp), TRIM(inpname))
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call prt_posthf_gms_inp(inpname, .false.)
  call submit_gms_job(gms_path, gms_scr_path, inpname, nproc)
  call read_mp2_e_from_gms_out(outname, ref_e, mp2_e)
  if(gen_no) then
   call copy_file(hf_fch, no_fch, .false.)
   if(.not. mo_rhf) call fch_u2r_wrap(no_fch, no_fch)
   call dump_gms_no_dat2fch(datname, no_fch)
  end if
 case('orca')
  inpname = hf_fch(1:i-1)//'_MP2.inp'
  old_inp = hf_fch(1:i-1)//'_o.inp'
  chkname = hf_fch(1:i-1)//'_MP2.engrad'
  mklname = hf_fch(1:i-1)//'_MP2.mkl'
  no_chk = hf_fch(1:i-1)//'_MP2.mp2nat'
  call fch2mkl_wrap(hf_fch, mklname)
  i = RENAME(TRIM(old_inp), TRIM(inpname))
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  ! if bgchg = .True., .inp and .mkl file will be updated
  call mkl2gbw(mklname)
  call delete_file(TRIM(mklname))
  call prt_posthf_orca_inp(inpname, .false.)
  call submit_orca_job(orca_path, inpname, .true., .false., .false.)
  call read_posthf_e_from_orca_out(outname, .false., rtmp, ref_e, mp2_e)
  if(gen_no) call dump_orca_no_gbw2fch(no_chk, hf_fch)
 case('molpro')
  inpname = hf_fch(1:i-1)//'_MP2.com'
  mklname = hf_fch(1:i-1)//'_MP2.xml'
  call fch2com_wrap(hf_fch, inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call prt_posthf_molpro_inp(inpname)
  call submit_molpro_job(inpname, mem, nproc)
  call read_mp2_e_from_molpro_out(outname, ref_e, mp2_e)
  ref_e = ref_e + ptchg_e
  mp2_e = mp2_e + ptchg_e
  if(gen_no) then
   call copy_file(hf_fch, no_fch, .false.)
   i = SYSTEM('xml2fch '//TRIM(mklname)//' '//TRIM(no_fch)//' -no')
  end if
 case('psi4')
  inpname = hf_fch(1:i-1)//'_MP2.inp'
  old_inp = hf_fch(1:i-1)//'_MP2.log'
  call fch2psi_wrap(hf_fch, inpname)
  call prt_posthf_psi4_inp(inpname, .false.)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call submit_psi4_job(psi4_path, inpname, nproc)
  call read_mp2_e_from_psi4_out(outname, ref_e, mp2_e)
  ref_e = ref_e + ptchg_e
  mp2_e = mp2_e + ptchg_e
  call delete_file(TRIM(old_inp))
  ! We cannot simply use the .fch file generated by PSI4, because:
  ! 1) The Cartesian coordinates in the PSI4-generated .fch may be translated.
  !    Although the MO coefficients and AO density are not affected by translation,
  !    always using the input coordinates should be adopted. So we need to copy
  !    density into a Gaussian-generated .fch file and generate NOs there.
  ! 2) The density in the PSI4-generated .fch are at MP2, but MOs are at HF, so
  !    we need to generate NOs using the MP2 density.
  ! 3) I hope the post-HF total/spin densities are stored in the 'Total SCF Density'/
  !    'Spin SCF Density' sections, and there is no 'Total MP2 Density'/'Spin MP2 Density'
  !    section.
  ! Considering these things, we call a subroutine to do the job.
  if(gen_no) call copy_dm_and_gen_no(no_fch, hf_fch, 5)
 case('qchem')
  inpname = hf_fch(1:i-1)//'_MP2.in'
  old_inp = hf_fch(1:i-1)//'_MP2.0.FChk'
  call fch2qchem_wrap(hf_fch, 0, inpname)
  call prt_posthf_qchem_inp(inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call submit_qchem_job(inpname, nproc)
  call read_mp2_e_from_qchem_out(outname, ref_e, mp2_e)
  if(force) then
   i = INDEX(hf_fch, '.fch', back=.true.)
   datname = hf_fch(1:i-1)//'_MP2/131.0'
   call getenv('QCSCRATCH', qcscratch)
   call sys_copy_file(TRIM(qcscratch)//'/'//TRIM(datname), '131.0', .false.)
  end if
  if(gen_no) then
   mklname = hf_fch(1:i-1)//'_MP2.FChk'
   i = RENAME(TRIM(mklname), TRIM(no_fch))
   call copy_dm_and_gen_no(no_fch, hf_fch, 1)
  end if
  call delete_file(TRIM(old_inp))
 case('cfour')
  inpname = 'ZMAT'
  call fch2cfour_wrap(hf_fch)
  call prt_posthf_cfour_inp()
  call submit_cfour_job(nproc, outname, .false.)
  call read_mp2_e_from_cfour_out(outname, ref_e, mp2_e)
 case('dalton')
  old_inp = hf_fch(1:i-1)//'_MP2'
  inpname = hf_fch(1:i-1)//'_MP2.dal'
  molname = hf_fch(1:i-1)//'_MP2.mol'
  call fch2dal_wrap(hf_fch, inpname)
  call prt_posthf_dalton_inp(inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(molname))
  call submit_dalton_job(old_inp,mem,nproc,dalton_mpi,.false.,.false.,.false.)
  call read_posthf_e_from_dalton_out(outname, .true., .false., .false., .false.,&
                                     rtmp, ref_e, mp2_e)
 case('openmolcas')
  inpname = hf_fch(1:i-1)//'_MP2.input'
  call fch2inporb_wrap(hf_fch, .false., inpname)
  call prt_posthf_molcas_inp(inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call submit_molcas_job(inpname, mem, nproc, molcas_omp)
  call read_mp2_e_from_molcas_out(outname, ref_e, mp2_e)
  ref_e = ref_e + ptchg_e
  mp2_e = mp2_e + ptchg_e
 case default
  write(6,'(/,A)') 'ERROR in subroutine do_mp2: invalid MP2_prog='//TRIM(mp2_prog)
  stop
 end select

 corr_e = mp2_e - ref_e
 write(6,'(/,A,F18.8,A)')'E(ref)  = ', ref_e, ' a.u.'
 write(6,'(A,F18.8,A)')  'E(corr) = ', corr_e,' a.u.'
 write(6,'(A,F18.8,A)')  'E(MP2)  = ', mp2_e, ' a.u.'

 if(force) then
  ! rename outname to use subroutine read_grad_from_output
  select case(TRIM(mp2_prog))
  case('gamess')
   if(.not. bgchg) outname = datname
  case('orca')
   outname = chkname
  end select
  allocate(grad(3*natom))
  call read_grad_from_output(mp2_prog, outname, natom, grad)
 end if

 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine do_mp2 at '//TRIM(data_string)
end subroutine do_mp2

subroutine do_cc()
 use sr_keyword, only: mem, nproc, ccd, ccsd, ccsd_t, cc_enabled, bgchg, hf_fch,&
  chgname, cc_prog, method, ccd_e, ccsd_e, ccsd_t_e, ref_e, corr_e, ptchg_e, &
  nuc_pt_e, t1diag, chem_core, ecp_core, customized_core, core_wish, RI, gau_path,&
  orca_path, psi4_path, dalton_mpi, gms_path, gms_scr_path, check_gms_path, &
  gen_no, relaxed_dm, no_fch, force, molcas_omp
 use util_wrapper, only: formchk, unfchk, fch2mkl_wrap, mkl2gbw, fch2com_wrap, &
  bas_fch2py_wrap, fch2psi_wrap, fch2inp_wrap, fch2qchem_wrap, fch2cfour_wrap, &
  fch2dal_wrap, fch2inporb_wrap
 use mol, only: natom, grad
 implicit none
 integer :: i, RENAME, SYSTEM
 real(kind=8) :: e = 0d0
 character(len=15) :: method0 = ' '
 character(len=24) :: data_string = ' '
 character(len=240) :: chkname, old_inp, inpname, molname, mklname, outname, &
  no_chk, qcscratch

 if(.not. cc_enabled) return
 write(6,'(//,A)') 'Enter subroutine do_cc...'
 write(6,'(A)') 'CC using program '//TRIM(cc_prog)
 call prt_fc_info(chem_core, ecp_core, customized_core, core_wish)

 if(force .and. chem_core>0) then
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

 if(gen_no .and. (.not.relaxed_dm) .and. RI .and. ccsd_t .and. &
    TRIM(CC_prog)=='psi4') then
  write(6,'(/,A)') 'ERROR in subroutine do_cc: DF-CCSD(T) unrelaxed density is &
                   &not supported in'
  write(6,'(A)') 'PSI4 currently. You can either calculate the relaxed density,&
                 & or change another'
  write(6,'(A)') 'CC_prog. For example, CC_prog=Molpro.'
  stop
 end if

 if((force .or. gen_no) .and. (chem_core>0) .and. (.not.RI) .and. relaxed_dm &
    .and. TRIM(CC_prog)=='psi4') then
  write(6,'(/,A)') 'ERROR in subroutine do_cc: conventional CCSD(T) relaxed den&
                   &sity or analytical'
  write(6,'(A)') 'gradients in PSI4 can only be used when no core orbitals are &
                 &frozen, i.e. FC=0.'
  stop
 end if

 method0 = method
 call strip_ip_ea_eom(method0)
 i = INDEX(hf_fch, '.fch', back=.true.)
 outname = hf_fch(1:i-1)//'_CC.out'
 no_fch = hf_fch(1:i-1)//'_CC_NO.fch'

 select case(TRIM(cc_prog))
 case('pyscf')
  inpname = hf_fch(1:i-1)//'_CC.py'
  call bas_fch2py_wrap(hf_fch, .false., inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call prt_posthf_pyscf_inp(inpname, .false.)
  call submit_pyscf_job(inpname, .true.)
  call read_cc_e_from_pyscf_out(outname, t1diag, ref_e, e)
  ref_e = ref_e + ptchg_e
  e = e + ptchg_e
 case('gaussian')
  inpname = hf_fch(1:i-1)//'_CC.gjf'
  chkname = hf_fch(1:i-1)//'_CC.chk'
  outname = hf_fch(1:i-1)//'_CC.log'
  call unfchk(hf_fch, chkname)
  call prt_posthf_gau_inp(inpname, .false.)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call submit_gau_job(gau_path, inpname, .true.)
  call read_cc_e_from_gau_log(outname, t1diag, ref_e, e)
  if(gen_no) then
   call formchk(chkname, no_fch)
   call copy_dm_and_gen_no(no_fch, hf_fch, 7)
  end if
  call delete_file(TRIM(chkname))
 case('gamess')
  call check_gms_path()
  inpname = hf_fch(1:i-1)//'_CC.inp'
  outname = hf_fch(1:i-1)//'_CC.gms'
  call fch2inp_wrap(hf_fch, .false., 0, 0)
  old_inp = hf_fch(1:i-1)//'.inp'
  i = RENAME(TRIM(old_inp), TRIM(inpname))
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call prt_posthf_gms_inp(inpname, .false.)
  i = nproc
  if(ccd) i = 1 ! CCD in GAMESS is serial
  call submit_gms_job(gms_path, gms_scr_path, inpname, i)
  call read_cc_e_from_gms_out(outname, ccd, ccsd, ccsd_t, t1diag, ref_e, e)
 case('orca')
  inpname = hf_fch(1:i-1)//'_CC.inp'
  old_inp = hf_fch(1:i-1)//'_o.inp'
  mklname = hf_fch(1:i-1)//'_CC.mkl'
  no_chk = hf_fch(1:i-1)//'_CC.mdci.nat'
  call fch2mkl_wrap(hf_fch, mklname)
  i = RENAME(TRIM(old_inp), TRIM(inpname))
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  ! if bgchg = .True., .inp and .mkl file will be updated
  call mkl2gbw(mklname)
  call delete_file(TRIM(mklname))
  call prt_posthf_orca_inp(inpname, .false.)
  call submit_orca_job(orca_path, inpname, .true., .false., .false.)
  call read_posthf_e_from_orca_out(outname, (.not.ccd), t1diag, ref_e, e)
  if(gen_no) call dump_orca_no_gbw2fch(no_chk, hf_fch)
 case('molpro')
  inpname = hf_fch(1:i-1)//'_CC.com'
  mklname = hf_fch(1:i-1)//'_CC.xml'
  call fch2com_wrap(hf_fch, inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call prt_posthf_molpro_inp(inpname)
  call submit_molpro_job(inpname, mem, nproc)
  call read_cc_e_from_molpro_out(outname, (.not.ccd), t1diag, ref_e, e)
  ref_e = ref_e + ptchg_e
  e = e + ptchg_e
  if(gen_no) then
   call copy_file(hf_fch, no_fch, .false.)
   i = SYSTEM('xml2fch '//TRIM(mklname)//' '//TRIM(no_fch)//' -no')
   call update_density_using_no_and_on(no_fch)
  end if
 case('psi4')
  inpname = hf_fch(1:i-1)//'_CC.inp'
  old_inp = hf_fch(1:i-1)//'_CC.log'
  call fch2psi_wrap(hf_fch, inpname)
  call prt_posthf_psi4_inp(inpname, .false.)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call submit_psi4_job(psi4_path, inpname, nproc)
  if(gen_no .and. ccsd_t) then
   call read_ccsd_t_e_from_psi4_grad_out(outname, RI, t1diag, ref_e, e)
  else
   call read_cc_e_from_psi4_out(outname, RI, ccd, ccsd, ccsd_t, t1diag, ref_e, e)
  end if
  ref_e = ref_e + ptchg_e
  e = e + ptchg_e + nuc_pt_e
  if(gen_no) call copy_dm_and_gen_no(no_fch, hf_fch, 7)
  call delete_file(TRIM(old_inp))
 case('qchem')
  inpname = hf_fch(1:i-1)//'_CC.in'
  old_inp = hf_fch(1:i-1)//'_CC.0.FChk'
  call fch2qchem_wrap(hf_fch, 0, inpname)
  call prt_posthf_qchem_inp(inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call submit_qchem_job(inpname, nproc)
  call read_cc_e_from_qchem_out(outname, ccd, ccsd, ccsd_t, t1diag, ref_e, e)
  i = INDEX(hf_fch, '.fch', back=.true.)
  if(force) then
   mklname = hf_fch(1:i-1)//'_CC/131.0'
   call getenv('QCSCRATCH', qcscratch)
   call sys_copy_file(TRIM(qcscratch)//'/'//TRIM(mklname), '131.0', .false.)
  end if
  if(gen_no) then
   mklname = hf_fch(1:i-1)//'_CC.FChk'
   i = RENAME(TRIM(mklname), TRIM(no_fch))
   call copy_dm_and_gen_no(no_fch, hf_fch, 1)
  end if
  call delete_file(TRIM(old_inp))
 case('cfour')
  inpname = 'ZMAT'
  call fch2cfour_wrap(hf_fch)
  call prt_posthf_cfour_inp()
  call submit_cfour_job(nproc, outname, .false.)
  call read_cc_e_from_cfour_out(outname, ccd, ccsd, ccsd_t, t1diag, ref_e, e)
 case('dalton')
  old_inp = hf_fch(1:i-1)//'_CC'
  inpname = hf_fch(1:i-1)//'_CC.dal'
  molname = hf_fch(1:i-1)//'_CC.mol'
  call fch2dal_wrap(hf_fch, inpname)
  call prt_posthf_dalton_inp(inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(molname))
  call submit_dalton_job(old_inp,mem,nproc,dalton_mpi,.false.,.false.,.false.)
  call read_posthf_e_from_dalton_out(outname, .false., ccd, ccsd, ccsd_t, t1diag,&
                                     ref_e, e)
 case('openmolcas')
  inpname = hf_fch(1:i-1)//'_CC.input'
  call fch2inporb_wrap(hf_fch, .false., inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call prt_posthf_molcas_inp(inpname)
  call submit_molcas_job(inpname, mem, nproc, molcas_omp)
  call read_cc_e_from_molcas_out(outname, ccsd_t, t1diag, ref_e, e)
  ref_e = ref_e + ptchg_e
  e = e + ptchg_e
 case default
  write(6,'(/,A)') 'ERROR in subroutine do_cc: invalid CC_prog='//TRIM(cc_prog)
  stop
 end select

 if(ccd) then
  ccd_e = e
 else if(ccsd) then
  ccsd_e = e
 else if(ccsd_t) then
  ccsd_t_e = e
 end if

 corr_e = e - ref_e
 write(6,'(/,A,F18.8,A)')'E(ref)  = ', ref_e, ' a.u.'
 write(6,'(A,F18.8,A)')  'E(corr) = ', corr_e, ' a.u.'
 write(6,'(A,F18.8,A)')  'E('//TRIM(method0)//') = ', e, ' a.u.'
 if(.not. ccd) write(6,'(A,F10.5,A)')  'T1_diag = ', t1diag

 if(force) then
  allocate(grad(3*natom))
  call read_grad_from_output(cc_prog, outname, natom, grad)
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
 use sr_keyword, only: nproc, eom, ip, ea, eom_prog, hf_fch, bgchg, chgname, &
  nstate, gau_path, orca_path, gms_path, gms_scr_path, check_gms_path, &
  psi4_path, ex_elec_e
 use mol, only: ci_mult, fosc
 use util_wrapper, only: bas_fch2py_wrap, unfchk, fch2inp_wrap, fch2mkl_wrap, &
  mkl2gbw, fch2psi_wrap, fch2qchem_wrap
 use phys_cons, only: au2ev
 implicit none
 integer :: i, RENAME, SYSTEM
 real(kind=8), allocatable :: e_ev(:)
 character(len=24) :: data_string = ' '
 character(len=240) :: chkname, inpname, outname, old_inp, mklname

 if(.not. eom) return
 write(6,'(//,A)') 'Enter subroutine do_eomcc...'
 i = INDEX(hf_fch, '.fch', back=.true.)

 if(nstate == 0) nstate = 3
 allocate(ex_elec_e(0:nstate), ci_mult(0:nstate), fosc(nstate))

 select case(TRIM(eom_prog))
 case('pyscf')
  inpname = hf_fch(1:i-1)//'_EOMCC.py'
  outname = hf_fch(1:i-1)//'_EOMCC.out'
  call bas_fch2py_wrap(hf_fch, .false., inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call prt_posthf_pyscf_inp(inpname, .true.)
  call submit_pyscf_job(inpname, .true.)
  call read_eomcc_e_from_pyscf_out(outname, (ip .or. ea), nstate, ex_elec_e, &
                                   ci_mult, fosc)
 case('gaussian')
  inpname = hf_fch(1:i-1)//'_EOMCC.gjf'
  chkname = hf_fch(1:i-1)//'_EOMCC.chk'
  outname = hf_fch(1:i-1)//'_EOMCC.log'
  call unfchk(hf_fch, chkname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call prt_posthf_gau_inp(inpname, .true.)
  call submit_gau_job(gau_path, inpname, .true.)
  call read_eomcc_e_from_gau_log(outname, nstate, ex_elec_e, ci_mult, fosc)
 case('gamess')
  call check_gms_path()
  inpname = hf_fch(1:i-1)//'_EOMCC.inp'
  outname = hf_fch(1:i-1)//'_EOMCC.gms'
  call fch2inp_wrap(hf_fch, .false., 0, 0)
  old_inp = hf_fch(1:i-1)//'.inp'
  i = RENAME(TRIM(old_inp), TRIM(inpname))
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call prt_posthf_gms_inp(inpname, .true.)
  call submit_gms_job(gms_path, gms_scr_path, inpname, 1)
  call read_eomcc_e_from_gms_out(outname, ip, ea, nstate, ex_elec_e, ci_mult, &
                                 fosc)
 case('orca')
  inpname = hf_fch(1:i-1)//'_EOMCC.inp'
  old_inp = hf_fch(1:i-1)//'_o.inp'
  mklname = hf_fch(1:i-1)//'_EOMCC.mkl'
  outname = hf_fch(1:i-1)//'_EOMCC.out'
  call fch2mkl_wrap(hf_fch, mklname)
  i = RENAME(TRIM(old_inp), TRIM(inpname))
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  ! if bgchg = .True., .inp and .mkl file will be updated
  call mkl2gbw(mklname)
  call delete_file(TRIM(mklname))
  call prt_posthf_orca_inp(inpname, .true.)
  call submit_orca_job(orca_path, inpname, .true., .false., .false.)
  call read_eomcc_e_from_orca_out(outname, (ip .or. ea), nstate, ex_elec_e, &
                                  ci_mult, fosc)
 !case('molpro')
 case('psi4')
  inpname = hf_fch(1:i-1)//'_EOMCC.inp'
  outname = hf_fch(1:i-1)//'_EOMCC.out'
  call fch2psi_wrap(hf_fch, inpname)
  call prt_posthf_psi4_inp(inpname, .true.)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call submit_psi4_job(psi4_path, inpname, nproc)
  call read_eomcc_e_from_psi4_out(outname, nstate, ex_elec_e, ci_mult, fosc)
 case('qchem')
  inpname = hf_fch(1:i-1)//'_EOMCC.in'
  outname = hf_fch(1:i-1)//'_EOMCC.out'
  call fch2qchem_wrap(hf_fch, 0, inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call prt_posthf_qchem_inp(inpname)
  stop
 case default
  write(6,'(/,A)') 'ERROR in subroutine do_eomcc: invalid EOM_prog='//TRIM(eom_prog)
  stop
 end select

 allocate(e_ev(nstate), source=0d0)
 forall(i = 1:nstate) e_ev(i) = (ex_elec_e(i) - ex_elec_e(0))*au2ev
 if(ip .or. ea) then
  if(ip) then
   write(6,'(A)') '(0: ground state CCSD, 1: -1e state, >=2: excited states of &
                  &-1e species)'
   write(6,'(54X,A)') 'E_ex/eV'
  else
   write(6,'(A)') '(0: ground state CCSD, 1: +1e state, >=2: excited states of &
                  &+1e species)'
   write(6,'(54X,A)') 'E_ex/eV'
  end if
  write(6,'(A,I3,A,F16.8,A,F6.3)') 'State ',0,', E =',ex_elec_e(0),' a.u. <S**2>&
                                   & =',ci_mult(0)
  do i = 1, nstate, 1
   write(6,'(A,I3,A,F16.8,A,F6.3,2X,F7.3)') 'State ', i, ', E =', ex_elec_e(i),&
                                           ' a.u. <S**2> =', ci_mult(i), e_ev(i)
  end do ! for i
 else
  write(6,'(A)') '(0 for ground state)                                  E_ex/eV &
                 &  fosc'
  write(6,'(A,I3,A,F16.8,A,F6.3)') 'State ',0,', E =',ex_elec_e(0),' a.u. <S**2>&
                                   & =',ci_mult(0)
  do i = 1, nstate, 1
   write(6,'(A,I3,A,F16.8,A,F6.3,2X,F7.3,2X,F7.4)') 'State ', i, ', E =', &
                ex_elec_e(i), ' a.u. <S**2> =', ci_mult(i), e_ev(i), fosc(i)
  end do ! for i
  if(TRIM(eom_prog) == 'pyscf') then
   write(6,'(A)') REPEAT('-', 79)
   write(6,'(A)') 'Warning from subroutine do_eomcc: oscillator strengths of EE&
                  &-EOM-CCSD cannot be'
   write(6,'(A)') 'calculated using EOM_prog=PySCF. So the oscillator strengths&
                  & above are all'
   write(6,'(A)') 'printed as zero. If you care about this, please use another &
                  &EOM_prog, e.g. ORCA.'
   write(6,'(A)') REPEAT('-', 79)
  end if
 end if

 deallocate(e_ev)
 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine do_eomcc at '//TRIM(data_string)
end subroutine do_eomcc

! print a Gaussian MP2/CCSD/CCSD/CCSD(T) gjf file
subroutine prt_posthf_gau_inp(gjfname, excited)
 use sr_keyword, only: mem, nproc, DKH2, mp2, ccd, ccsd, ccsd_t, chem_core, &
  eom, nstate, gen_no, relaxed_dm, force
 implicit none
 integer :: i, fid
 character(len=240), intent(in) :: gjfname
 logical, intent(in) :: excited

 i = INDEX(gjfname, '.gjf', back=.true.)
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
  if(excited) then
   if(eom) write(fid,'(A,I0,A)',advance='no') 'EOMCCSD(nstates=',nstate,')'
  else
   write(fid,'(A)',advance='no') 'CCSD(T1diag)'
  end if
 else if(ccsd_t) then
  write(fid,'(A)',advance='no') 'CCSD(T1diag,T)'
 end if

 write(fid,'(A)',advance='no') ' chkbasis nosymm guess=read geom=allcheck'

 if(DKH2) then
  write(fid,'(A)',advance='no') ' int(nobasistransform,DKH2) iop(3/93=1)'
 else
  write(fid,'(A)',advance='no') ' int=nobasistransform'
 end if

 write(fid,'(A,I0,A)',advance='no') ' window=(', chem_core+1 ,',0)'
 if(force) write(fid,'(A)',advance='no') ' force'

 if(gen_no) then
  write(fid,'(A)',advance='no') ' density'
  if(.not. relaxed_dm) write(fid,'(A)',advance='no') '=Rho2 trans=conventional'
 end if

 write(fid,'(/)',advance='no')
 close(fid)
end subroutine prt_posthf_gau_inp

! add CC keywords into a ORCA input file
subroutine prt_posthf_orca_inp(inpname, excited)
 use sr_keyword, only: mem, nproc, RI, DLPNO, F12, RIJK_bas, RIC_bas, F12_cabs,&
  mo_rhf, lin_dep, mp2, ccd, ccsd, ccsd_t, iterative_t, gen_no, relaxed_dm, ip,&
  ea, nstate, chem_core, force
 use mol, only: mult
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname
 logical, intent(in) :: excited

 i = INDEX(inpname, '.inp', back=.true.)
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
  if(TRIM(RIJK_bas) == 'autoaux') then
   write(fid1,'(A)',advance='no') ' RIJK AutoAux'
  else
   write(fid1,'(A)',advance='no') ' RIJK '//TRIM(RIJK_bas)//' '//TRIM(RIC_bas)
  end if
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
  else if(ccsd) then
   write(fid1,'(A)',advance='no') 'CCSD'
  else if(mp2) then
   write(fid1,'(A)',advance='no') 'MP2'
  end if
 else
  if(RI) then
   write(fid1,'(A)',advance='no') ' RI-'
  else
   write(fid1,'(A)',advance='no') ' '
  end if
  if(ccsd_t) then
   write(fid1,'(A)',advance='no') 'CCSD(T)'
  else if(ccsd) then
   if(excited) then
    if(ip) then
     write(fid1,'(A)',advance='no') 'IP-EOM-CCSD'
    else if(ea) then
     write(fid1,'(A)',advance='no') 'EA-EOM-CCSD'
    else
     write(fid1,'(A)',advance='no') 'EOM-CCSD'
    end if
   else
    write(fid1,'(A)',advance='no') 'CCSD'
   end if
  else if(ccd) then
   write(fid1,'(A)',advance='no') 'CCD'
  else if(mp2) then
   write(fid1,'(A)',advance='no') 'MP2'
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

 ! keep the integral accuracy the same as Gaussian
 write(fid1,'(A)') '%scf'
 write(fid1,'(A)') ' Thresh 1e-12'
 write(fid1,'(A)') ' Tcut 1e-14'
 write(fid1,'(A)') ' CNVDamp False'
 if(lin_dep) write(fid1,'(A)') ' sthresh 1e-6'
 write(fid1,'(A)') 'end'

 if(.not. mp2) then
  write(fid1,'(A)') '%mdci'
  write(fid1,'(A)') ' MaxIter 300'
  if(excited) then
   write(fid1,'(A,I0)') ' nroots ',nstate
   if(.not. (ip .or. ea)) then
    write(fid1,'(A)') ' DoLeft True'
    write(fid1,'(A)') ' DoTDM True'
   end if
  end if
  if(gen_no .and. (ccd .or. ccsd)) then
   write(fid1,'(A)') ' density unrelaxed'
   write(fid1,'(A)') ' NatOrbs True'
  end if
  write(fid1,'(A)') 'end'
 end if

 if(gen_no .and. mp2) then
  write(fid1,'(A)') '%mp2'
  if(relaxed_dm) then
   write(fid1,'(A)') ' density relaxed'
  else
   write(fid1,'(A)') ' density unrelaxed'
  end if
  write(fid1,'(A)') ' NatOrbs True'
  write(fid1,'(A)') 'end'
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
 use sr_keyword, only: mp2, ccd, ccsd_t, RI, F12, mo_rhf, basis, RIJK_bas, &
  RIC_bas, chem_core, gen_no, relaxed_dm, force
 use mol, only: mult
 implicit none
 integer :: i, fid
 character(len=21) :: MP2FIT_bas, RIJK_bas1
 character(len=240) :: wfuname
 character(len=240), intent(in) :: inpname
 logical :: alive

 if((.not.mo_rhf) .and. mp2 .and. gen_no) then
  write(6,'(/,A)') 'ERROR in subroutine prt_posthf_molpro_inp: UMP2 NOs are not&
                  & supported in'
  write(6,'(A)') 'Molpro. Please change another MP2_prog.'
  stop
 end if

 i = INDEX(inpname, '.com', back=.true.)
 wfuname = inpname(1:i-1)//'.wfu'
 call lower(wfuname)
 inquire(file=TRIM(wfuname),exist=alive)

 if(RI) then
  call auxbas_convert(RIJK_bas, RIJK_bas1, 2)
  i = INDEX(RIJK_bas, '/')
  RIC_bas = RIJK_bas(1:i)//'optri'
  MP2FIT_bas = RIJK_bas(1:i)//'mp2fit'
 end if

 open(newunit=fid,file=TRIM(inpname),status='old',position='append')
 BACKSPACE(fid)

 if(.not. mp2) then
  if(alive) then
   write(fid,'(A)') 'file,3,'//TRIM(wfuname)
  else
   write(fid,'(A)') 'file,3,'//TRIM(wfuname)//',new'
  end if
 end if

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
  if(mo_rhf) then
    if(mult /= 1) write(fid,'(A)',advance='no') 'U' ! ROHF-UCCSD
    ! ROHF-ROCCSD interface is currently not supported in autosr
  else
   write(6,'(/,A)') 'ERROR in subroutine prt_posthf_molpro_inp: UHF-based CC i&
                    &s not supported in Molpro.'
   close(fid)
   stop
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
  write(fid,'(A)',advance='no') '{'
  if(.not. mo_rhf) write(fid,'(A)',advance='no') 'U'
  if(mp2) then
   write(fid,'(A)',advance='no') 'MP2'
  else
   write(fid,'(A)',advance='no') 'CCSD'
  end if
  if(ccsd_t) write(fid,'(A)',advance='no') '(T)'
  if(F12) write(fid,'(A)',advance='no') '-F12'
 end if

 if(.not. mp2) write(fid,'(A)',advance='no') ',thrden=1d-7'
 write(fid,'(A,I0,A)',advance='no') ';core,',chem_core
 if(ccd) write(fid,'(A)',advance='no') ';NoSING'
 if(gen_no) then
  if(.not. relaxed_dm) write(fid,'(A)',advance='no') ';expec,norelax'
  write(fid,'(A)',advance='no') ';natorb,5200.2'
 end if

 if(mp2) then
  write(fid,'(A)') '}'
 else
  if(alive) then
   write(fid,'(A)') ';start,4000.3}'
  else
   write(fid,'(A)') ';save,4000.3}'
  end if
 end if

 if(gen_no) write(fid,'(A)') '{put,xml;orbital,5200.2;keepspherical}'
 if(force) write(fid,'(A)') 'Force'
 close(fid)
end subroutine prt_posthf_molpro_inp

subroutine prt_posthf_pyscf_inp(pyname, excited)
 use sr_keyword, only: mem, nproc, mp2, cc_enabled, ccd, ccsd_t, chem_core, hf_fch,&
  force, gen_no, relaxed_dm, mo_rhf, RI, RIJK_bas, RIC_bas, ip, ea, nstate
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=21) :: RIJK_bas1, RIC_bas1
 character(len=240) :: buf, pyname1, no_fch, mo_npy, t1_npy, t2_npy
 character(len=240), intent(in) :: pyname
 logical, intent(in) :: excited

 i = INDEX(pyname, '.py', back=.true.)
 pyname1 = pyname(1:i-1)//'.t'
 no_fch = pyname(1:i-1)//'_NO.fch'
 mo_npy = pyname(1:i-1)//'_mo.npy'
 t1_npy = pyname(1:i-1)//'_t1.npy'
 t2_npy = pyname(1:i-1)//'_t2.npy'

 open(newunit=fid,file=TRIM(pyname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(pyname1),status='replace')

 if(mp2) then
  write(fid1,'(A)') 'from pyscf import mp'
  if(RI) then
   if(mo_rhf) then
    write(fid1,'(A)') 'from pyscf.mp.dfmp2_native import DFMP2'
   else
    write(fid1,'(A)') 'from pyscf.mp.dfump2_native import DFMP2'
   end if
  end if
 else
  write(fid1,'(A)',advance='no') 'from pyscf import '
  if(RI) then
   write(fid1,'(A)') 'cc, df'
  else
   write(fid1,'(A)') 'cc'
  end if
  if(ccsd_t) write(fid1,'(A)') 'from pyscf.cc import ccsd_t'
 end if

 if(gen_no) then
  write(fid1,'(A)') 'from shutil import copyfile'
  write(fid1,'(A)') 'from mokit.lib.py2fch import write_pyscf_dm_into_fch'
  write(fid1,'(A)') 'from mokit.lib.lo import gen_no_using_density_in_fch'
  if(ccsd_t) then
   write(fid1,'(A)') 'from pyscf.cc import ccsd_t_lambda_slow as ccsd_t_lambda'
   write(fid1,'(A)') 'from pyscf.cc import ccsd_t_rdm_slow as ccsd_t_rdm'
  end if
 end if
 if(.not. mp2) then
  write(fid1,'(A,/,A)') 'import numpy as np', 'import os'
  write(fid1,'(A)') 'from mokit.lib.mo_svd import update_amp_from_mo'
 end if

 do while(.true.)
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 write(fid1,'(/,A,I0)') 'nproc = ', nproc
 write(fid1,'(A)') 'lib.num_threads(nproc)'

 read(fid,'(A)') buf
 if(buf(1:15) /= 'lib.num_threads') write(fid1,'(A)') TRIM(buf)

 if(RI) then
  do while(.true.)
   read(fid,'(A)') buf
   if(buf(1:2) == 'mf') exit
   write(fid1,'(A)') TRIM(buf)
  end do ! for while

  call auxbas_convert(RIJK_bas, RIJK_bas1, 1)
  write(fid1,'(A)') TRIM(buf)//".density_fit(auxbasis='"//TRIM(RIJK_bas1)//"')"
 end if

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:3) == '#dm') exit
  if(buf(1:13) == 'mf.max_memory') cycle
  if(buf(1:9) == 'mf.kernel') then
   write(fid1,'(A,I0,A)') 'mf.max_memory = ',mem*1000,' # MB'
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
   write(fid1,'(/,A)') 'mc = DFMP2(mf)'
   write(fid1,'(A)') 'mc.max_memory = mf.max_memory'
  else
   write(fid1,'(/,A)') 'mc = mp.MP2(mf)'
  end if
 else
  write(fid1,'(/,A)') 'mc = cc.CCSD(mf)'
 end if

 write(fid1,'(A)') 'mc.verbose = 4'
 write(fid1,'(A,I0)') 'mc.frozen = ', chem_core
 if(RI .and. cc_enabled) then
  i = LEN_TRIM(RIC_bas)
  RIC_bas1 = RIC_bas(1:i-2)//'-ri'
  write(fid1,'(A)') "mc.with_df = df.DF(mol, auxbasis='"//TRIM(RIC_bas1)//"')"
 end if

 if(ccd) then
  write(fid1,'(/,A)') 'old_update_amps = mc.update_amps'
  write(fid1,'(A)') 'def update_amps(t1, t2, eris):'
  write(fid1,'(A)') '  t1, t2 = old_update_amps(t1, t2, eris)'
  write(fid1,'(A)') '  return t1*0, t2'
  write(fid1,'(A)') 'mc.update_amps = update_amps'
 end if

 if(mp2) then
  write(fid1,'(A)') 'mc.kernel()'
 else
  write(fid1,'(/,A)') "alive1 = os.path.exists('"//TRIM(t1_npy)//"')"
  write(fid1,'(A)') "alive2 = os.path.exists('"//TRIM(t2_npy)//"')"
  write(fid1,'(A)') 'if alive1 and alive2:'
  write(fid1,'(A)') '  nocc = np.count_nonzero(mf.mo_occ > 0)'
  write(fid1,'(A)') '  nvir = nif - nocc'
  write(fid1,'(A)') '  nocc = nocc - mc.frozen'
  write(fid1,'(A)') "  old_mo = np.load('"//TRIM(mo_npy)//"')"
  write(fid1,'(A)') "  t1 = np.load('"//TRIM(t1_npy)//"')"
  write(fid1,'(A)') "  t2 = np.load('"//TRIM(t2_npy)//"')"
  write(fid1,'(A)') '  new_t1, new_t2 = update_amp_from_mo(nbf,nif,nocc,nvir,ol&
                    &d_mo,mc.mo_coeff,t1,t2)'
  write(fid1,'(A)') '  mc.kernel(t1=new_t1, t2=new_t2)'
  write(fid1,'(A)') 'else:'
  write(fid1,'(A)') '  mc.kernel()'
  write(fid1,'(A)') "np.save('"//TRIM(mo_npy)//"', mc.mo_coeff)"
  write(fid1,'(A)') "np.save('"//TRIM(t1_npy)//"', mc.t1)"
  write(fid1,'(A)') "np.save('"//TRIM(t2_npy)//"', mc.t2)"
 end if

 if((.not.mp2) .and. (.not.ccd)) then
  write(fid1,'(/,A)') 'T1diag = mc.get_t1_diagnostic()'
  write(fid1,'(A)') "print('T1_diag =', T1diag)"
 end if

 if(ccsd_t) write(fid1,'(/,A)') 'mc.ccsd_t()'

 if(excited) then
  if(ip) then      ! IP-EOM-CCSD
   write(fid1,'(A,I0,A)') 'mc.ipccsd(nroots=', nstate, ')'
  else if(ea) then ! EA-EOM-CCSD
   write(fid1,'(A,I0,A)') 'mc.eaccsd(nroots=', nstate, ')'
  else             ! EE-EOM-CCSD
   write(fid1,'(A,I0,A)') 'mc.eomee_ccsd_singlet(nroots=', nstate, ')'
  end if
 end if

 if(gen_no) then
  write(fid1,'(/,A)') "copyfile('"//TRIM(hf_fch)//"', '"//TRIM(no_fch)//"')"
  if(RI .and. mp2) then ! DF-MP2
   if(relaxed_dm) then
    write(fid1,'(A)') 'dm1 = mc.make_rdm1_relaxed(ao_repr=True)'
   else
    write(fid1,'(A)') 'dm1 = mc.make_rdm1_unrelaxed(ao_repr=True)'
   end if
  else                  ! not DF-MP2
   if(ccsd_t) then
    write(fid1,'(A)') 'eris = mc.ao2mo()'
    write(fid1,'(A)') 'conv, l1, l2 = ccsd_t_lambda.kernel(mc, eris, mc.t1, mc.t2)'
    write(fid1,'(A)') 'dm1 = ccsd_t_rdm.make_rdm1(mc, mc.t1, mc.t2, l1, l2, eri&
                      &s=eris, ao_repr=True)'
   else
    write(fid1,'(A)') 'dm1 = mc.make_rdm1(ao_repr=True)'
   end if
  end if
  write(fid1,'(A)') "write_pyscf_dm_into_fch('"//TRIM(no_fch)//"', nbf, dm1, 1,&
                    & True)"
  write(fid1,'(A)') "gen_no_using_density_in_fch('"//TRIM(no_fch)//"', 1)"
 end if

 close(fid1)
 i = RENAME(TRIM(pyname1), TRIM(pyname))
 if(force) call add_force_key2py_script(mem, pyname, ccsd_t)
end subroutine prt_posthf_pyscf_inp

subroutine prt_posthf_psi4_inp(inpname, excited)
 use sr_keyword, only: mem, mp2, ccd, ccsd, ccsd_t, cc_enabled, chem_core, RI, &
  RIJK_bas, RIC_bas, gen_no, relaxed_dm, force, no_fch, nstate
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=10) :: psi4_ver
 character(len=21) :: RIJK_bas1, RIC_bas1
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname
 logical, intent(in) :: excited

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
  call auxbas_convert(RIJK_bas, RIJK_bas1, 1)
  write(fid1,'(A)') ' df_basis_scf '//TRIM(RIJK_bas1)
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
 write(fid1,'(A,I0)') 'set num_frozen_docc ', chem_core
 if(gen_no .and. (.not.relaxed_dm)) write(fid1,'(A)') 'set opdm_relax False'

 if(RI) then
  i = LEN_TRIM(RIC_bas)
  RIC_bas1 = RIC_bas(1:i-2)//'-ri'
  if(cc_enabled) then
   write(fid1,'(A)') 'set cc_type df'
   write(fid1,'(A)') 'set df_basis_cc '//TRIM(RIC_bas1)
  end if
  if(mp2) write(fid1,'(A)') 'set df_basis_mp2 '//TRIM(RIC_bas1)
 else
  if(mp2) write(fid1,'(A)') 'set mp2_type conv'
  if(cc_enabled) write(fid1,'(A)') 'set cc_type conv'
 end if
 if(excited) write(fid1,'(A,I0)') 'set roots_per_irrep ', nstate

 if(ccsd_t .and. (.not.RI)) call get_psi4_version(psi4_ver)

 if(force) then
  ! PSI4 >=1.8, extra keywords are needed to calculate CCSD(T) analytical gradients
  if(ccsd_t .and. (.not.RI) .and. psi4_ver(1:3)=='1.8') then
   write(fid1,'(A)') 'set qc_module ccenergy'
  end if
  write(fid1,'(A)',advance='no') "gradient('"
 else
  if(excited) then
   write(fid1,'(A)',advance='no') "properties('"
  else
   if(gen_no) then
    if(ccsd_t .and. (.not.RI) .and. psi4_ver(1:3)=='1.8') then
     write(fid1,'(A)') 'set qc_module ccenergy'
    end if
    write(fid1,'(A)',advance='no') "grad, wfn = gradient('"
   else
    write(fid1,'(A)',advance='no') "energy('"
   end if
  end if
 end if

 if(excited) then
  write(fid1,'(A)') "eom-ccsd',properties=['oscillator_strength'])"
 else
  if(mp2) then
   write(fid1,'(A)',advance='no') "mp2'"
  else if(ccd) then
   write(fid1,'(A)',advance='no') "ccd'"
  else if(ccsd) then
   write(fid1,'(A)',advance='no') "ccsd'"
  else if(ccsd_t) then
   write(fid1,'(A)',advance='no') "ccsd(t)'"
  end if
  if(gen_no) write(fid1,'(A)',advance='no') ', return_wfn=True'
  write(fid1,'(A)') ')'
 end if

 if(gen_no) then
  write(fid1,'(A)') "fchk(wfn, '"//TRIM(no_fch)//"')"
 end if

 if(.not. (mp2 .or. ccd .or. excited)) write(fid1,'(A)') "print_variables()"
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine prt_posthf_psi4_inp

subroutine prt_posthf_gms_inp(inpname, excited)
 use sr_keyword, only: mem, nproc, mp2, ccd, ccsd, ccsd_t, cc_enabled, force, &
  chem_core, RI, mo_rhf, gen_no, ip, ea, nstate
 implicit none
 integer :: i, j, fid, fid1, RENAME
 character(len=14) :: method
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname
 logical, intent(in) :: excited

 if(mp2) then
  method = 'MPLEVL=2'
 else if(ccd) then
  method = 'CCTYP=CCD'
 else if(ccsd) then
  if(excited) then
   if(ip) then
    method = 'CCTYP=IP-EOM2'
   else if(ea) then
    method = 'CCTYP=EA-EOM2'
   else
    method = 'CCTYP=EOM-CCSD'
   end if
  else ! ground state
   method = 'CCTYP=CCSD'
  end if
 else if(ccsd_t) then
  method = 'CCTYP=CCSD(T)'
 end if

 inpname1 = TRIM(inpname)//'.t'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 read(fid,'(A)') buf
 if(force) then
  i = INDEX(buf, 'ENERGY')
  buf = buf(1:i-1)//'GRADIENT'//TRIM(buf(i+6:))
 end if
 write(fid1,'(A)') TRIM(buf)

 read(fid,'(A)') buf
 i = INDEX(buf, '$END')
 write(fid1,'(A)') buf(1:i-1)//TRIM(method)//' $END'

 read(fid,'(A)') buf ! assuming this is $SYSTEM
 i = FLOOR(DBLE(mem)*475d0/(4d0*DBLE(nproc)))
 j = FLOOR(DBLE(mem)*25d0/4d0)
 write(fid1,'(2(A,I0),A)') " $SYSTEM MWORDS=",i,' MEMDDI=',j," $END"

 if(mp2) then
  write(fid1,'(A,I0)',advance='no') " $MP2 NACORE=", chem_core
  if(.not. mo_rhf) write(fid1,'(A,I0)',advance='no') " NBCORE=", chem_core
  if(gen_no) write(fid1,'(A)',advance='no') ' MP2PRP=.T.'
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

 if(excited) then
  if(ip .or. ea) then
   write(fid1,'(A,I0,A)') ' $EOMINP NSTATE(1)=', nstate, ' $END'
  else
   write(fid1,'(A,I0,A)') ' $EOMINP CCPRPE=.T. NSTATE(1)=', nstate, ' $END'
  end if
 end if

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

! print post-HF keywords into a Q-Chem input(.in) file
subroutine prt_posthf_qchem_inp(inpname)
 use sr_keyword, only: mem, mp2, ccd, ccsd, ccsd_t, cc_enabled, chem_core, RI, &
  force, gen_no, relaxed_dm, basis
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 inpname1 = TRIM(inpname)//'.t'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:6) == 'method') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(mp2) then
  if(RI) then
   write(fid1,'(A)') 'method RIMP2'
   write(fid1,'(A)') 'aux_basis rimp2-'//TRIM(basis)
  else
   write(fid1,'(A)') 'method MP2'
  end if
  if(gen_no .and. (.not.force)) write(fid1,'(A)') 'jobtype force'
 else if(cc_enabled) then
  if(ccd) then
   write(fid1,'(A)') 'method CCD'
  else if(ccsd) then
   write(fid1,'(A)') 'method CCSD'
  else if(ccsd_t) then
   write(fid1,'(A)') 'method CCSD(T)'
  end if
  if(gen_no) then
   write(fid1,'(A)') 'cc_ref_prop true'
   if(relaxed_dm) write(fid1,'(A)') 'cc_fullresponse true'
  end if
 end if
 if(force) write(fid1,'(A)') 'jobtype force'

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:7) == 'mem_tot') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 write(fid1,'(A,I0)') 'mem_total ', mem*1000 ! MB
 i = min(2000, mem*200)
 write(fid1,'(A,I0)') 'mem_static ', i ! MB
 if(cc_enabled) write(fid1,'(A,I0)') 'cc_memory ', mem*1000-i
 write(fid1,'(A,I0)') 'n_frozen_core ', chem_core

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine prt_posthf_qchem_inp

! print post-HF keywords into a Q-Chem input(.in) file
subroutine prt_posthf_cfour_inp()
 use sr_keyword, only: mem, mo_rhf, mp2, ccd, ccsd, ccsd_t, cc_enabled, force, &
  chem_core
 use mol, only: mult
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=240) :: buf

 open(newunit=fid,file='ZMAT',status='old',position='rewind')
 open(newunit=fid1,file='ZMAT1',status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:6) == '*CFOUR') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine prt_posthf_cfour_inp: keyword '*CFOUR' &
                   &not found in file ZMAT."
  close(fid1,status='delete')
  close(fid)
  stop
 end if

 write(fid1,'(A)',advance='no') '*CFOUR(CALC='

 if(mp2) then
  write(fid1,'(A)',advance='no') 'MP2'
 else if(cc_enabled) then
  if(ccd) then
   write(fid1,'(A)',advance='no') 'CCD'
  else if(ccsd) then
   write(fid1,'(A)',advance='no') 'CCSD'
  else if(ccsd_t) then
   write(fid1,'(A)',advance='no') 'CCSD(T)'
  end if
  write(fid1,'(A)',advance='no') ',CC_PROGRAM='
  if(mult==1 .and. mo_rhf) then
   write(fid1,'(A)',advance='no') 'NCC'
  else
   write(fid1,'(A)',advance='no') 'ECC'
  endif
 end if

 i = INDEX(buf, ',')
 write(fid1,'(A)') TRIM(buf(i:))

 if(.not. ccd) write(fid1,'(A)') 'ABCDTYPE=AOBASIS'

 if(chem_core > 0) then
  if(chem_core == 1) then
   write(fid1,'(A)') 'DROPMO=1'
  else
   write(fid1,'(A,I0)') 'DROPMO=1>', chem_core
  end if
 end if

 write(fid1,'(A,I0,A)') 'MEMORY=', mem ,',MEM_UNIT=GB'
 if(force) write(fid1,'(A)') 'DERIV_LEVEL=1,FCGRADNEW=NEW'

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME('ZMAT1', 'ZMAT')
end subroutine prt_posthf_cfour_inp

! print post-HF keywords into a Dalton input file (.dal)
subroutine prt_posthf_dalton_inp(inpname)
 use sr_keyword, only: mp2, ccd, ccsd, ccsd_t, cc_enabled, chem_core, force
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 inpname1 = TRIM(inpname)//'.t'
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 write(fid1,'(A,/,A)') '**DALTON INPUT', '.RUN WAVE FUNCTIONS'
 if(force) write(fid1,'(A,/,A,/,A,/,A)') '**INTEGRAL','.DIPLEN','.DEROVL','.DERHAM'
 write(fid1,'(A)') '**WAVE FUNCTIONS'

 write(fid1,'(A,/,A)') '.CC', '*CC INPUT'
 if(mp2) then
  write(fid1,'(A)') '.MP2'
 else if(cc_enabled) then
  if(ccd) then
   write(fid1,'(A)') '.CCD'
  else if(ccsd) then
   write(fid1,'(A)') '.CCSD'
  else if(ccsd_t) then
   write(fid1,'(A)') '.CC(T)'
  end if
 end if
 if(chem_core > 0) write(fid1,'(A,/,I0,A)') '.FREEZE', chem_core, ' 0'

 if(force) write(fid1,'(A)') '*DERIVATIVES'

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:8) == '*ORBITAL') exit
 end do ! for while

 write(fid1,'(A)') '*ORBITAL INPUT'
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid)
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine prt_posthf_dalton_inp

! print post-HF keywords into a OpenMolcas input file (.input)
subroutine prt_posthf_molcas_inp(inpname)
 use sr_keyword, only: mp2, ccsd, ccsd_t, cc_enabled, chem_core, force, RI
 use mol, only: charge, mult
 implicit none
 integer :: i, fid
 character(len=240) :: buf, orbname
 character(len=240), intent(in) :: inpname

 if(RI) call add_RI_kywd_into_molcas_inp(inpname, .true.)

 call find_specified_suffix(inpname, '.input', i)
 orbname = inpname(1:i-1)//'.INPORB'
 open(newunit=fid,file=TRIM(inpname),status='old',position='append')

 ! the RASSCF module is required to perform SCF for conventional CC
 if((.not.RI) .and. cc_enabled) then
  do while(.true.)
   BACKSPACE(fid)
   BACKSPACE(fid)
   read(fid,'(A)') buf
   if(buf(1:4)=="&SCF" .or. buf(1:7)=="&RASSCF") exit
  end do ! for while
  BACKSPACE(fid)
  write(fid,'(A)') "&RASSCF"
  write(fid,'(A,I0)') 'Spin= ', mult
  write(fid,'(A,I0)') 'Charge= ', charge
  write(fid,'(A,I0,/,A,I0)') 'nActEl= ',mult-1, 'RAS2= ',mult-1
  write(fid,'(A)') 'FILEORB= '//TRIM(orbname)
  write(fid,'(A)') 'OutOrbitals= Canonical'
  write(fid,'(/,A,/,A)') "&MOTRA", 'JobIph'
  write(fid,'(A,I0)') 'Frozen= ', chem_core
 end if

 if(mp2) then
  write(fid,'(/,A)') "&MBPT2"
  write(fid,'(A,I0)') 'Frozen= ', chem_core
  if(force) write(fid,'(A,//,A)') 'GRDT', "&ALASKA"
 else if(cc_enabled) then
  if(RI) then
   write(fid,'(/,A,/,A)') "&CHCC", 'MAXITERATIONS= 300'
   write(fid,'(A,I0)') 'Frozen= ', chem_core
   if(ccsd_t) write(fid,'(/,A,/,A,I0)') "&CHT3", 'Frozen= ', chem_core
  else
   write(fid,'(/,A)') "&CCSDT"
   if(ccsd) then
    write(fid,'(A)') 'CCSD'
   else if(ccsd_t) then
    write(fid,'(A,/,A)') 'CCT','TRIPLES= 2'
   end if
   write(fid,'(A)') 'ITERATIONS= 300'
  end if
 end if

 write(fid,'(/)',advance='no')
 close(fid)
end subroutine prt_posthf_molcas_inp

subroutine read_mp2_e_from_gau_log(outname, ref_e, tot_e)
 implicit none
 integer :: i, fid
 real(kind=8) :: ss
 real(kind=8), intent(out) :: ref_e, tot_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 ref_e = 0d0; tot_e = 0d0
 call read_hf_e_and_ss_from_gau_log(outname, ref_e, ss)

 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(28:34)=='EUMP2 =' .or. buf(35:40)=='EUMP2=') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_mp2_e_from_gau_log: 'EUMP2 =' not &
                   &found in"
  write(6,'(A)') 'file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, '=', back=.true.)
 read(buf(i+1:),*) tot_e ! MP2 total energy
 close(fid)
end subroutine read_mp2_e_from_gau_log

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
  if(INDEX(buf,'MP2') > 0) exit 
 end do ! for while

 read(fid,*) tot_e, ref_e
 close(fid)
end subroutine read_mp2_e_from_molpro_out

subroutine read_mp2_e_from_pyscf_out(outname, ref_e, tot_e)
 implicit none
 integer :: i, icase, fid
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

  if(buf(1:8)=='E(MP2) =') then
   icase = 1
   exit
  else if(buf(1:10)=='DF-MP2 cor') then
   icase = 2
   exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_mp2_e_from_pyscf_out: MP2 energy k&
                   &eywords not found'
  write(6,'(A)') 'in file '//TRIM(outname)
  close(fid)
  stop
 end if

 select case(icase)
 case(1) ! MP2
  i = INDEX(buf, '=')
  read(buf(i+1:),*) tot_e
  i = INDEX(buf, '=', back=.true.)
  read(buf(i+1:),*) ref_e ! here is MP2 correlation energy
  ref_e = tot_e - ref_e
 case(2) ! DF-MP2
  i = INDEX(buf, ':')
  read(buf(i+1:),*) tot_e ! here is DF-MP2 correlation energy

  do while(.true.)
   BACKSPACE(fid)
   BACKSPACE(fid)
   read(fid,'(A)') buf
   if(buf(1:15) == 'converged SCF e') exit
  end do ! for while
  i = INDEX(buf, '=')
  read(buf(i+1:),*) ref_e
  tot_e = tot_e + ref_e
 case default
  write(6,'(/,A)') 'ERROR in subroutine read_mp2_e_from_pyscf_out: icase out of&
                   & range.'
  stop
 end select

 close(fid)
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

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_scf_e_from_gms_out: ' FINAL ' not &
                   &found in file "//TRIM(outname)
  stop
 end if

 i = INDEX(buf, 'IS')
 buf = ADJUSTL(buf(i+2:))
 read(buf,*) scf_e
end subroutine read_scf_e_from_gms_out

subroutine read_mp2_e_from_gms_out(outname, ref_e, tot_e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: ref_e, tot_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 call read_scf_e_from_gms_out(outname, ref_e)
 tot_e = 0d0

 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(INDEX(buf,'E(MP2)=') > 0) exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_mp2_e_from_gms_out: no 'E(MP2)=' f&
                   &ound in file "//TRIM(outname)
  stop
 end if

 i = INDEX(buf,'E(MP2)=')
 read(buf(i+7:),*) tot_e
end subroutine read_mp2_e_from_gms_out

subroutine read_mp2_e_from_psi4_out(outname, ref_e, mp2_e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: ref_e, mp2_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical :: ri

 ri = .false.
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(INDEX(buf,'REF Energy') > 0) exit
  if(INDEX(buf,'Reference Energy') > 0) then
   ri = .true.
   exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_mp2_e_from_psi4_out: SCF energy no&
                   &t found in file "//TRIM(outname)
  close(fid)
  stop
 end if

 if(ri) then ! RI
  i = INDEX(buf, '=')
  read(buf(i+1:),*) ref_e
  do while(.true.)
   read(fid,'(A)') buf
   if(INDEX(buf,'Total Energy') > 0) exit
  end do ! for while
  i = INDEX(buf, '=')
 else        ! noRI
  i = INDEX(buf, ':')
  read(buf(i+1:),*) ref_e
  do while(.true.)
   read(fid,'(A)') buf
   if(INDEX(buf,'=====') > 0) exit
  end do ! for while
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  i = INDEX(buf, ':')
 end if

 read(buf(i+1:),*) mp2_e
 close(fid)
end subroutine read_mp2_e_from_psi4_out

subroutine read_mp2_e_from_qchem_out(outname, ref_e, mp2_e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: ref_e, mp2_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 ref_e = 0d0; mp2_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(2:11)=='Total  MP2' .or. buf(2:13)=='Total  RIMP2' .or. &
     buf(2:11)=='RI-MP2 COR') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "EEROR in subroutine read_mp2_e_from_qchem_out: MP2 energy k&
                   &eywords not found"
  write(6,'(A)') 'in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, '=')
 read(buf(i+1:),*) ref_e ! correlation energy here
 read(fid,'(A)') buf
 close(fid)
 i = INDEX(buf, '=')
 read(buf(i+1:),*) mp2_e
 ref_e = mp2_e - ref_e
end subroutine read_mp2_e_from_qchem_out

subroutine read_mp2_e_from_cfour_out(outname, ref_e, mp2_e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: ref_e, mp2_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 ref_e = 0d0; mp2_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(3:15) == 'Total MBPT(2)') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "EEROR in subroutine read_mp2_e_from_cfour_out: MP2 energy k&
                   &eywords not found"
  write(6,'(A)') 'in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, ':')
 read(buf(i+1:),*) mp2_e

 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(buf(4:17) == 'The total corr') exit
 end do ! for while

 i = INDEX(buf, 'is')
 read(buf(i+2:),*) ref_e
 ref_e = mp2_e - ref_e
end subroutine read_mp2_e_from_cfour_out

subroutine read_mp2_e_from_molcas_out(outname, ref_e, mp2_e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: ref_e, mp2_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 ref_e = 0d0; mp2_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(7:19) == 'Total MBPT2 e') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "EEROR in subroutine read_mp2_e_from_molcas_out: MP2 energy &
                   &keywords not found"
  write(6,'(A)') 'in file '//TRIM(outname)
  close(fid)
  stop
 end if

 read(buf(25:),*) mp2_e

 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(buf(8:12) == 'SCF e') exit
 end do ! for while

 i = INDEX(buf, '=')
 read(buf(i+1:),*) ref_e
end subroutine read_mp2_e_from_molcas_out

subroutine read_cc_e_from_gau_log(outname, t1diag, ref_e, tot_e)
 implicit none
 integer :: i, fid
 real(kind=8) :: ss
 real(kind=8), intent(out) :: t1diag, ref_e, tot_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 t1diag = 0d0; ref_e = 0d0; tot_e = 0d0
 call read_hf_e_and_ss_from_gau_log(outname, ref_e, ss)

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

 i = INDEX(buf, '=')
 read(buf(i+1:),*) tot_e ! CCD/CCSD total energy

 read(fid,'(A)') buf
 if(buf(2:8) == 'T1 Diag') then
  i = INDEX(buf, '=')
  read(buf(i+1:),*) t1diag
 end if

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:14) == 'Discarding MO') exit
  if(buf(2:8) == 'CCSD(T)') then
   i = INDEX(buf, '=')
   read(buf(i+1:),*) tot_e ! CCSD(T) total energy
   exit
  end if
 end do ! for while

 close(fid)
end subroutine read_cc_e_from_gau_log

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

 i = INDEX(buf, '.', back=.true.)
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
  if(INDEX(buf,'HF-SCF') > 0) exit
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
   if(INDEX(buf,'T1 diag') > 0) exit
  end do ! for while
  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_orca_out: 'T1 diag' not&
                   & found in file "//TRIM(outname)
   close(fid)
   stop
  end if
  i = INDEX(buf, 'T1 diagnostic:')
  read(buf(i+14:),*) t1diag
 end if

 close(fid)
end subroutine read_cc_e_from_molpro_out

subroutine read_cc_e_from_pyscf_out(outname, t1diag, ref_e, tot_e)
 implicit none
 integer :: i, k, fid
 real(kind=8) :: ccsd_t_corr
 real(kind=8), intent(out) :: t1diag, ref_e, tot_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 t1diag = 0d0; ref_e = 0d0; tot_e = 0d0; ccsd_t_corr = 0d0
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
   k = INDEX(buf, '=')
   read(buf(k+1:),*) ccsd_t_corr
  case('E(CCSD) =')
   read(buf(10:),*) tot_e
  case('E(RCCSD) ')
   read(buf(11:),*) tot_e
  case('T1_diag =')
   read(buf(10:),*) t1diag
  case('converged')
   exit
  end select
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_cc_e_from_pyscf_out: incomplete fi&
                   &le '//TRIM(outname)
  stop
 end if

 k = INDEX(buf, '=')
 read(buf(k+1:),*) ref_e
 tot_e = tot_e + ccsd_t_corr
end subroutine read_cc_e_from_pyscf_out

! Read CCSD(T) energy from a PSI4 output file. Usually used for the job:
!  grad, wfn = gradient('ccsd(t)', return_wfn=True)
subroutine read_ccsd_t_e_from_psi4_grad_out(outname, RI, t1diag, ref_e, tot_e)
 implicit none
 integer :: i, fid, ntimes
 real(kind=8) :: r1, r2
 real(kind=8), intent(out) :: t1diag, ref_e, tot_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: RI

 t1diag = 0d0; ref_e = 0d0; tot_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 ! read SCF energy
 ntimes = 0
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(5:18) == 'Total Energy =') then
   ntimes = ntimes + 1
   if(ntimes == 2) exit
  end if
 end do ! for while
 read(buf(19:),*) ref_e

 if(RI) then
  ! read T1diag
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(INDEX(buf,'T1 diagnostic') > 0) exit
  end do ! for while
  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_ccsd_t_e_from_psi4_grad_out: no '&
                    &T1 diagnostic' found"
   write(6,'(A)') 'in file '//TRIM(outname)
   close(fid)
   stop
  end if
  read(fid,'(A)') buf
  i = INDEX(buf, ':')
  read(buf(i+1:),*) t1diag

  ! read CCSD(T) energy
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(INDEX(buf,'DF-CCSD(T) Tot') > 0) exit
  end do ! for while
  close(fid)
  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_ccsd_t_e_from_psi4_grad_out: no '&
                    &DF-CCSD(T) Tot'"
   write(6,'(A)') 'found in file '//TRIM(outname)
   stop
  end if
  i = INDEX(buf, ':', back=.true.)

 else ! no RI
  ! read T1diag
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(51:56) == 'T1Diag') exit
  end do ! for while
  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_ccsd_t_e_from_psi4_grad_out: no 'T&
                    &1Diag' found in"
   write(6,'(A)') 'file '//TRIM(outname)
   close(fid)
   stop
  end if
  read(fid,'(A)') ! skip 1 line
  do while(.true.)
   read(fid,'(A)') buf
   if(LEN_TRIM(buf) == 0) exit
   read(buf,*) i, r1, r2, t1diag
  end do ! for while

  ! read CCSD(T) energy
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(9:23) == 'CCSD(T) total e') exit
  end do ! for while
  close(fid)
  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_ccsd_t_e_from_psi4_grad_out: no 'C&
                    &CSD(T) total e'"
   write(6,'(A)') 'found in file '//TRIM(outname)
   stop
  end if
  i = INDEX(buf, '=', back=.true.)
 end if

 read(buf(i+1:),*) tot_e
end subroutine read_ccsd_t_e_from_psi4_grad_out

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
   if(INDEX(buf,'T1 DIAG')>0 .or. INDEX(buf,'T1 diag')>0) exit
  end do ! for while
  i = INDEX(buf, '=>')
  if(i == 0) i = INDEX(buf, ':')
  read(buf(i+2:),*) t1diag
 end if

 rewind(fid)
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(INDEX(buf, key) > 0) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_psi4_out: key '"//key//&
                   "'"
  write(6,'(A)') 'not found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, '=')
 if(i == 0) i = INDEX(buf, ':')
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
  if(INDEX(buf,'SCF E')>0 .or. INDEX(buf,'SCF e')>0) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_psi4_out: 'SCF E' or 'SC&
                   &F e' not found in"
  write(6,'(A)') 'file '//TRIM(outname)
  close(fid)
  stop
 end if

 close(fid)
 i = INDEX(buf, '=')
 if(i == 0) i = INDEX(buf, ':')
 read(buf(i+1:),*) ref_e
end subroutine read_cc_e_from_psi4_out

! read CC energy from GAMESS output
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
  i = INDEX(buf, '=')
  read(buf(i+1:),*) t1diag
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(INDEX(buf, key) > 0) exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_gms_out: key '"//key//"'"
  write(6,'(A)') 'not found in file '//TRIM(outname)
  stop
 end if

 i = INDEX(buf, ':')
 read(buf(i+1:),*) tot_e
end subroutine read_cc_e_from_gms_out

subroutine read_cc_e_from_qchem_out(outname, ccd, ccsd, ccsd_t, t1diag, ref_e, tot_e)
 implicit none
 integer :: i, na, nb, fid
 real(kind=8), intent(out) :: t1diag, ref_e, tot_e
 character(len=10) :: key
 character(len=10), parameter :: key0(3) = ['CCD total ','CCSD total','CCSD(T) to']
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: ccd, ccsd, ccsd_t

 t1diag = 0d0; ref_e = 0d0;  tot_e = 0d0
 if(ccd) then
  key = key0(1)
 else if(ccsd) then
  key = key0(2)
 else if(ccsd_t) then
  key = key0(3)
 end if

 open(newunit=fid,file=TRIM(outname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid, iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid, iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(2:8) == 'SCF ene') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_qchem_out: 'SCF ene' not&
                   & found in file "//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, '=')
 read(buf(i+1:),*) ref_e ! SCF energy

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:11) == key) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_qchem_out: key '"//key//"'"
  write(6,'(A)') 'not found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, '=')
 read(buf(i+1:),*) tot_e

 if(.not. ccd) then ! read T1diag
  read(fid,'(A)') buf
  read(fid,'(A)') buf
  i = INDEX(buf, '=')
  read(buf(i+1:),*) t1diag
  if(t1diag < 5d-3) then
   write(6,'(A)') REPEAT('-',79)
   write(6,'(A)') 'Warning from subroutine read_cc_e_from_qchem_out: T1^2 in Q-&
                  &Chem output file'
   write(6,'(A)') 'has too few digits. This does not affect the electronic ener&
                  &gy at all, but will'
   write(6,'(A)') 'lead to inaccurate T1diag value. If you care about T1diag, y&
                  &ou should use'
   write(6,'(A)') 'another CC_prog.'
   write(6,'(A)') REPEAT('-',79)
  end if
  do while(.true.)
   BACKSPACE(fid)
   BACKSPACE(fid)
   read(fid,'(A)') buf
   if(buf(2:10) == 'Alpha orb') exit
  end do ! for while
  read(fid,'(A)') buf
  read(fid,'(A)') buf
  read(buf(20:),*) i, na
  do while(.true.)
   read(fid,'(A)') buf
   if(buf(2:9) == 'Beta orb') exit
  end do ! for while
  read(fid,'(A)') buf
  read(fid,'(A)') buf
  read(buf(20:),*) i, nb
  t1diag = DSQRT(t1diag/DBLE(2*(na+nb)))
 end if

 close(fid)
end subroutine read_cc_e_from_qchem_out

subroutine read_cc_e_from_cfour_out(outname, ccd, ccsd, ccsd_t, t1diag, ref_e, &
                                    tot_e)
 implicit none
 integer :: i, nocc, fid
 real(kind=8) :: rtmp
 real(kind=8), intent(out) :: t1diag, ref_e, tot_e
 character(len=13) :: key
 character(len=13), parameter :: key0(3) = ['Total CCD ene','Total CCSD en','To&
                                            &tal CCSD(T)']
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: ccd, ccsd, ccsd_t

 ref_e = 0d0; tot_e = 0d0
 if(ccd) then
  key = key0(1)
 else if(ccsd) then
  key = key0(2)
 else if(ccsd_t) then
  key = key0(3)
 end if

 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(1:13) == key) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_cfour_out: key '"//key//"'"
  write(6,'(A)') 'not found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, ':')
 read(buf(i+1:),*) tot_e ! total electronic energy

 do i = 1, 5
  BACKSPACE(fid)
 end do

 if(ccsd_t) then ! read (T) correlation energy
  read(fid,'(A)') buf
  i = INDEX(buf, ':')
  read(buf(i+1:),*) rtmp
  do i = 1, 8
   BACKSPACE(fid)
  end do
 end if

 read(fid,*) i, ref_e ! read CCD/CCSD correlation energy
 close(fid)
 if(ccsd_t) ref_e = ref_e + rtmp
 ref_e = tot_e - ref_e

 if(.not. ccd) then
  call read_nocc_from_cfour_out(outname, nocc)
  call read_t1diag_from_cfour_t1(nocc, t1diag)
 end if
end subroutine read_cc_e_from_cfour_out

subroutine read_posthf_e_from_dalton_out(outname, mp2, ccd, ccsd, ccsd_t, &
                                         t1diag, ref_e, tot_e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: t1diag, ref_e, tot_e
 character(len=14) :: key
 character(len=14), parameter :: key0(5) = ['Total SCF   en','Total MP2   en',&
                           'Total CCD   en','Total CCSD  en','Total energy C']
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: mp2, ccd, ccsd, ccsd_t

 t1diag = 0d0; ref_e = 0d0; tot_e = 0d0
 if(mp2) then
  key = key0(2)
 else if(ccd) then
  key = key0(3)
 else if(ccsd) then
  key = key0(4)
 else if(ccsd_t) then
  key = key0(5)
 end if

 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(13:26) == key) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_dalton_out: key '"//key//"'"
  write(6,'(A)') 'not found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, ':')
 read(buf(i+1:),*) tot_e ! total electronic energy

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(13:26) == key0(1)) exit
 end do ! for while

 i = INDEX(buf, ':')
 read(buf(i+1:),*) ref_e ! SCF electronic energy
 close(fid)

 if(ccsd .or. ccsd_t) then
  write(6,'(/,A)') REPEAT('-',79)
  write(6,'(A)') 'Warning: t1 diagnostic value cannot be printed when CC_prog=D&
                 &alton. Below it'
  write(6,'(A)') 'will be set as zero.'
  write(6,'(A)') REPEAT('-',79)
 end if
end subroutine read_posthf_e_from_dalton_out

subroutine read_cc_e_from_molcas_out(outname, ccsd_t, t1diag, ref_e, tot_e)
 implicit none
 integer :: fid
 real(kind=8), intent(out) :: t1diag, ref_e, tot_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical :: ricd, cholesky
 logical, intent(in) :: ccsd_t

 ricd = .false.; cholesky = .false.   ! initialization

 ! determine RICD/CHOLESKY/noRI
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:4) == 'RICD') then
   ricd = .true.
   exit
  end if
  if(buf(1:8) == 'CHOLESKY') then
   cholesky = .true.
   exit
  end if
  if(buf(1:4) == '()()') exit
 end do ! for while
 close(fid)

 if(ricd .or. cholesky) then
  call read_cc_e_from_molcas_out1(outname, ccsd_t, t1diag, ref_e, tot_e)
 else
  call read_cc_e_from_molcas_out2(outname, ccsd_t, t1diag, ref_e, tot_e)
 end if
end subroutine read_cc_e_from_molcas_out

! read CCSD/CCSD(T) energies from a (Open)Molcas output file (for RICD)
subroutine read_cc_e_from_molcas_out1(outname, ccsd_t, t1diag, ref_e, tot_e)
 implicit none
 integer :: i, j, nocc, fid
 real(kind=8) :: r(2), corr_e
 real(kind=8), intent(out) :: t1diag, ref_e, tot_e
 character(len=240) :: buf, proname
 character(len=240), intent(in) :: outname
 logical, intent(in) :: ccsd_t

 t1diag = 0d0; ref_e = 0d0; tot_e = 0d0
 call read_nocc_from_molcas_out(outname, nocc)
 call find_specified_suffix(outname, '.out', i)
 proname = outname(1:i-1)
 call read_t1diag_from_rstfil(proname, nocc, t1diag)

 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 if(ccsd_t) then
  do while(.true.)
   BACKSPACE(fid,iostat=i)
   if(i /= 0) exit
   BACKSPACE(fid,iostat=i)
   if(i /= 0) exit
   read(fid,'(A)') buf
   if(buf(4:15) == 'CCSD(T) corr') exit
  end do ! for while

  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_molcas_out1: keyword 'CC&
                    &SD(T) corr' not"
   write(6,'(A)') 'found in file '//TRIM(outname)
   close(fid)
   stop
  end if

  i = INDEX(buf, '=')
  read(buf(i+1:),*) corr_e ! CCSD(T) corr e
 end if

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(3:16) == 'Final CCSD ene') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_molcas_out1: keyword 'Fi&
                   &nal CCSD energy' not"
  write(6,'(A)') 'found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 do i = 1, 2   ! read CCSD correlation e
  read(fid,'(A)') buf
  j = INDEX(buf, ':')
  read(buf(j+1:),*) r(i)
 end do ! for i

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(7:17) == 'Total SCF e') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_molcas_out1: keyword 'To&
                   &tal SCF e' not found"
  write(6,'(A)') 'in file '//TRIM(outname)
  stop
 end if

 read(buf(23:),*) ref_e ! SCF e

 if(ccsd_t) then ! CCSD(T)
  tot_e = ref_e + corr_e
 else            ! CCSD
  tot_e = ref_e + r(1) + r(2)
 end if
end subroutine read_cc_e_from_molcas_out1

! read CCSD/CCSD(T) energies from a (Open)Molcas output file (no RI)
subroutine read_cc_e_from_molcas_out2(outname, ccsd_t, t1diag, ref_e, tot_e)
 implicit none
 integer :: i, nocc, fid
 real(kind=8), intent(out) :: t1diag, ref_e, tot_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: ccsd_t

 t1diag = 0d0; ref_e = 0d0; tot_e = 0d0
 call read_nocc_from_molcas_out(outname, nocc)
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(7:18) == 'Total energy') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_molcas_out2: keyword 'To&
                   &tal energy' not"
  write(6,'(A)') 'found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 ! read CCSD energy
 i = INDEX(buf, ':')
 read(buf(i+1:),*) tot_e

 ! read SCF energy
 read(fid,'(A)') buf
 read(fid,'(A)') buf
 i = INDEX(buf, ':')
 read(buf(i+1:),*) ref_e

 ! read t1 diag
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:15) == 'Euclidian norm') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_molcas_out2: keyword 'Eu&
                   &clidian norm' not"
  write(6,'(A)') 'found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, ':')
 read(buf(i+1:),*) t1diag

 t1diag = DSQRT(t1diag*t1diag/DBLE(2*nocc))

 if(ccsd_t) then ! read CCSD(T) energy
  do while(.true.)
   read(fid,'(A)') buf
   if(buf(8:16) == 'CCSD + T3') exit
  end do ! for while
  i = INDEX(buf, '=')
  read(buf(i+1:),*) tot_e
 end if

 close(fid)
end subroutine read_cc_e_from_molcas_out2

subroutine read_eomcc_e_from_pyscf_out(outname, ip_or_ea, nstate, e, mult, fosc)
 implicit none
 integer :: i, j, fid
 integer, intent(in) :: nstate
 real(kind=8), intent(out) :: e(0:nstate), mult(0:nstate), fosc(nstate)
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: ip_or_ea

 buf = ' '; e = 0d0; mult = 0d0; fosc = 0d0
 if(ip_or_ea) mult(1:nstate) = 0.75d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:9) == 'E(CCSD) =') exit
 end do ! for while

 read(buf(10:),*) e(0)

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:15) == 'EOM-CCSD root 0') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_eomcc_e_from_pyscf_out: 'EOM-CCSD &
                   &root 0' not"
  write(6,'(A)') 'located in output file '//TRIM(outname)
  close(fid)
  stop
 end if

 BACKSPACE(fid)
 do i = 1, nstate, 1
  read(fid,'(A)') buf
  j = INDEX(buf, '=')
  read(buf(j+1:),*) e(i)
 end do ! for i

 close(fid)
 forall(i = 1:nstate) e(i) = e(i) + e(0)
end subroutine read_eomcc_e_from_pyscf_out

! read EOM-CCSD energies, spin multiplities and oscillator strengths from Gaussian
! output
subroutine read_eomcc_e_from_gau_log(outname, nstate, e, mult, fosc)
 implicit none
 integer :: i, j, fid
 integer, intent(in) :: nstate ! the number of excited states
 real(kind=8) :: rtmp(4)
 real(kind=8), intent(out) :: e(0:nstate), mult(0:nstate), fosc(nstate)
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 buf = ' '; e = 0d0; mult = 0d0; fosc = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:20) == 'Wavefunction amplit') exit
 end do ! for while

 if(i /= 0) then
  close(fid)
  write(6,'(/,A)') 'ERROR in subroutine read_eomcc_e_from_gau_log: CCSD energy &
                   &is not found'
  write(6,'(A)') 'in file '//TRIM(outname)
  stop
 end if

 i = INDEX(buf, '=', back=.true.)
 read(buf(i+1:),*) e(0)

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:18) == 'Final Eigenvalues') exit
 end do ! for while

 read(fid,'(A)') buf
 do i = 1, nstate, 1
  read(fid,*) j, e(i)
 end do ! for i
 forall(i = 1:nstate) e(i) = e(i) + e(0)

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:18) == 'Ground to excited') exit
 end do ! for while

 read(fid,'(A)') buf
 do i = 1, nstate, 1
  read(fid,*) j, rtmp(1:4), fosc(i)
 end do ! for i

 close(fid)
end subroutine read_eomcc_e_from_gau_log

! read EOM-CCSD energies, spin multiplities and oscillator strengths from ORCA
! output
subroutine read_eomcc_e_from_orca_out(outname, ip_or_ea, nstate, e, mult, fosc)
 implicit none
 integer :: i, j, fid
 integer, intent(in) :: nstate ! the number of excited states
 real(kind=8) :: rtmp(2)
 real(kind=8), intent(out) :: e(0:nstate), mult(0:nstate), fosc(nstate)
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: ip_or_ea

 buf = ' '; e = 0d0; mult = 0d0; fosc = 0d0
 if(ip_or_ea) mult(1:nstate) = 0.75d0

 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:6) == 'E(TOT)') exit
 end do ! for while

 if(i /= 0) then
  close(fid)
  write(6,'(/,A)') "ERROR in subroutine read_eomcc_e_from_orca_out: 'E(TOT)' no&
                   &t found"
  write(6,'(A)') 'in file '//TRIM(outname)
  stop
 end if
 read(buf(47:),*) e(0)

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:12) == 'EOM-CCSD RES') exit
 end do ! for while

 if(i /= 0) then
  close(fid)
  write(6,'(/,A)') "ERROR in subroutine read_eomcc_e_from_orca_out: 'EOM-CCSD R&
                   &ES' not found"
  write(6,'(A)') 'in file '//TRIM(outname)
  stop
 end if

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:6) == 'IROOT=') then
   i = INDEX(buf, ':')
   read(buf(7:i-1),*) j
   read(buf(i+1:),*) e(j)
   if(j == nstate) exit
  end if
 end do ! for while
 forall(i = 1:nstate) e(i) = e(i) + e(0)

 if(.not. ip_or_ea) then
  do while(.true.)
   read(fid,'(A)') buf
   if(buf(10:24) == 'ABSORPTION SPEC') exit
  end do ! for while
  do i = 1, 4
   read(fid,'(A)') buf
  end do ! for i
  do i = 1, nstate, 1
   read(fid,*) j, rtmp(1:2), fosc(i)
  end do ! for i
 end if

 close(fid)
end subroutine read_eomcc_e_from_orca_out

! read EOM-CCSD energies, spin multiplities and oscillator strengths from GAMESS
! output
subroutine read_eomcc_e_from_gms_out(outname, ip, ea, nstate, e, mult, fosc)
 implicit none
 integer :: i, fid
 integer, intent(in) :: nstate ! the number of excited states
 real(kind=8) :: rtmp, r(3)
 real(kind=8), intent(out) :: e(0:nstate), mult(0:nstate), fosc(nstate)
 character(len=1) :: str1
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: ip, ea

 buf = ' '; e = 0d0; mult = 0d0; fosc = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(ip) then
   if(buf(15:31) == 'SUMMARY OF IP-EOM') exit
  else if(ea) then
   if(buf(15:31) == 'SUMMARY OF EA-EOM') exit
  else ! (EE-)EOM-CCSD
   if(buf(21:34) == 'SUMMARY OF EOM') exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_eomcc_e_from_gms_out: EOM-related &
                   &keywords not'
  write(6,'(A)') 'located in file '//TRIM(outname) 
  close(fid)
  stop
 end if

 read(fid,'(A)') buf
 read(fid,'(A)') buf
 if(.not. ip) read(fid,'(A)') buf

 if(ip .or. ea) then
  do i = 1, nstate, 1
   read(fid,*) str1, mult(i), rtmp, e(i)
  end do ! for i
  forall(i = 1:nstate)
   mult(i) = (mult(i)-1d0)*0.5d0
   mult(i) = mult(i)*(mult(i)+1d0)
  end forall
 else   ! (EE-)EOM-CCSD
  do i = 1, nstate, 1
   read(fid,*) str1, rtmp, rtmp, e(i)
  end do ! for i
 end if

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(11:21)=='CCSD ENERGY' .or. buf(9:22)=='CCSD    ENERGY') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_eomcc_e_from_gms_out: ground state&
                  & CCSD energy not'
  write(6,'(A)') 'found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, ':')
 read(buf(i+1:),*) e(0)

 ! if EOM-CCSD, read oscillator strengths
 if(.not. (ip .or. ea)) then
  rewind(fid)

  do i = 1, nstate, 1
   do while(.true.)
    read(fid,'(A)') buf
    if(buf(21:43) == 'EXCITED STATE EOMCCSD P') exit
   end do ! for while

   do while(.true.)
    read(fid,'(A)') buf
    if(buf(2:22) == 'OSCILLATOR   STRENGTH') exit
   end do ! for while

   read(buf(24:),*) r(1:3)
   fosc(i) = SUM(r)
  end do ! for i
 end if

 close(fid)
end subroutine read_eomcc_e_from_gms_out

subroutine read_eomcc_e_from_psi4_out(outname, nstate, e, mult, fosc)
 implicit none
 integer :: i, j, fid
 integer, intent(in) :: nstate
 real(kind=8) :: rtmp(4)
 real(kind=8), intent(out) :: e(0:nstate), mult(0:nstate), fosc(nstate)
 character(len=1) :: str1
 character(len=12) :: str12
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 buf = ' '; e = 0d0; mult = 0d0; fosc = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 ! read CCSD energy and store it into e(0)
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(9:20) == 'CCSD total e') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_eomcc_e_from_psi4_out: 'Ground Sta&
                   &te -> Excited State T'"
  write(6,'(A)') 'not located in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, '=')
 read(buf(i+1:),*) e(0)

 ! read EOM-CCSD electronic energies
 do i = 1, nstate, 1
  write(str12,'(A10,I0)') 'EOM State ', i
  j = LEN_TRIM(str12)

  do while(.true.)
   read(fid,'(A)') buf
   if(buf(1:j) == str12(1:j)) exit
  end do ! for while

  read(buf(13:),*) rtmp(1:3), e(i)
 end do ! for i

 ! read CCSD -> EOM-CCSD oscillator strengths
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(INDEX(buf,'Ground State -> Excited State T') > 0) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_eomcc_e_from_psi4_out: 'Ground Sta&
                   &te -> Excited State T'"
  write(6,'(A)') 'not located in file '//TRIM(outname)
  close(fid)
  stop
 end if

 do i = 1, 3
  read(fid,'(A)') buf
 end do ! for i

 do i = 1, nstate, 1
  read(fid,*) j, str1, rtmp(1:4), fosc(i)
 end do ! for i

 close(fid)
end subroutine read_eomcc_e_from_psi4_out

! dump GAMESS MP2/CC NOs from .dat to .fch, and make sure that Total SCF Density
! in .fch is corresponding (un)relaxed density matrix
subroutine dump_gms_no_dat2fch(datname, no_fch)
 implicit none
 integer :: i, SYSTEM
 character(len=240) :: gmsname
 character(len=240), intent(in) :: datname, no_fch
 character(len=500) :: buf

 call del_vec_in_dat(datname)

 buf = 'dat2fch '//TRIM(datname)//' '//TRIM(no_fch)
 write(6,'(A)') '$'//TRIM(buf)
 i = SYSTEM(TRIM(buf))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine dump_gms_no_dat2fch: failed to call uti&
                   &lity dat2fch.'
  write(6,'(A)') 'datname='//TRIM(datname)
  write(6,'(A)') 'no_fch='//TRIM(no_fch)
  stop
 end if

 call find_specified_suffix(datname, '.dat', i)
 gmsname = datname(1:i-1)//'.gms'

 buf = 'extract_noon2fch '//TRIM(gmsname)//' '//TRIM(no_fch)//' -mp2'
 write(6,'(A)') '$'//TRIM(buf)
 i = SYSTEM(TRIM(buf))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine dump_gms_no_dat2fch: failed to call uti&
                   &lity extract_noon2fch.'
  write(6,'(A)') 'Related files: '//TRIM(gmsname)//', '//TRIM(no_fch)
  stop
 end if

 call update_density_using_no_and_on(no_fch)
end subroutine dump_gms_no_dat2fch

! dump ORCA MP2/CC NOs from .gbw to .fch, and make sure that Total SCF Density
! in .fch is corresponding (un)relaxed density matrix
subroutine dump_orca_no_gbw2fch(no_gbw, fchname)
 use mkl_content, only: check_uhf_in_mkl
 use util_wrapper, only: mkl2fch_wrap, fch_u2r_wrap
 use fch_content, only: check_uhf_in_fch
 implicit none
 integer :: i, SYSTEM, RENAME
 character(len=240) :: mklname, tmp_out, no_fch, nso_fch, nso_no
 character(len=240), intent(in) :: no_gbw, fchname
 logical :: uhf, uhf1

 i = INDEX(no_gbw, '.mdci.nat', back=.true.)
 if(i == 0) i = INDEX(no_gbw, '.', back=.true.)
 mklname = no_gbw(1:i-1)//'.mkl'
 tmp_out = no_gbw(1:i-1)//'.t'
 no_fch = no_gbw(1:i-1)//'_NO.fch'
 nso_fch = no_gbw(1:i-1)//'_NSO.fch'
 nso_no = no_gbw(1:i-1)//'_NSO_NO.fch'

 i = SYSTEM('orca_2mkl '//TRIM(no_gbw)//' '//TRIM(mklname)//' -mkl -anyorbs >'&
          //TRIM(tmp_out)//" 2>&1")
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine dump_orca_no_into_fch: failed to call u&
                   &tility orca_2mkl.'
  stop
 end if

 call check_uhf_in_mkl(mklname, uhf)
 if(uhf) then
  call copy_file(fchname, nso_fch, .false.)
  call mkl2fch_wrap(mklname, nso_fch, 2) ! '-nso'
  call gen_no_from_nso(nso_fch)
  i = RENAME(TRIM(nso_no), TRIM(no_fch))
 else
  call copy_file(fchname, no_fch, .false.)
  call check_uhf_in_fch(no_fch, uhf1)
  if(uhf1) call fch_u2r_wrap(no_fch, no_fch)
  call mkl2fch_wrap(mklname, no_fch, 1) ! '-no'
  call update_density_using_no_and_on(no_fch)
 end if

 call delete_files(3, [no_gbw, mklname, tmp_out])
end subroutine dump_orca_no_gbw2fch

! copy specified density from no_fch to a copy of hf_fch, and generate NOs in
! the copied file
subroutine copy_dm_and_gen_no(no_fch, hf_fch, itype)
 use fch_content, only: check_uhf_in_fch
 use util_wrapper, only: fch_u2r_wrap
 implicit none
 integer :: i, RENAME
 integer, intent(in) :: itype
 ! itype in subroutine read_density_from_fch, e.g.
 ! 5/6/7/8: 'Total MP2 D'/'Spin MP2 De'/'Total CC De'/'Spin CC Den'
 character(len=240) :: no_fch1, nso_fch
 character(len=240), intent(in) :: no_fch, hf_fch
 ! files no_fch and hf_fch both exist before calling this subroutine
 logical :: uhf

 call check_uhf_in_fch(hf_fch, uhf)
 i = LEN_TRIM(no_fch)

 if(uhf) then
  if(no_fch(i-6:i) /= '_NO.fch') then
   write(6,'(/,A)') "ERROR in subroutine copy_dm_and_gen_no: '_NO.fch' suffix n&
                    &ot found in filename "
   write(6,'(A)') TRIM(no_fch)
   stop
  end if
  nso_fch = no_fch(1:i-7)//'_NSO.fch'
  call copy_file(hf_fch, nso_fch, .false.)
 end if

 ! generate spatial NOs using total density
 no_fch1 = no_fch(1:i-4)//'1.fch'
 call copy_file(hf_fch, no_fch1, .false.)
 if(uhf) call fch_u2r_wrap(no_fch1, no_fch1)
 call copy_dm_between_fch(no_fch, no_fch1, itype, .true.)
 call gen_no_using_density_in_fch(no_fch1, 1)

 ! generate natural spin orbitals (NSOs) using alpha/beta density, respectively
 if(uhf) then
  call copy_dm_between_fch(no_fch, nso_fch, itype, .true.)
  call copy_dm_between_fch(no_fch, nso_fch, itype+1, .false.)
  call gen_no_using_density_in_fch(nso_fch, 0)
 end if

 i = RENAME(TRIM(no_fch1), TRIM(no_fch))
end subroutine copy_dm_and_gen_no

! read the number of doubly occupied orbitals (frozen core orbitals not
!  included) from CFOUR output file
subroutine read_nocc_from_cfour_out(outname, nocc)
 implicit none
 integer :: i, na, nf, fid
 integer, intent(out) :: nocc
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical :: frozen

 nocc = 0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(11:26) == 'Alpha population') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_nocc_from_cfour_out: keyword 'Alph&
                   &a population' not"
  write(6,'(A)') 'found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, ':')
 read(buf(i+1:),*) na ! number of alpha spin orbitals

 frozen = .true.
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(3:15) == 'The following') exit
  if(buf(3:12) == 'Summary of') then
   frozen = .false.
   exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_nocc_from_cfour_out: keyword 'The &
                   &following' not"
  write(6,'(A)') 'found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 if(frozen) then
  read(buf(16:),*) nf ! number of frozen core orbitals
 else
  nf = 0
 end if
 close(fid)

 nocc = na - nf
end subroutine read_nocc_from_cfour_out

! read T1 diagnostic value from a CFOUR T1 binary file
subroutine read_t1diag_from_cfour_t1(nocc, t1diag)
 implicit none
 integer :: i, nvir, fid
 integer(kind=4) :: ibuf(13)
 integer, intent(in) :: nocc
! nocc: the number of doubly occupied orbitals (frozen core not included)
 real(kind=8), intent(out) :: t1diag
 real(kind=8), allocatable :: t1(:)
 logical :: alive

 t1diag = 0d0
 inquire(file='T1',exist=alive)
 if(.not. alive) then
  write(6,'(/,A)') 'ERROR in subroutine read_t1diag_from_cfour_t1: file T1 does&
                  & not exist.'
  stop
 end if

 ! determine nvir from the file size of T1
 call STAT('T1', ibuf)
 nvir = ibuf(8)/(8*nocc)

 allocate(t1(nocc*nvir), source=0d0)
 open(newunit=fid,file='T1',status='old',access='stream')
 read(fid,iostat=i) t1
 close(fid)

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_t1diag_from_cfour_t1: failed to re&
                   &ad t1 amplitudes'
  write(6,'(A)') 'from CFOUR T1 file.'
  write(6,'(A,2I6)') 'nocc, nvir=', nocc, nvir
  deallocate(t1)
  stop
 else
  t1diag = DSQRT(DOT_PRODUCT(t1,t1)/DBLE(nocc*2))
  deallocate(t1)
 end if
end subroutine read_t1diag_from_cfour_t1

! read T1 diagnostic value from an (Open)Molcas .RstFil unformatted file
subroutine read_t1diag_from_rstfil(proname, nocc, t1diag)
 implicit none
 integer :: i, nvir, fid
 integer(kind=4) :: ibuf(13)
 integer, intent(in) :: nocc
! nocc: the number of doubly occupied orbitals (frozen core not included)
 real(kind=8), intent(out) :: t1diag
 real(kind=8), allocatable :: t1(:)
 character(len=840) :: fullpath
 character(len=240), intent(in) :: proname

 t1diag = 0d0
 call getenv('MOLCAS_WORKDIR', fullpath)
 fullpath = TRIM(fullpath)//'/'//TRIM(proname)//'/'//TRIM(proname)//'.RstFil'

 ! Determine nvir from the file size of .RstFil. Here 40 means E1 CCSD energy,
 !  E2 CCSD energy and the number of CCSD iterations, so it is excluded.
 call STAT(TRIM(fullpath), ibuf)
 nvir = (ibuf(8)-40)/(8*nocc)

 allocate(t1(nocc*nvir), source=0d0)
 open(newunit=fid,file=TRIM(fullpath),status='old',form='unformatted')
 read(fid,iostat=i) t1
 close(fid)

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_t1diag_from_rstfil: failed to read&
                   & t1 amplitudes'
  write(6,'(A)') 'from Molcas .RstFil file.'
  write(6,'(A)') 'File prefix='//TRIM(proname)
  write(6,'(A,2I6)') 'nocc, nvir=', nocc, nvir
  deallocate(t1)
  stop
 else
  t1diag = DSQRT(DOT_PRODUCT(t1,t1)/DBLE(nocc*2))
  deallocate(t1)
 end if
end subroutine read_t1diag_from_rstfil

! find the number of occupied orbtials from a (Open)Molcas output file (frozen
! core orbitals not included)
subroutine read_nocc_from_molcas_out(outname, nocc)
 implicit none
 integer :: i, fid
 integer, intent(out) :: nocc
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 nocc = 0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(7:37)=='No. of occupied orbitals with a' .or. buf(2:9)=='Occupied') exit
 end do ! for while

 close(fid)

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_nocc_from_molcas_out: nocc cannot &
                   &be found in file"
  write(6,'(A)') TRIM(outname)
  stop
 end if

 i = INDEX(buf, ':')
 if(i > 0) then
  read(buf(i+1:),*) nocc
 else
  read(buf(47:),*) nocc
 end if
end subroutine read_nocc_from_molcas_out

