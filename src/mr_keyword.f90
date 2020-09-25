! written by jxzou at 20200510
! updated by jxzou at 20200805: add background point charges related subroutines
! updated by jxzou at 20200807: add formchk, unfchk, orca_2mkl wrappers

! file unit of printing
module print_id
 implicit none
 integer, parameter :: iout = 6
end module print_id

! molecular information
module mol
 implicit none
 integer :: charge = 0   ! charge
 integer :: mult   = 1   ! spin multiplicity
 integer :: nbf    = 0   ! number of AO basis functions
 integer :: nif    = 0   ! number of independent functions, i.e., nmo
 integer :: ndb    = 0   ! number of doubly occupied orbitals
 integer :: nopen  = 0   ! number of singly occupied orbitals
 integer :: npair  = 0   ! number of pairs in GVB
 integer :: npair0 = 0   ! number of active pairs in GVB (|C2| > 0.1)
 integer :: nacto  = 0   ! number of active orbitals in active space
 integer :: nacte  = 0   ! number of active electrons in active space
 integer :: nacta  = 0   ! number of alpha active electrons in active space
 integer :: nactb  = 0   ! number of beta active electrons in active space
 ! nacte = nacta + nactb
 ! nacta = npair0 + nopen
 ! nactb = npair0
 integer :: natom = 0    ! number of atoms
 integer :: nbgchg = 0   ! number of background point charges
 integer, allocatable :: nuc(:) ! nuclear charge number

 logical :: lin_dep = .false. ! whether basis set linear dependence exists
 ! (1) nbf = nif, lin_dep = .False.;
 ! (2) nbf > nif, lin_dep = .True. ;
 ! (3) nbf < nif is impossible.

 real(kind=8) :: rhf_e    = 0.0d0 ! RHF (electronic) energy
 real(kind=8) :: uhf_e    = 0.0d0 ! UHF energy
 real(kind=8) :: gvb_e    = 0.0d0 ! GVB energy
 real(kind=8) :: casci_e  = 0.0d0 ! CASCI/DMRG-CASCI energy
 real(kind=8) :: casscf_e = 0.0d0 ! CASSCF/DMRG-CASSCF energy
 real(kind=8) :: caspt2_e = 0.0d0 ! CASPT2/DMRG-CASPT2 energy
 real(kind=8) :: nevpt2_e = 0.0d0 ! CASSCF-NEVPT2/DMRG-NEVPT2 energy
 real(kind=8) :: davidson_e=0.0d0 ! Davidson correction energy
 real(kind=8) :: mrcisd_e = 0.0d0 ! MRCISD+Q energy
 real(kind=8) :: ptchg_e  = 0.0d0 ! Coulomb energy of background point charges
 real(kind=8) :: nuc_pt_e = 0.0d0 ! nuclear-point_charge interaction energy
 real(kind=8), allocatable :: coor(:,:)     ! Cartesian coordinates of this molecule
 real(kind=8), allocatable :: grad(:)       ! Cartesian gradient of this molecule, 3*natom
 real(kind=8), allocatable :: bgcharge(:,:) ! background point charges
 character(len=2), allocatable :: elem(:)   ! element symbols

end module mol

! keywords information (default values are set)
module mr_keyword
 use print_id, only: iout
 implicit none
 integer :: mem = 4                 ! memory, default 4 GB
 integer :: nproc = 4               ! number of processors, default 4
 integer :: npair_wish = 0          ! number of GVB pairs specified by user
 integer :: nacto_wish = 0          ! number of active orbitals specified by user
 integer :: nacte_wish = 0          ! number of active electrons specified by user

 integer :: ist = 0              ! the i-th strategy
 ! 0: if RHF wfn is stable, use strategy 3; otherwise use strategy 1
 ! 1: UHF -> UNO -> associated rotation -> GVB -> CASCI/CASSCF -> ...
 ! 2: UHF -> UNO -> (associated rotation ->) CASCI/CASSCF -> ...
 ! 3: RHF -> virtual orbital projection -> localization -> pairing -> GVB -> CASCI/CASSCF -> ...
 ! 4: RHF -> virtual orbital projection -> CASCI/CASSCF -> ...
 ! 5: NOs -> CASCI/CASSCF -> ...

 integer :: CtrType = 0             ! 1/2/3 for Uncontracted-/ic-/FIC- MRCI

 logical :: mo_rhf  = .false.       ! whether the initial wfn is RHF/UHF for True/False
 ! mo_rhf will be set as .True. in the follwing 3 cases:
 ! (1) the computed RHF wfn is stable; (2) readrhf = .True.; (3) readno = .True.
 ! the rhf/uhf variable/flag will be used in utilities like fch2inp

 logical :: cart     = .true.      ! Cartesian/spherical harmonic functions
 logical :: bgchg    = .false.     ! wthether there is back ground charge(s)
 logical :: tencycle = .true.      ! whether to perform 10 cycles of HF after transferring MOs
 logical :: readrhf  = .false.     ! read RHF MOs from a given .fch(k)
 logical :: readuhf  = .false.     ! read UHF MOs from a given .fch(k)
 logical :: readno   = .false.     ! read NOs from a given .fch(k), useful for MP2 and CCSD NOs
 logical :: skiphf   = .false.     ! (readrhf .or. readuhf .or. readno)
 logical :: hardwfn  = .false.     ! whether difficult wavefunction cases
 logical :: crazywfn = .false.     ! whether crazywfn wavefunction cases (e.g. Cr2 at 5 Anstrom)
 ! if hardwfn is .True., AutoMR will add additional keywords to ensure convergence or correct spin
 ! Note: SCF/CASSCF/CASCI will sometimes (rare cases for CASCI) stuck in a saddle point/local minimum
 ! if crazywfn is .True., AutoMR will add more keywords (than hardwfn) to ensure convergence or correct spin

 character(len=4)   :: localm = 'boys' ! localization method: boys/pm
 character(len=240) :: gjfname = ' '   ! filename of the input .gjf file
 character(len=240) :: chgname = ' '   ! filename of the .chg file (background point charges)
 character(len=240) :: hf_fch = ' '    ! filename of the given .fch(k) file
 character(len=240) :: datname = ' '   ! filename of GAMESS GVB .dat file
 character(len=240) :: casnofch = ' '  ! .fch(k) file of CASCI or CASSCF job

 integer :: maxM = 1000             ! bond-dimension in DMRG computation
 logical :: vir_proj = .false.      ! virtual orbitals projection onto those of STO-6G
 logical :: uno = .false.           ! generate UNOs

 logical :: gvb     = .false.
 logical :: casci   = .false.
 logical :: casscf  = .false.
 logical :: dmrgci  = .false.
 logical :: dmrgscf = .false.
 logical :: caspt2  = .false.
 logical :: nevpt2  = .false.
 logical :: mrcisd  = .false.
 logical :: CIonly  = .false.       ! whether to optimize orbitals before caspt2/nevpt2/mrcisd
 logical :: casscf_force = .false.  ! whether to calculate CASSCF force

 character(len=10) :: gvb_prog     = 'gamess'
 character(len=10) :: casci_prog   = 'pyscf'
 character(len=10) :: casscf_prog  = 'pyscf'
 character(len=10) :: dmrgci_prog  = 'pyscf'
 character(len=10) :: dmrgscf_prog = 'pyscf'
 character(len=10) :: caspt2_prog  = 'openmolcas'
 character(len=10) :: nevpt2_prog  = 'pyscf'
 character(len=10) :: mrcisd_prog  = 'openmolcas'
 character(len=10) :: ic_mrcc_prog = ' '

 character(len=240) :: gau_path = ' '
 character(len=240) :: gms_path = ' '
 character(len=240) :: gms_scr_path = ' '
 character(len=240) :: molcas_path = ' '
 character(len=240) :: orca_path = ' '

 character(len=7) :: method = ' '    ! model chemistry, theoretical method
 character(len=13) :: basis = ' '     ! basis set (gen and genecp supported)

contains

 ! repalce variables like '$USER' in path into real path
 subroutine replace_env_in_path(path)
  implicit none
  integer :: i, j, k
  character(len=100) :: str
  character(len=240) :: buf
  character(len=240), intent(inout) :: path

  i = index(path, '$')
  if(i == 0) return

  j = index(path(i+1:),'/')
  if(j == 0) then
   j = LEN_TRIM(path) + 1
  else
   j = j + i
  end if

  str = ' '
  call getenv(path(i+1:j-1), str)
  buf(1:i-1) = path(1:i-1)

  str = ADJUSTL(str)
  k = LEN_TRIM(str)
  buf(i:i+k-1) = TRIM(str)

  if(j > 0) buf(i+k:) = path(j:)

  path = TRIM(buf)
  return
 end subroutine replace_env_in_path

 ! read paths of programs in $MOKIT_ROOT/program.info
 subroutine read_program_path()
  implicit none
  integer :: i, j, fid
  integer(kind=4) :: hostnm
  character(len=8) :: hostname
  character(len=24) :: data_string
  character(len=240) :: buf, path
  logical :: alive

  write(iout,'(A)') 'Output information of AutoMR of MOKIT(Molecular Orbital Kit):'
  hostname = ' '
  data_string = ' '
  j = hostnm(hostname)
  call fdate(data_string)
  write(iout,'(A)') 'HOST '//TRIM(hostname)//', '//TRIM(data_string)

  path = ' '
  call getenv('MOKIT_ROOT', path)
  path = TRIM(path)//'/program.info'
  open(newunit=fid,file=TRIM(path),status='old',position='rewind')

  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:1)=='#' .or. LEN_TRIM(buf)==0) cycle

   i = index(buf,'='); j = index(buf,'#')
   if(j == 0) j = LEN_TRIM(buf) + 1

   if(buf(1:3) == 'GAU') then
    gau_path = buf(i+1:j-1)
   else if(buf(1:3) == 'GMS') then
    gms_path = buf(i+1:j-1)
   else if(buf(1:6) == 'MOLCAS') then
    molcas_path = buf(i+1:j-1)
   else if(buf(1:4) == 'ORCA') then
    orca_path = buf(i+1:j-1)
   end if
  end do

  close(fid)

  inquire(file=TRIM(gms_path),exist=alive)
  if(.not. alive) then
   write(iout,'(A)') 'ERROR in subroutine read_scr_from_rungms: rungms does not exist.'
   write(iout,'(A)') 'gms_path='//TRIM(gms_path)
   stop
  end if

  open(newunit=fid,file=TRIM(gms_path),status='old',position='rewind')
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:7) == 'set SCR') exit
  end do ! for while
  close(fid)

  i = index(buf,'=')
  gms_scr_path = buf(i+1:)

  call replace_env_in_path(gau_path)
  call replace_env_in_path(gms_path)
  call replace_env_in_path(gms_scr_path)
  call replace_env_in_path(molcas_path)
  call replace_env_in_path(orca_path)
  write(iout,'(/,A)') 'Read program paths from $MOKIT_ROOT/program.info:'
  write(iout,'(A)') 'gau_path     = '//TRIM(gau_path)
  write(iout,'(A)') 'gms_path     = '//TRIM(gms_path)
  write(iout,'(A)') 'gms_scr_path = '//TRIM(gms_scr_path)
  write(iout,'(A)') 'molcas_path  = '//TRIM(molcas_path)
  write(iout,'(A)') 'orca_path    = '//TRIM(orca_path)
  return
 end subroutine read_program_path

 subroutine parse_keyword()
  implicit none
  integer :: i, j, k, ifail, nblank, fid
  character(len=24) :: method0 = ' '
  character(len=240) :: buf = ' '
  character(len=1000) :: longbuf = ' '
  logical :: alive(2), alive1(4)

  open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')

  do while(.true.)
   read(fid,'(A)',iostat=ifail) buf
   if(ifail /= 0) exit
   if(buf(1:1) == '#') exit

   call lower(buf)
   i = index(buf,'mem=')
   if(i /= 0) then
    j = index(buf,'gb')
    read(buf(i+4:j-1),*) mem
    cycle
   end if

   i = index(buf,'nproc')
   if(i /= 0) then
    j = index(buf,'=')
    read(buf(j+1:),*) nproc
    cycle
   end if
  end do ! for while

  if(ifail /= 0) then
   write(iout,'(A)') 'ERROR in subroutine parse_keyword: end-of-file detected.'
   write(iout,'(A)') 'The provided .gjf file may be incomplete.'
   stop
  end if

  write(iout,'(/,2(A,I0))') 'memory = ', mem, 'GB     nproc = ', nproc

  call lower(buf)
  i = index(buf,'/')
  if(i == 0) then
   write(iout,'(A)') "ERROR in subroutine parse_keyword: no '/' symbol detected&
                    & in keyword line."
   write(iout,'(A)') 'The method and basis set must be specified via method/basis.'
   stop
  end if

  j = index(buf(1:i-1),' ', back=.true.)
  if(j == 0) then
   write(iout,'(A)') 'ERROR in subroutine parse_keyword: syntax error detected in:'
   write(iout,'(A)') "'"//TRIM(buf)//"'"
   stop
  end if
  method0 = buf(j+1:i-1)

  i = index(method0, '('); j = index(method0, ','); k = index(method0, ')')
  alive = [(i/=0 .and. k/=0), (i==0 .and. k==0)]
  if(.not. (alive(1) .or. alive(2)) ) then
   write(iout,'(A)') 'ERROR in subroutine parse_keyword: incomplete method specified.'
   write(iout,'(A)') 'method = '//TRIM(method0)
   stop
  end if

  if(i /= 0) then
   method = method0(1:i-1)

   if(index(method,'cas')/=0 .or. index(method,'dmrg')/=0) then ! e.g. CAS(6,6) is specified
    read(method0(i+1:j-1),*) nacte_wish
    read(method0(j+1:k-1),*) nacto_wish
    if(nacte_wish<1 .or. nacto_wish<1 .or. nacte_wish/=nacto_wish) then
     write(iout,'(A)') 'ERROR in subroutine parse_keyword: wrong number of active&
                      & electrons/orbitals specified.'
     write(iout,'(2(A,I0))') 'nacte/nacto=', nacte_wish, '/', nacto_wish
     stop
    end if
   else if(index(method,'gvb') /= 0) then ! e.g. GVB(6) is specified
    read(method0(i+1:k-1),*) npair_wish
    if(npair_wish < 1) then
     write(iout,'(A)') 'ERROR in subroutine parse_keyword: wrong number of pairs specified.'
     write(iout,'(A,I0)') 'npair_wish=', npair_wish
     stop
    end if
   else                                   ! not CAS, GVB or DMRG
    write(iout,'(A)') 'ERROR in subroutine parse_keyword: unsupported method '//TRIM(method)
    stop
   end if
  else ! i = 0
   method = TRIM(method0)
  end if

  select case(TRIM(method))
  case('mrcisd','caspt2','nevpt2','casscf','dmrgscf','casci','dmrgci','gvb')
   uno = .true.; gvb = .true.
  case default
   write(iout,'(A)') "ERROR in subroutine parse_keyword: specified method '"//&
                     TRIM(method)//"' not supported."
   stop
  end select

  select case(TRIM(method))
  case('mrcisd')
   mrcisd = .true.
   casscf = .true.
  case('caspt2')
   caspt2 = .true.
   casscf = .true.
  case('nevpt2')
   nevpt2 = .true.
   casscf = .true.
  case('casscf')
   casscf = .true.
  case('dmrgscf')
   dmrgscf = .true.
  case('casci')
   casci = .true.
  case('dmrgci')
   dmrgci = .true.
  case('gvb')
  end select

  i = index(buf,'/'); j = index(buf(i+1:),' ')
  if(j == 0) j = LEN_TRIM(buf)
  basis = buf(i+1:i+j)
  basis = ADJUSTL(basis)
  write(iout,'(A)') 'method/basis = '//TRIM(method)//'/'//TRIM(basis)
  if(npair_wish > 0) write(iout,'(A,I0)') 'User specified GVB npair = ', npair_wish
  if(nacte_wish>0 .and. nacto_wish>0) write(iout,'(2(A,I0))') 'User specified&
                            & CAS nacte/nacto = ',nacte_wish,'/',nacto_wish

  read(fid,'(A)') buf ! skip a black line
  read(fid,'(A)') buf ! the 1st line of keywords
  call lower(buf)

  i = index(buf,'{')
  if(i == 0) then
   write(fid,'(A)') "ERROR in subroutine parse_keyword: 'mokit{}' not detected&
                   & in input file "//TRIM(gjfname)//'.'
   stop
  end if

  j = index(buf,'}')
  if(j == 0) then
   j = LEN_TRIM(buf)
  else
   j = j - 1
  end if
  if(j == i) then ! {}, {
   if(buf(i+1:i+1) == '}') return ! no keyword specified
  else if(j > i) then ! { }
   if(LEN_TRIM(buf(i+1:j)) == 0) return ! no keyword specified
  else if(j < i) then
   write(fid,'(A)') 'ERROR in subroutine parse_keyword: unexpected error j<i.'
   write(fid,'(A)') "It seems that the '}' symbol is before '{'."
   stop
  end if

  k = 1
  if(j > i) then
   longbuf(1:j-i) = buf(i+1:j) ! some keywords specified
   k = j - i + 1 ! the beginning index for next keyword in longbuf
  end if

  if(index(buf,'}') == 0) then ! keywords are written in at least two lines
   do while(.true.)
    read(fid,'(A)',iostat=ifail) buf
    if(ifail /= 0) exit

    call lower(buf)
    i = LEN_TRIM(buf)
    if(index(buf,'}') /= 0) i = i - 1
    longbuf(k:k+i-1) = buf(1:i)
    k = k + i

    if(index(buf,'}') /= 0) exit
   end do ! for while

   if(ifail /= 0) then
    write(iout,'(A)') 'ERROR in subroutine parse_keyword: end-of-file detected.'
    write(iout,'(A)') 'The provided .gjf file may be incomplete.'
    stop
   end if
  end if

  close(fid)
  ! now all keywords are stored in longbuf

  longbuf = ADJUSTL(longbuf)
  write(iout,'(/,A)') 'The keywords in MOKIT{} are merged and shown as follows:'
  write(iout,'(A)') TRIM(longbuf)

  alive1(1:3) = [(index(longbuf,'caspt2_prog')/=0), (index(longbuf,'nevpt2_prog')/=0),&
                 (index(longbuf,'mrcisd_prog')/=0)]
  if(COUNT(alive1(1:3) .eqv. .true.) > 1) then
   write(iout,'(/,A)') "ERROR in subroutine parse_keyword: more than one keyword of&
                      & 'caspt2_prog', 'nevpt2_prog', 'mrcisd_prog' are detected."
   write(iout,'(A)') 'Only one can be specified in a job.'
   stop
  end if

  if(index(longbuf,'casci_prog')/=0 .and. index(longbuf,'casscf_prog')/=0) then
   write(iout,'(/,A)') 'ERROR in subroutine parse_keyword: both casci_prog and&
                    & casscf_prog are detected.'
   write(iout,'(A)') 'Only one can be specified in a job.'
   stop
  else if(index(longbuf,'casci_prog')/=0 .and. index(longbuf,'cionly')==0 .and. casscf) then
   write(iout,'(/,A)') 'ERROR in subroutine parse_keyword: CASSCF activated, but&
                      & you specify the casci_prog.'
   stop
  else if(index(longbuf,'casscf_prog')/=0 .and. casci) then
   write(iout,'(/,A)') 'ERROR in subroutine parse_keyword: CASCI activated, but&
                      & you specify the casscf_prog.'
   stop
  end if

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
   case('readno')
    mo_rhf = .true.
    readno = .true.
    read(longbuf(j+1:i-1),*) hf_fch
   case('nocart')    ! use spherical harmonic functions
    cart = .false.
   case('localm')    ! localization method
    read(longbuf(j+1:i-1),*) localm
    localm = ADJUSTL(localm)
   case('ist')       ! the i-th strategy
    read(longbuf(j+1:i-1),*) ist
   case('ctrtype')   ! unconctracted-/ic-/FIC- MRCI
    read(longbuf(j+1:i-1),*) CtrType
   case('cionly')    ! do CASPT2/NEVPT2 after CASCI, skip CASSCF
    CIonly = .true.
    casscf = .false.; casci = .true.
    dmrgscf = .false.; dmrgci = .true.
   case('no10cycle') ! skip 10 cycle HF
    tencycle = .false.
   case('maxm')
    read(longbuf(j+1:i-1),*) maxM
   case('hardwfn')
    hardwfn = .true.
   case('crazywfn')
    crazywfn = .true.
   case('gvb_prog')
    read(longbuf(j+1:i-1),*) gvb_prog
   case('casci_prog')
    read(longbuf(j+1:i-1),*) casci_prog
   case('casscf_prog')
    read(longbuf(j+1:i-1),*) casscf_prog
   case('dmrgci_prog')
    read(longbuf(j+1:i-1),*) dmrgci_prog
   case('dmrgscf_prog')
    read(longbuf(j+1:i-1),*) dmrgscf_prog
   case('caspt2_prog')
    read(longbuf(j+1:i-1),*) caspt2_prog
   case('nevpt2_prog')
    read(longbuf(j+1:i-1),*) nevpt2_prog
   case('mrcisd_prog')
    read(longbuf(j+1:i-1),*) mrcisd_prog
   case('ic_mrcc_prog')
    read(longbuf(j+1:i-1),*) ic_mrcc_prog
   case('force')
    casscf_force = .true.
   case('charge')
    bgchg = .true.
   case default
    write(iout,'(A)') "ERROR in subroutine parse_keyword: keyword '"//longbuf(1:j-1)&
                      //"' not recognized."
    stop
   end select

   ! delete the keyword which has been specified
   longbuf(1:i) = ' '
   longbuf = ADJUSTL(longbuf)
   if(LEN_TRIM(longbuf) == 0) exit
  end do ! for while

  if(readrhf .or. readuhf .or. readno) then
   inquire(file=TRIM(hf_fch),exist=alive(1))
   if(.not. alive(1)) then
    write(iout,'(A)') 'ERROR in subroutine parse_keyword: file '//TRIM(hf_fch)&
                     //' does not exist.'
    stop
   end if

   skiphf = .true.
   i = index(hf_fch, '.fchk', back=.true.)
   if(i /= 0) then
    buf = hf_fch(1:i-1)//'.fch'
    call copy_file(hf_fch, buf, .false.)
    hf_fch = buf
   end if
  end if

  select case(ist)
  case(0) ! to be determined after RHF and UHF completed
  case(1)
   uno = .true.; gvb = .true.
  case(2)
   uno = .true.; gvb = .false.
  case(3)
   vir_proj = .true.; gvb = .true.
  case(4)
   vir_proj = .true.; gvb = .false.
  case(5)
   gvb = .false.
  case default
   write(iout,'(A)') 'ERROR in subroutine parse_keyword: ist out of range.'
   stop
  end select

  call prt_strategy()
  return
 end subroutine parse_keyword

 subroutine prt_strategy()
  implicit none
  write(iout,'(/,A,I0)') 'No. Strategy = ', ist

  write(iout,'(4(A,L1,3X))') 'readRHF = ', readrhf, 'readUHF = ', readuhf,&
                             'readNO  = ', readno , 'skipHF  = ', skiphf

  write(iout,'(4(A,L1,3X))') 'Cart    = ', cart   , 'Vir_Proj= ', vir_proj,&
                             'UNO     = ', uno    , 'GVB     = ', gvb

  write(iout,'(4(A,L1,3X))') 'CASCI   = ', casci  , 'CASSCF  = ', casscf  ,&
                             'DMRGCI  = ', dmrgci , 'DMRGSCF = ', dmrgscf

  write(iout,'(4(A,L1,3X))') 'CASPT2  = ', caspt2, 'NEVPT2  = ', nevpt2 ,&
                             'MRCISD  = ', mrcisd, 'CIonly  = ', CIonly

  write(iout,'(2(A,L1,3X),A,I1,3X,A)') 'BgCharge= ', bgchg, 'Ana_Grad= ', casscf_force, &
                                       'CtrType = ',CtrType,'LocalM  = '//TRIM(localm)

  if(skiphf) then
   write(iout,'(A)') 'HF_fch = '//TRIM(hf_fch)
  else
   write(iout,'(A)') 'HF_fch = NONE'
  end if

  return
 end subroutine prt_strategy

 subroutine check_kywd_compatible()
  implicit none
  integer :: i
  logical :: alive(3)
  character(len=43), parameter :: error_warn = 'ERROR in subroutine check_kywd_compatible: '

  write(iout,'(/,A)') 'Check if the keywords are compatible with each other...'

  if(hardwfn .and. crazywfn) then
   write(iout,'(A)') error_warn//"only one of 'hardwfn' or 'crazywfn' can be specified."
   stop
  end if

  if(TRIM(localm)/='pm' .and. TRIM(localm)/='boys') then
   write(iout,'(A)') error_warn//"only 'PM' or 'Boys' is supported."
   write(iout,'(A)') 'localm='//TRIM(localm)
   stop
  end if

  alive = [readrhf, readuhf, readno]
  i = COUNT(alive .eqv. .true.)
  if(i > 1) then
   write(iout,'(A)') error_warn//"more than one of 'readrhf', 'readuhf', and 'readno'&
                   & are set .True."
   write(iout,'(A)') 'These three keywords are mutually exclusive.'
   stop
  end if

  if(readrhf .and. .not.(ist==3 .or. ist==4)) then
   write(iout,'(A)') error_warn//"'readrhf' is only compatible with ist=3 or 4."
   stop
  end if

  if(.not.readrhf .and. (ist==3 .or. ist==4)) then
   write(iout,'(A)') error_warn//"ist=3 or 4 specified, it must be used combined with 'readrhf'."
   stop
  end if

  if(readuhf .and. .not.(ist==1 .or. ist==2)) then
   write(iout,'(A)') error_warn//"'readuhf' is only compatible with ist=1 or 2."
   stop
  end if

  if(.not.readuhf .and. (ist==1 .or. ist==2)) then
   write(iout,'(A)') error_warn//"ist=1 or 2 specified, it must be used combined with 'readuhf'."
   stop
  end if

  if(readno .and. ist/=5) then
   write(iout,'(A)') error_warn//"'readno' is only compatible with ist=5."
   stop
  end if

  if(.not.readno .and. ist==5) then
   write(iout,'(A)') error_warn//"ist=5 specified, it must be used combined with 'readno'."
   stop
  end if

  if(gvb .and. (.not.cart)) then
   write(iout,'(A)') error_warn//'GVB computations will be done using GAMESS.'
   write(iout,'(A)') 'This is incompatible with cart=.false.'
   stop
  end if

  if(CIonly .and. (.not.caspt2) .and. (.not.nevpt2) .and. (.not.mrcisd)) then
   write(iout,'(A)') error_warn//"keyword 'CIonly' can only be used in CASPT2/NEVPT2/MRCISD computations."
   write(iout,'(A)') 'But none of CASPT2/NEVPT2/MRCISD is specified.'
   stop
  end if

  if(gvb_prog /= 'gamess') then
   write(iout,'(A)') error_warn//"only 'gamess' is supported for the GVB computation."
   write(iout,'(A)') 'User specified GVB program cannot be identified: '//TRIM(gvb_prog)
   stop
  end if

  if(casci_prog/='gaussian' .and. casci_prog/='gamess' .and. casci_prog/='openmolcas'&
     .and. casci_prog/='orca' .and. casci_prog/='pyscf') then
   write(iout,'(A)') error_warn
   write(iout,'(A)') 'User specified CASCI program cannot be identified: '//TRIM(casci_prog)
   stop
  end if

  if(casscf_prog/='gaussian' .and. casscf_prog/='gamess' .and. casscf_prog/='openmolcas'&
     .and. casscf_prog/='orca' .and. casscf_prog/='pyscf') then
   write(iout,'(A)') error_warn
   write(iout,'(A)') 'User specified CASSCF program cannot be identified: '//TRIM(casscf_prog)
   stop
  end if

  if(mrcisd_prog/='gaussian' .and. mrcisd_prog/='orca' .and. mrcisd_prog/='openmolcas') then
   write(iout,'(A)') error_warn
   write(iout,'(A)') 'User specified MRCISD program cannot be identified:'//TRIM(mrcisd_prog)
   stop
  end if

  if(mrcisd) then
   select case(CtrType)
   case(1) ! uncontracted MRCISD
   case(2) ! ic-MRCISD
    if(mrcisd_prog/='openmolcas') then
     write(iout,'(A)') error_warn
     write(iout,'(A)') 'The ic-MRCISD is only supported by OpenMolcas. But you&
                      & specify mrcisd_prog='//TRIM(mrcisd_prog)
     stop
    end if
   case(3) ! FIC-MRCISD
    if(mrcisd_prog/='orca') then
     write(iout,'(A)') error_warn
     write(iout,'(A)') 'The FIC-MRCISD is only supported by ORCA. But you&
                      & specify mrcisd_prog='//TRIM(mrcisd_prog)
     stop
    end if
   case default
    write(iout,'(A)') error_warn
    write(iout,'(A)') 'You should specify a valid CtrType=1/2/3 for uncontracted/ic-/FIC- MRCISD.'
    stop
   end select

   if(mrcisd_prog=='gaussian' .and. CtrType/=1) then
    write(iout,'(A)') error_warn
    write(iout,'(A,I0)') 'Gaussian can only perform uncontracted MRCISD. But you&
                        & specify CtrType=', CtrType
    stop
   end if

   if(mrcisd_prog=='orca' .and. cart) then
    write(iout,'(A)') error_warn//'conflict settings.'
    write(iout,'(A)') 'ORCA is set as the MRCI_prog, and it only supports spherical&
                     & harmonic functions,'
    write(iout,'(A)') 'but Cartesian functions are used.'
    stop
   end if
  end if

  if(casci_prog=='orca' .or. casscf_prog=='orca') then
   if(gvb .or. cart) then
    write(iout,'(A)') error_warn
    if(gvb) then
     write(iout,'(A)') 'ORCA and GAMESS both activated, but the former uses spherical functions,'
     write(iout,'(A)') 'while the latter uses Cartesian functions. Conflict.'
     write(iout,'(A)') 'You may try ist = 2, where GVB is skipped: UNO -> CASCI/CASSCF'
    end if
    if(cart) then
     write(iout,'(A)') "ORCA is specified. But ORCA uses spherical&
                      & functions, please also specify 'nocart'."
     stop
    end if
   end if
  end if

  if((casci_prog=='gamess' .or. casscf_prog=='gamess') .and. (.not.cart)) then
   write(iout,'(A)') error_warn
   write(iout,'(A)') "GAMESS is specified. But GAMESS uses Cartesian functions, do not specify 'nocart'."
   stop
  end if

  if(caspt2_prog /= 'openmolcas') then
   write(iout,'(A)') error_warn
   write(iout,'(A)') 'User specified CASPT2 program cannot be identified: '//TRIM(caspt2_prog)
   stop
  end if

  if(nevpt2_prog /= 'pyscf') then
   write(iout,'(A)') error_warn
   write(iout,'(A)') 'User specified NEVPT2 program cannot be identified: '//TRIM(nevpt2_prog)
   stop
  end if

  if(ic_mrcc_prog /= ' ') then
   write(iout,'(A)') error_warn
   write(iout,'(A)') 'Currently AutoMR of MOKIT does not support the ic-MRCC method.'
   stop
  end if

  if(casscf_force .and. (.not.casscf)) then
   write(iout,'(A)') error_warn//"'force' keyword is only avaible for CASSCF."
   write(iout,'(A)') 'But CASSCF is not activated.'
   stop
  end if

  if(casscf_force .and. cart .and. casscf_prog=='pyscf') then
   write(iout,'(A)') error_warn//"current version of PySCF can only compute force"
   write(iout,'(A)') 'using spherical harmonic basis fucntions.'
   stop
  end if

  write(iout,'(A)') 'Check done. All keywords are compatible.'
  return
 end subroutine check_kywd_compatible

 ! turn letters in buf into lower case, except those in symbol ''
 subroutine lower(buf)
  implicit none
  integer :: i, j, k1, k2
  character(len=240), intent(inout) :: buf

  k1 = index(buf,"'")
  k2 = index(buf,"'",back=.true.)

  do i = 1, LEN_TRIM(buf), 1
   if(i>k1 .and. i<k2) cycle

   j = IACHAR(buf(i:i))
   if(j>64 .and. j<91) buf(i:i) = ACHAR(j+32)
  end do ! for i

  return
 end subroutine lower

 ! read background point charge(s) from .gjf file
 subroutine read_bgchg_from_gjf(no_coor)
  use print_id, only: iout
  use mol, only: natom, nbgchg, bgcharge, ptchg_e, nuc_pt_e, nuc, coor
  implicit none
  integer :: i, fid, nblank, nblank0
  character(len=240) :: buf
  character(len=41), parameter :: error_warn='ERROR in subroutine read_bgchg_from_gjf: '
  logical, intent(in) :: no_coor

  nblank = 0
  if(no_coor) then ! no Cartesian Coordinates
   nblank0 = 2
  else             ! there exists Cartesian Coordinates
   nblank0 = 3
  end if

  open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(LEN_TRIM(buf) == 0) nblank = nblank + 1
   if(nblank == nblank0) exit
  end do ! for while

  if(i /= 0) then
   write(iout,'(A)') error_warn//'wrong format of background point charges.'
   stop
  end if

  nbgchg = 0

  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i/=0 .or. LEN_TRIM(buf)==0) exit
   nbgchg = nbgchg + 1
  end do ! for while

  if(nbgchg == 0) then
   write(iout,'(A)') error_warn//'no background point charge(s) found in'
   write(iout,'(A)') 'file: '//TRIM(gjfname)
   stop
  end if

  write(iout,'(A,I0)') 'Background point charge specified: nbgchg = ', nbgchg
  allocate(bgcharge(4,nbgchg), source=0.0d0)

  rewind(fid)   ! jump to the 1st line of the file
  nblank = 0
  do while(.true.)
   read(fid,'(A)') buf
   if(LEN_TRIM(buf) == 0) nblank = nblank + 1
   if(nblank == nblank0) exit
  end do ! for while

  do i = 1, nbgchg, 1
   read(fid,*) bgcharge(1:4,i)
  end do ! for i

  close(fid)

  call calc_Coulomb_energy_of_charges(nbgchg, bgcharge, ptchg_e)
  write(iout,'(A,F18.8,A)') 'Coulomb interaction energy of background point&
                           & charges:', ptchg_e, ' a.u.'
  write(iout,'(A)') 'This energy are taken into account for all energies below.'

  i = index(gjfname, '.gjf', back=.true.)
  chgname = gjfname(1:i-1)//'.chg'
  call write_charge_into_chg(nbgchg, bgcharge, chgname)

  call calc_nuc_pt_e(nbgchg, bgcharge, natom, nuc, coor, nuc_pt_e)
  return
 end subroutine read_bgchg_from_gjf

end module mr_keyword

! delete the specified file (if not exist, return)
subroutine delete_file(fname)
 implicit none
 integer :: fid
 character(len=240), intent(in) :: fname
 logical :: alive

 inquire(file=TRIM(fname),exist=alive)

 if(.not. alive) then
  return
 else
  open(newunit=fid,file=TRIM(fname),status='old')
  close(fid,status='delete')
 end if

 return
end subroutine delete_file

! copy file fname1 to fname2 (if delete=.True., delete fname1)
subroutine copy_file(fname1, fname2, delete)
 use print_id, only: iout
 implicit none
 integer :: i, fid1, fid2
 character(len=240) :: buf
 character(len=240), intent(in) :: fname1, fname2
 logical, intent(in) :: delete
 logical :: alive

 inquire(file=TRIM(fname1),exist=alive)
 if(.not. alive) then
  write(iout,'(A)') 'ERROR in subroutine copy_file: fname1 does not exist.'
  write(iout,'(A)') 'fname1: '//TRIM(fname1)
  stop
 end if

 open(newunit=fid1,file=TRIM(fname1),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(fname2),status='replace')

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 if(delete) then
  close(fid1, status='delete')
 else
  close(fid1)
 end if

 close(fid2)
 return
end subroutine copy_file

! calculate the Coulomb interaction energy of point charges
subroutine calc_Coulomb_energy_of_charges(n, charge, e)
 use print_id, only: iout
 implicit none
 integer :: i, j
 integer, intent(in) :: n
 real(kind=8) :: q1, q2, rtmp1(3), rtmp2(3)
 real(kind=8), intent(in) :: charge(4,n)
 ! charge(1:3,i) is the Cartesian coordinates of the i-th point charge
 ! charge(4,i) is the electronic charge of the i-th point charge
 real(kind=8), intent(out) :: e
 real(kind=8), parameter :: zero1 = 1.0d-2, zero2 = 1.0d-3
 real(kind=8), parameter :: Bohr_const = 0.52917721092d0
 real(kind=8), allocatable :: r(:,:)

 e = 0.0d0
 allocate(r(n,n))

 do i = 1, n-1, 1
 rtmp1 = charge(1:3,i)

  do j = i+1, n, 1
   rtmp2 = (rtmp1 - charge(1:3,j))/Bohr_const

   r(j,i) = DSQRT(DOT_PRODUCT(rtmp2, rtmp2))

   if(r(j,i) < zero2) then
    write(iout,'(A)') 'ERROR in subroutine calc_Coulomb_energy_of_charges:'
    write(iout,'(2(A,I0))') 'There exists two point charges too close: i=',i,',j=',j
    write(iout,'(A,F15.8)') 'r(j,i) = ', r(j,i)
    stop
   else if(r(j,i) < zero1) then
    write(iout,'(A)') 'Warning in subroutine calc_Coulomb_energy_of_charges:'
    write(iout,'(2(A,I0))') 'There exists two point charges very close: i=',i,',j=',j
    write(iout,'(A,F15.8)') 'r(j,i) = ', r(j,i)
   end if

  end do ! for j
 end do ! for i

 do i = 1, n-1, 1
  q1 = charge(4,i)

  do j = i+1, n, 1
   q2 = charge(4,j)

   e = e + q1*q2/r(j,i)
  end do ! for j
 end do ! for i

 deallocate(r)
 return
end subroutine calc_Coulomb_energy_of_charges

! write point charges into a .chg file
subroutine write_charge_into_chg(n, charge, chgname)
 implicit none
 integer :: i, fid
 integer, intent(in) :: n
 real(kind=8), intent(in) :: charge(4,n)
 character(len=240), intent(in) :: chgname

 open(newunit=fid,file=TRIM(chgname),status='replace')
 write(fid,'(I0)') n

 do i = 1, n, 1
  write(fid,'(4(1X,F18.8))') charge(1:4,i)
 end do ! for i

 close(fid)
 return
end subroutine write_charge_into_chg

! calculate nuclear-point_charge interaction energy
subroutine calc_nuc_pt_e(nbgchg, bgcharge, natom, nuc, coor, nuc_pt_e)
 implicit none
 integer :: i, j
 integer, intent(in) :: nbgchg, natom
 integer, intent(in) :: nuc(natom)
 real(kind=8) :: rtmp1(3), rtmp2(3), dis, pt_e
 real(kind=8), parameter :: Bohr_const = 0.52917721092d0
 real(kind=8), intent(in) :: bgcharge(4,nbgchg), coor(3,natom)
 real(kind=8), intent(out) :: nuc_pt_e

 nuc_pt_e = 0.0d0

 do i = 1, nbgchg, 1
  pt_e = bgcharge(4,i)
  rtmp1 = bgcharge(1:3,i)

  do j = 1, natom, 1
   rtmp2 = (rtmp1 - coor(:,j))/Bohr_const
   dis = DSQRT(DOT_PRODUCT(rtmp2, rtmp2))

   nuc_pt_e = nuc_pt_e + pt_e*DBLE(nuc(j))/dis
  end do ! for j
 end do ! for i

 return
end subroutine calc_nuc_pt_e

! read nuclear charge number from a given .fch file
subroutine read_nuc_from_fch(natom, nuc, fchname)
 implicit none
 integer :: i, fid
 integer, intent(in) :: natom
 integer, intent(out) :: nuc(natom)
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:14) == 'Atomic numbers') exit
 end do ! for while
 read(fid,'(6(1X,I11))') (nuc(i),i=1,natom)

 close(fid)
 return
end subroutine read_nuc_from_fch

