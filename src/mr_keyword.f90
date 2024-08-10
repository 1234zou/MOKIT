! written by jxzou at 20200510
! updated by jxzou at 20200805: add background point charges related subroutines
! updated by jxzou at 20200807: add formchk, unfchk, orca_2mkl wrappers
! updated by jxzou at 20201211: add Molpro interfaces
! updated by jxzou at 20201218: add 'back = .true.' for detecting GAUSS_EXEDIR
! updated by jxzou at 20210129: add RI options and aux_basis sets auto-determination

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
 integer :: natom = 0      ! number of atoms
 integer :: nbgchg = 0     ! number of background point charges
 integer :: nfrag = 0      ! number of fragments
 integer :: chem_core = 0  ! ECP plus core orbitals of the molecule
 integer :: ecp_core = 0   ! core orbitals replaced by PP/ECP during computing
 integer :: scan_itype = 0 ! the type of scanning, 1/2/3 for bond/angle/dihedral
 integer :: scan_atoms(4)  ! 2/3/4 atoms to be scanned
 integer, allocatable :: nuc(:) ! nuclear charge number

 integer, allocatable :: frag_char_mult(:,:)
 ! charge and spin multiplicity of each fragment, size (2,nfrag)
 integer, allocatable :: atom2frag(:) ! atom to fragment map, size natom

 logical :: lin_dep = .false.    ! whether basis set linear dependence exists
 ! (1) nbf = nif, lin_dep = .False.;
 ! (2) nbf > nif, lin_dep = .True. ;
 ! (3) nbf < nif is impossible.

 logical :: beyond_xe = .false. ! whether there is any element > Xe
 ! Used in ist=3. H~Xe is the range of STO-6G in Gaussian.

 real(kind=8) :: rhf_e     = 0d0 ! RHF (electronic) energy
 real(kind=8) :: uhf_e     = 0d0 ! UHF energy
 real(kind=8) :: gvb_e     = 0d0 ! GVB energy
 real(kind=8) :: casci_e   = 0d0 ! CASCI/DMRG-CASCI energy
 real(kind=8) :: casscf_e  = 0d0 ! CASSCF/DMRG-CASSCF energy
 real(kind=8) :: caspt2_e  = 0d0 ! CASPT2/DMRG-CASPT2 energy
 real(kind=8) :: caspt3_e  = 0d0 ! CASPT3 energy
 real(kind=8) :: nevpt2_e  = 0d0 ! CASSCF-NEVPT2/DMRG-NEVPT2 energy
 real(kind=8) :: nevpt3_e  = 0d0 ! CASSCF-NEVPT3 energy
 real(kind=8) :: mrmp2_e   = 0d0 ! MRMP2 energy
 real(kind=8) :: ovbmp2_e  = 0d0 ! OVB-MP2 energy
 real(kind=8) :: sdspt2_e  = 0d0 ! SDSPT2 energy
 real(kind=8) :: davidson_e= 0d0 ! Davidson correction energy
 real(kind=8) :: mrcisd_e  = 0d0 ! MRCISD+Q energy
 real(kind=8) :: mcpdft_e  = 0d0 ! MC-PDFT energy
 real(kind=8) :: mrcc_e    = 0d0 ! FIC-MRCC/MkMRCC/BWMRCC/BCCC energy
 real(kind=8) :: ptchg_e   = 0d0 ! Coulomb energy of background point charges
 real(kind=8) :: nuc_pt_e  = 0d0 ! nuclear-point_charge interaction energy
 real(kind=8), allocatable :: coor(:,:)     ! Cartesian coordinates of this molecule
 real(kind=8), allocatable :: grad(:)       ! Cartesian gradient of this molecule, 3*natom
 real(kind=8), allocatable :: bgcharge(:,:) ! background point charges

 real(kind=8), allocatable :: sa_cas_e(:) ! multi-root CASCI energies in SA-CASSCF
 real(kind=8), allocatable :: ci_mult(:)  ! spin multiplicities of excited states
 ! size 0:nstate, the ground states is included as 0
 real(kind=8), allocatable :: fosc(:)     ! oscillator strengths, size nstate

 character(len=2), allocatable :: elem(:)   ! element symbols
end module mol

! keywords information (default values are set)
module mr_keyword
 use mol, only: nfrag
 implicit none
 integer :: mem = 4        ! memory, default 4 GB
 integer :: nproc = 4      ! number of processors, default 4
 integer :: npair_wish = 0 ! number of GVB pairs specified by user
 integer :: nacto_wish = 0 ! number of active orbitals specified by user
 integer :: nacte_wish = 0 ! number of active electrons specified by user

 integer :: ist = 0        ! the i-th strategy
 ! 0: if RHF wfn is stable, use strategy 3; otherwise use strategy 1
 ! 1: UHF -> UNO -> associated rotation -> GVB -> CASCI/CASSCF -> ...
 ! 2: UHF -> UNO -> (associated rotation ->) CASCI/CASSCF -> ...
 ! 3: RHF -> virtual orbital projection -> localization -> pairing -> GVB -> CASCI/CASSCF -> ...
 ! 4: RHF -> virtual orbital projection -> CASCI/CASSCF -> ...
 ! 5: NOs -> CASCI/CASSCF -> ...

 integer :: eist = 0       ! the i-th strategy for excited state calculation
 ! 0: if RHF wfn is stable, use strategy 2; otherwise use strategy 1
 ! 1: RHF -> CIS NTO -> CASCI/CASSCF -> ...
 ! 2: UHF -> UCIS NTO -> CASCI/CASSCF -> ...

 integer :: nskip_uno = 0  ! the number of pair of UNO to be skipped when
 ! performing orbital localization
 ! 1 for skipping HONO and LUNO
 ! 2 for skipping HONO-1, HONO, LUNO, and LUNO+1
 ! This keyword is used to try to obtain biradical GVB-SCF solution

 integer :: CtrType = 0    ! 1/2/3 for Uncontracted-/ic-/FIC- MRCI
 integer :: mrcc_type = 1
 ! 1~8 for FIC-MRCC/MkMRCCSD/MkMRCCSD(T)/BWMRCCSD/BWMRCCSD(T)/BCCC2b,3b,4b
 integer :: maxM = 1000    ! bond-dimension in DMRG computation
 integer :: scan_nstep = 0 ! number of steps to scan

 integer :: nstate = 0 ! number of excited states in SA-CASSCF
 ! ground state not included, so nstate=0 means only ground state

 integer :: iroot = 0 ! the state you are interested in SS-CASSCF
 ! 0/1 for ground state/the 1st excited state. Related variable: ss_opt
 integer :: target_root = 0 ! the real i-th number of the interested state
 integer :: xmult = 1 ! spin multiplicity of the target excited state

 real(kind=8) :: uno_thres = 1d-5 ! threshold for UNO occupation number
 ! uno_thres is used for ist = 1,2
 real(kind=8) :: on_thres = 2d-2 ! threshold for NO occupation number, ist=5
 real(kind=8), allocatable :: scan_val(:) ! values of scanned variables

 character(len=4) :: GVB_conv = '5d-4'
 ! density matrix convergence criterion of GVB (1d-5 default in GAMESS)
 character(len=4) :: localm = 'pm'    ! localization method: boys/pm
 character(len=9) :: otpdf = 'tPBE'   ! on-top pair density functional
 character(len=240) :: gjfname = ' '  ! filename of the input .gjf file
 character(len=240) :: chgname = ' '  ! filename of the .chg file (background point charges)
 character(len=240) :: hf_fch = ' '   ! filename of the given .fch(k) file
 character(len=240) :: datname = ' '  ! filename of GAMESS GVB .dat file
 character(len=240) :: casnofch = ' ' ! .fch(k) file of CASCI or CASSCF job
 character(len=240) :: basname = ' '  ! file to store gen/genecp data

 logical :: molcas_omp = .true.  ! OpenMP/MPI version of OpenMolcas
 logical :: dalton_mpi = .false. ! MKL/MPI version of Dalton

 logical :: mo_rhf  = .false.     ! whether the initial wfn is RHF/UHF for True/False
 ! mo_rhf will be set as .True. in the follwing 3 cases:
 ! (1) the computed RHF wfn is stable; (2) readrhf = .True.; (3) readno = .True.
 ! the rhf/uhf variable/flag will be used in utilities like fch2inp

 logical :: cart     = .false.    ! Cartesian/spherical harmonic functions
 logical :: DKH2     = .false.    ! scalar relativistic Douglas-Kroll-Hess 2nd order correction
 logical :: X2C      = .false.    ! scalar relativistic eXact-Two-Component correction
 logical :: dkh2_or_x2c = .false. ! (DKH2 .or. X2C)
 logical :: bgchg    = .false.    ! wthether there is back ground charge(s)
 logical :: readrhf  = .false.    ! read RHF MOs from a given .fch(k)
 logical :: readuhf  = .false.    ! read UHF MOs from a given .fch(k)
 logical :: readno   = .false.    ! read NOs from a given .fch(k), useful for MP2 and CCSD NOs
 logical :: skiphf   = .false.    ! (readrhf .or. readuhf .or. readno)
 logical :: hardwfn  = .false.    ! whether difficult wavefunction cases
 logical :: crazywfn = .false.    ! whether crazywfn wavefunction cases (e.g. Cr2 at 5 Anstrom)
 ! If hardwfn is .True., AutoMR will add additional keywords to ensure convergence
 ! or correct spin. SCF/CASSCF/CASCI will sometimes (rare cases for CASCI) stuck
 ! in a saddle point/local minimum. If crazywfn is .True., AutoMR will add more
 ! keywords (than hardwfn)

 logical :: vir_proj = .false.    ! virtual orbitals projection onto those of STO-6G
 logical :: uno = .false.         ! generate UNOs
 logical :: inherit = .false.     ! whether to inherit keywords in GVB/STO-6G
 logical :: frag_guess = .false.
 ! whether to perform UHF using initial guess constructed from fragments

 logical :: gvb     = .false.
 logical :: casci   = .false.
 logical :: casscf  = .false.
 logical :: dmrgci  = .false.
 logical :: dmrgscf = .false.
 logical :: caspt2  = .false.
 logical :: caspt2k = .false.
 logical :: caspt3  = .false.
 logical :: nevpt2  = .false.
 logical :: nevpt3  = .false.
 logical :: mrmp2   = .false.
 logical :: ovbmp2  = .false.
 logical :: sdspt2  = .false.
 logical :: mrcisd  = .false.
 logical :: mrcisdt = .false. ! uncontracted MRCISDT
 logical :: mcpdft  = .false.
 logical :: mrcc    = .false.
 logical :: fcgvb   = .true.  ! GVB with all doubly occupied orbitals frozen
 logical :: c_fcgvb = .false. ! whether the user has changed defaultFcGVB
 logical :: c_gvb_conv = .false. ! whether the user has changed default GVB_conv
 logical :: HFonly  = .false. ! stop after the HF calculations
 logical :: CIonly  = .false.     ! whether to optimize orbitals before caspt2/nevpt2/mrcisd
 logical :: dyn_corr= .false.     ! dynamic correlation, post-GVB or post-CAS
 logical :: force = .false.       ! whether this is a force calculation
 logical :: casscf_force = .false.! whether to calculate CASSCF force
 logical :: caspt2_force = .false.! whether to calculate CASPT2 force
 logical :: nevpt2_force = .false.! whether to calculate NEVPT2 force
 logical :: mcpdft_force = .false.! whether to calculate MC-PDFT force
 logical :: dmrg_no = .true.      ! whether to generate NO in a DMRG-CASCI job
 logical :: block_mpi = .false.   ! OpenMP or MPI calling the Block program
 logical :: FIC = .false.         ! False/True for FIC-/SC-NEVPT2
 logical :: RI = .false.          ! whether to RI approximation in CASSCF, NEVPT2
 logical :: F12 = .false.         ! whether F12 used in NEVPT2, MRCI
 logical :: DLPNO = .false.       ! whether to turn on DLPNO-NEVPT2
 logical :: pop = .false.         ! whether to perform population analysis
 logical :: nmr = .false.         ! whether to calcuate nuclear shielding
 logical :: ICSS = .false.        ! whether to calcuate ICSS
 logical :: soc = .false.         ! whether to calcuate spin-orbit coupling (SOC)
 logical :: excludeXH = .false.   ! whether to exclude inactive X-H bonds from GVB
 logical :: onlyXH = .false.      ! whether to keep only X-H bonds in GVB
 logical :: LocDocc = .false.     ! whether to localize GVB doubly occupied orb
 logical :: rigid_scan = .false.  ! rigid/unrelaxed PES scan
 logical :: relaxed_scan = .false.! relaxed PES scan
 logical :: ss_opt = .false.      ! State-specific orbital optimization
 logical :: excited = .false.     ! whether to perform excited states calculations
 logical :: mixed_spin = .false.  ! allow multiple spin in SA-CASSCF, e.g. S0/T1
 logical :: TDHF = .false.        ! True/False for CIS/TDHF
 logical :: sa_cas = .false.      ! State-Averaged CASSCF
 logical :: QD = .false.
 ! False: multi-root CASCI-based NEVPT2 or CASPT2
 ! True : QD-NEVPT2 or (X)MS-CASPT2

 character(len=10) :: hf_prog      = 'gaussian'
 character(len=10) :: gvb_prog     = 'gamess'
 character(len=10) :: casci_prog   = 'pyscf'
 character(len=10) :: casscf_prog  = 'pyscf'
 character(len=10) :: dmrgci_prog  = 'pyscf'
 character(len=10) :: dmrgscf_prog = 'pyscf'
 character(len=10) :: caspt2_prog  = 'openmolcas'
 character(len=10) :: nevpt2_prog  = 'pyscf'
 character(len=10) :: mrmp2_prog   = 'gamess'
 character(len=10) :: mrcisd_prog  = 'openmolcas'
 character(len=10) :: mrcisdt_prog = 'openmolcas' ! uncontracted MRCISDT
 character(len=10) :: mcpdft_prog  = 'openmolcas'
 character(len=10) :: mrcc_prog    = 'orca'
 character(len=10) :: cis_prog     = 'gaussian'

 character(len=240) :: mokit_root = ' '
 character(len=240) :: gau_path = ' '
 character(len=240) :: gms_path = ' '
 character(len=240) :: gms_scr_path = ' '
 character(len=240) :: molcas_path = ' '
 character(len=240) :: molpro_path = ' '
 character(len=240) :: orca_path = ' '
 character(len=240) :: psi4_path = ' '
 character(len=240) :: dalton_path = ' '
 character(len=240) :: bdf_path = ' '

 character(len=15) :: method = ' ' ! model chemistry, theoretical method
 character(len=21) :: basis = ' '  ! basis set (gen and genecp supported)
 character(len=21) :: RIJK_bas = 'NONE' ! cc-pVTZ/JK, def2/JK, etc for CASSCF
 character(len=21) :: RIC_bas  = 'NONE' ! cc-pVTZ/C, def2-TZVP/C, etc for NEVPT2
 character(len=21) :: F12_cabs = 'NONE' ! F12 cabs
contains

subroutine get_molcas_path()
 implicit none
 integer :: i, fid, system

 i = SYSTEM("which pymolcas >mokit.pymolcas 2>&1")
 if(i /= 0) then
  molcas_path = 'mokit.pymolcas'
  call delete_file(TRIM(molcas_path))
  molcas_path = 'NOT FOUND'
  return
 end if

 open(newunit=fid,file='mokit.pymolcas',status='old',position='rewind')
 read(fid,'(A)',iostat=i) molcas_path
 close(fid,status='delete')

 if(i /= 0) then
  molcas_path = 'NOT FOUND'
 else
  if(LEN_TRIM(molcas_path) == 0) then
   molcas_path = 'NOT FOUND'
  else if(index(molcas_path,'no pymolcas') > 0) then
   molcas_path = 'NOT FOUND'
  end if
 end if
end subroutine get_molcas_path

 ! repalce variables like '$USER' in path into real path
subroutine replace_env_in_path(path)
 implicit none
 integer :: i, j, k
 character(len=100) :: str
 character(len=240) :: buf
 character(len=240), intent(inout) :: path

 i = INDEX(path, '$')
 if(i == 0) return

 j = INDEX(path(i+1:),'/')
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
end subroutine replace_env_in_path

 ! read paths of various programs from environment variables
subroutine read_program_path()
 implicit none
 integer :: i
 integer(kind=4) :: hostnm
 character(len=8) :: hostname
 character(len=24) :: data_string
 character(len=240), external :: get_mokit_root 

 write(6,'(A)') '------ Output of AutoMR of MOKIT(Molecular Orbital Kit) ------'
 write(6,'(A)') '       GitLab page: https://gitlab.com/jxzou/mokit'
 write(6,'(A)') '     Documentation: https://jeanwsr.gitlab.io/mokit-doc-mdbook'
 write(6,'(A)') '           Version: 1.2.6rc37 (2024-Aug-6)'
 write(6,'(A)') '       How to cite: see README.md or $MOKIT_ROOT/doc/'

 hostname = ' '
 data_string = ' '
 i = hostnm(hostname)
 call fdate(data_string)
 write(6,'(/,A)') 'HOST '//TRIM(hostname)//', '//TRIM(data_string)

 write(6,'(/,A)') 'Read program paths from environment variables:'
 !call getenv('MOKIT_ROOT', mokit_root)
 mokit_root = get_mokit_root()
 write(6,'(A)') 'MOKIT_ROOT  = '//TRIM(mokit_root)

 call get_gau_path(gau_path)
 call get_molcas_path()
 call check_molcas_is_omp(molcas_omp)
 call get_molpro_path(molpro_path)
 call get_orca_path(orca_path)
 call get_psi4_path(psi4_path)
 call get_dalton_path(dalton_path)
 if(TRIM(dalton_path) /= 'NOT FOUND') call check_dalton_is_mpi(dalton_mpi)
 call getenv('GMS', gms_path)
 call getenv('BDF', bdf_path)
 if(LEN_TRIM(gms_path) == 0) gms_path = 'NOT FOUND'
 if(LEN_TRIM(bdf_path) == 0) bdf_path = 'NOT FOUND'

 write(6,'(A)') 'gau_path    = '//TRIM(gau_path)
 write(6,'(A)') 'gms_path    = '//TRIM(gms_path)
 write(6,'(A)') 'orca_path   = '//TRIM(orca_path)
 write(6,'(A)') 'molpro_path = '//TRIM(molpro_path)
 write(6,'(A)') 'molcas_path = '//TRIM(molcas_path)
 write(6,'(A)') 'psi4_path   = '//TRIM(psi4_path)
 write(6,'(A)') 'dalton_path = '//TRIM(dalton_path)
 write(6,'(A)') 'bdf_path    = '//TRIM(bdf_path)
end subroutine read_program_path

 ! check whether GAMESS path exists
subroutine check_gms_path()
 implicit none
 integer :: i, fid
 character(len=240) :: buf
 logical :: alive

 inquire(file=TRIM(gms_path),exist=alive)
 if(.not. alive) then
  write(6,'(/,A)') 'ERROR in subroutine check_gms_path: rungms does not exist.'
  write(6,'(A)') 'gms_path='//TRIM(gms_path)
  stop
 end if

 open(newunit=fid,file=TRIM(gms_path),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:7) == 'set SCR') exit
 end do ! for while
 close(fid)

 i = INDEX(buf, '=')
 if(i > 0) then
  gms_scr_path = buf(i+1:)
 else
  write(6,'(/,A)') "ERROR in subroutine check_gms_path: '=' not found in path "&
                  //TRIM(buf)
  stop
 end if
 i = INDEX(gms_scr_path, '#')
 if(i > 0) gms_scr_path = gms_scr_path(1:i-1)
 call replace_env_in_path(gms_scr_path)
 write(6,'(A)') 'gms_scr_path = '//TRIM(gms_scr_path)
end subroutine check_gms_path

 subroutine parse_keyword()
  implicit none
  integer :: i, j, k, ifail, fid
  character(len=24) :: method0 = ' '
  character(len=240) :: buf = ' '
  character(len=1000) :: longbuf = ' '
  logical :: alive(2), alive1(5)

  open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
  call read_mem_nproc_route(fid, mem, nproc, buf)

  i = INDEX(buf,'/')
  if(i == 0) then
   write(6,'(/,A)') "ERROR in subroutine parse_keyword: no '/' symbol detected &
                    &in keyword line."
   write(6,'(A)') "The method and basis set must be specified via '/' symbol, e&
                  &.g. CASSCF/cc-pVDZ."
   close(fid)
   stop
  end if

  if(DBLE(mem) < DBLE(nproc)*0.8d0) then
   write(6,'(/,A)') 'ERROR in subroutine parse_keyword: please specify larger m&
                    &emory. Multi-'
   write(6,'(A)') 'reference calculations usually requires large memory.'
   stop
  end if
  if(INDEX(buf,'scan') > 0) rigid_scan = .true.

  j = INDEX(buf(1:i-1),' ', back=.true.)
  if(j == 0) then
   write(6,'(/,A)') 'ERROR in subroutine parse_keyword: syntax error detected in'
   write(6,'(A)') "the current line '"//TRIM(buf)//"'"
   stop
  end if
  method0 = buf(j+1:i-1)

  i = INDEX(method0, '('); j = INDEX(method0, ','); k = INDEX(method0, ')')
  alive = [(i/=0 .and. k/=0), (i==0 .and. k==0)]
  if(.not. (alive(1) .or. alive(2)) ) then
   write(6,'(/,A)') 'ERROR in subroutine parse_keyword: incomplete method speci&
                    &fied.'
   write(6,'(A)') 'method = '//TRIM(method0)
   stop
  end if

  if(i /= 0) then
   method = method0(1:i-1)

   select case(TRIM(method))
   case('mcpdft','mc-pdft','mrcisd','mrcisdt','sdspt2','mrmp2','ovbmp2','caspt3',&
        'caspt2','caspt2k','caspt2-k','nevpt3','nevpt2','casscf','dmrgscf', &
        'casci','dmrgci','ficmrccsd','mkmrccsd','mkmrccsd(t)','bwmrccsd', &
        'bwmrccsd(t)')
    read(method0(i+1:j-1),*) nacte_wish
    read(method0(j+1:k-1),*) nacto_wish
    if(nacte_wish<1 .or. nacto_wish<1 .or. nacte_wish/=nacto_wish) then
     write(6,'(A)') 'ERROR in subroutine parse_keyword: wrong number of active&
                   & electrons/orbitals specified.'
     write(6,'(2(A,I0))') 'nacte/nacto=', nacte_wish, '/', nacto_wish
     write(6,'(A)') 'Currently only NactE = NactO > 0 such as (6,6) supported.'
     stop
    end if
   case('gvb','bccc2b','bccc3b') ! e.g. GVB(6), BCCC2b(6) (it means GVB(6)-BCCC2b)
    if(j /= 0) then
     write(6,'(/,A)') 'ERROR in subroutine parse_keyword: GVB active space shou&
                      &ld be specified like GVB(3),'
     write(6,'(A)') 'where 3 is the number of pairs. Did you specify (6,6) like&
                    & CAS?'
     stop
    end if
    read(method0(i+1:k-1),*) npair_wish
    if(npair_wish < 0) then
     write(6,'(/,A)') 'ERROR in subroutine parse_keyword: wrong number of pairs&
                      & specified.'
     write(6,'(A,I0)') 'npair_wish=', npair_wish
     stop
    end if
   case default
    write(6,'(/,A)') 'ERROR in subroutine parse_keyword: unsupported method '//&
                     TRIM(method)
    stop
   end select
  else ! i = 0
   method = TRIM(method0)
  end if

  select case(TRIM(method))
  case('mcpdft','mc-pdft','mrcisd','mrcisdt','sdspt2','mrmp2','ovbmp2','caspt3',&
       'caspt2','caspt2k','caspt2-k','nevpt3','nevpt2','casscf','dmrgscf', &
       'casci','dmrgci','gvb','ficmrccsd','mkmrccsd','mkmrccsd(t)','bwmrccsd', &
       'bwmrccsd(t)','bccc2b','bccc3b')
   uno = .true.; gvb = .true.
  case default
   write(6,'(/,A)') "ERROR in subroutine parse_keyword: specified method '"//&
                     TRIM(method)//"' not supported."
   write(6,'(/,A)') 'All supported methods are GVB, CASCI, CASSCF, DMRGCI, &
                    &DMRGSCF, NEVPT2,'
   write(6,'(A)') 'NEVPT3, CASPT2, CASPT2K, CASPT3, MRMP2, OVBMP2, MRCISD, &
                  &MRCISDT, MCPDFT,'
   write(6,'(A)') 'FICMRCCSD, MkMRCCSD, MkMRCCSD(T), BWMRCCSD, BWMRCCSD(T), &
                  &BCCC2b, BCCC3b.'
   stop
  end select

  select case(TRIM(method))
  case('mcpdft','mc-pdft')
   mcpdft = .true.
   casscf = .true.
  case('mrcisd')
   mrcisd = .true.
   casscf = .true.
  case('mrcisdt')
   mrcisdt= .true.
   casscf = .true.
  case('sdspt2')
   sdspt2 = .true.
   casscf = .true.
  case('mrmp2')
   mrmp2 = .true.
   casscf = .true.
  case('ovbmp2')
   ovbmp2 = .true.
   casscf = .true.
  case('caspt2')
   caspt2 = .true.
   casscf = .true.
  case('caspt2k','caspt2-k')
   caspt2 = .true.; casscf = .true.
   caspt2k = .true.; caspt2_prog = 'orca'
  case('caspt3')
   caspt3 = .true.
   casscf = .true.
  case('nevpt2')
   nevpt2 = .true.
   casscf = .true.
  case('nevpt3')
   nevpt3 = .true.
   casscf = .true.
  case('casscf')
   casscf = .true.
  case('dmrgscf')
   dmrgscf = .true.
  case('casci')
   casci = .true.
  case('dmrgci')
   dmrgci = .true.
  ! 1~8 for FIC-MRCC/MkMRCCSD/MkMRCCSD(T)/BWMRCCSD/BWMRCCSD(T)/BCCC2b,3b,4b
  case('ficmrccsd','mkmrccsd','mkmrccsd(t)','bwmrccsd','bwmrcccsd(t)')
   mrcc = .true.; casscf = .true.
   select case(TRIM(method))
   case('ficmrccsd')
    mrcc_type = 1
   case('mkmrccsd')
    mrcc_type = 2
   case('mkmrccsd(t)')
    mrcc_type = 3; mrcc_prog = 'nwchem'
   case('bwmrccsd')
    mrcc_type = 4
   case('bwmrccsd(t)')
    mrcc_type = 5; mrcc_prog = 'nwchem'
   end select
  case('bccc2b')
   mrcc = .true.; mrcc_type = 6
  case('bccc3b')
   mrcc = .true.; mrcc_type = 7
  case('gvb')
  end select

  i = INDEX(buf,'/')
  j = i - 1 + index(buf(i+1:),' ')
  if(j == 0) j = LEN_TRIM(buf)
  basis = buf(i+1:j)
  write(6,'(/,2(A,I4))',advance='no') 'memory =', mem, 'GB, nproc =', nproc
  write(6,'(A)') ', method/basis = '//TRIM(method)//'/'//TRIM(basis)

  if(basis(1:5) == 'def2-') then
   write(6,'(/,A)') "ERROR in subroutine parse_keyword: 'def2-' prefix detected&
                    & in given basis set."
   write(6,'(A)') 'Basis set in Gaussian syntax should be like def2TZVP, not de&
                  &f2-TZVP.'
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

  if(npair_wish > 0) write(6,'(A,I0)') 'User specified GVB npair = ',npair_wish
  if(nacte_wish>0 .and. nacto_wish>0) write(6,'(2(A,I0))') 'User specified&
                            & CAS nacte/nacto = ',nacte_wish,'/',nacto_wish

  i = INDEX(buf, 'guess=')
  if(i > 0) then
   write(6,'(/,A)') "ERROR in subroutine parse_keyword: 'guess=' syntax not sup&
                    &ported in automr."
   write(6,'(A)') "You can use 'guess()' syntax instead."
   stop
  end if

  i = INDEX(buf, 'guess(')
  if(i > 0) then
   frag_guess = .true.
   j = INDEX(buf(i+6:),'='); k = INDEX(buf(i+6:),')')
   if(j*k == 0)then
    write(6,'(A)') "ERROR in subroutine parse_keyword: 'guess(fragment=N)' synt&
                   &ax is wrong in file "//TRIM(gjfname)
    close(fid)
    stop
   end if
   read(buf(j+i+6:k+i+4),*) nfrag
  end if

  read(fid,'(A)') buf ! skip a blank line
  read(fid,'(A)') buf ! Title Card, the 1st line of keywords
  call lower(buf)
  ! Note: this lower() is not from string_manipulate.f90, but defined in this
  ! module

  if(buf(1:6) /= 'mokit{') then
   write(6,'(/,A)') "ERROR in subroutine parse_keyword: 'mokit{' not detected i&
                    &n file "//TRIM(gjfname)
   write(6,'(A)') "Syntax error. You must put 'mokit{' in leading position of t&
                  &he Title Card line."
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
    write(6,'(/,A)') 'ERROR in subroutine parse_keyword: end-of-file detected.'
    write(6,'(A)') 'The provided .gjf file may be incomplete.'
    close(fid)
    stop
   end if
  end if

  close(fid)
  ! now all keywords are stored in longbuf

  write(6,'(/,A)') 'Keywords in MOKIT{} are merged and shown as follows:'
  write(6,'(A)') TRIM(longbuf)

  alive1(1:3) = [(INDEX(longbuf,'hf_prog')>0), (INDEX(longbuf,'readuhf')>0), &
                 (INDEX(longbuf,'readno')>0)]
  if(alive1(1) .and. (alive1(2) .or. alive1(3))) then
   write(6,'(/,A)') "ERROR in subroutine parse_keyword: 'HF_prog' is useless wh&
                    &en 'readuhf'"
   write(6,'(A)') "or 'readno' is specified."
   stop
  end if

  alive1(1:5) = [(INDEX(longbuf,'caspt2_prog')/=0), (INDEX(longbuf,'nevpt2_prog')/=0),&
                 (INDEX(longbuf,'mrcisd_prog')/=0), (INDEX(longbuf,'mrmp2_prog')/=0), &
                 (INDEX(longbuf,'mcpdft_prog')/=0)]
  if(COUNT(alive1(1:5) .eqv. .true.) > 1) then
   write(6,'(/,A)') "ERROR in subroutine parse_keyword: more than one keyword o&
                    &f 'caspt2_prog',"
   write(6,'(A)') "'nevpt2_prog', 'mrmp2_prog', 'mrcisd_prog', 'mcpdft_prog' ar&
                  &e detected. Only one"
   write(6,'(A)') "can be specified in a job."
   stop
  end if

  alive1(1:4)= [(INDEX(longbuf,'casci_prog')/=0),(INDEX(longbuf,'casscf_prog')/=0),&
                (INDEX(longbuf,'dmrgci_prog')/=0),(INDEX(longbuf,'dmrgscf_prog')/=0)]
  if(alive1(1) .and. alive1(2)) then
   write(6,'(/,A)') 'ERROR in subroutine parse_keyword: both CASCI_prog and CAS&
                   &SCF_prog are detected.'
   write(6,'(A)') 'Only one can be specified in a job.'
   stop
  end if

  if(alive1(3) .and. alive1(4)) then
   write(6,'(/,A)') 'ERROR in subroutine parse_keyword: both DMRGCI_prog and DM&
                    &RGSCF_prog are detected.'
   write(6,'(A)') 'Only one can be specified in a job.'
   stop
  end if

  if(casscf .and. (alive1(1).or.alive1(3))) then
   write(6,'(/,A)') 'ERROR in subroutine parse_keyword: CASSCF activated, but y&
                    &ou specify the'
   write(6,'(A)') 'CASCI_prog or DMRGCI_prog. You should specify CASSCF_prog or&
                  & DMRGSCF_prog.'
   stop
  end if

  if(casci .and. (alive1(2).or.alive1(4))) then
   write(6,'(/,A)') 'ERROR in subroutine parse_keyword: CASCI activated, but yo&
                    &u specify the'
   write(6,'(A)') 'CASSCF_prog or DMRGSCF_prog. You should specify CASCI_prog o&
                  &r DMRGCI_prog.'
   stop
  end if

  if(dmrgscf .and. alive1(3)) then
   write(6,'(/,A)') 'ERROR in subroutine parse_keyword: DMRG-CASSCF activated,&
                   & but you specify the DMRGCI_prog.'
   stop
  end if

  if(dmrgci .and. alive1(4)) then
   write(6,'(/,A)') 'ERROR in subroutine parse_keyword: DMRG-CASCI activated,&
                   & but you specify the DMRGSCF_prog.'
   stop
  end if

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
    readuhf = .true.
    read(longbuf(j+1:i-1),*) hf_fch
   case('readno')
    mo_rhf = .true.
    readno = .true.
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
   case('ist')       ! the i-th strategy
    read(longbuf(j+1:i-1),*) ist
   case('ctrtype')   ! unconctracted-/ic-/FIC- MRCI
    read(longbuf(j+1:i-1),*) CtrType
   case('cionly')    ! do CASPT2/NEVPT2 after CASCI, skip CASSCF
    CIonly = .true.
    casscf = .false.; casci = .true.
    dmrgscf = .false.; dmrgci = .false.
   case('maxm')
    read(longbuf(j+1:i-1),*) maxM
   case('block_mpi')
    block_mpi = .true.
   case('hardwfn')
    hardwfn = .true.
   case('crazywfn')
    crazywfn = .true.
   case('hf_prog')
    hf_prog = longbuf(j+1:i-1)
   case('gvb_prog')
    gvb_prog = longbuf(j+1:i-1)
   case('casci_prog')
    casci_prog = longbuf(j+1:i-1)
   case('casscf_prog')
    casscf_prog = longbuf(j+1:i-1)
   case('dmrgci_prog')
    dmrgci_prog = longbuf(j+1:i-1)
   case('dmrgscf_prog')
    dmrgscf_prog = longbuf(j+1:i-1)
   case('caspt2_prog')
    caspt2_prog = longbuf(j+1:i-1)
   case('nevpt2_prog')
    nevpt2_prog = longbuf(j+1:i-1)
    if(TRIM(nevpt2_prog) == 'bdf') FIC = .true.
   case('mrmp2_prog')
    mrmp2_prog = longbuf(j+1:i-1)
   case('mrcisd_prog')
    mrcisd_prog = longbuf(j+1:i-1)
   case('mrcisdt_prog')
    mrcisdt_prog = longbuf(j+1:i-1)
   case('mcpdft_prog')
    mcpdft_prog = longbuf(j+1:i-1)
   case('mrcc_prog')
    mrcc_prog = longbuf(j+1:i-1)
   case('gvb_conv')
    GVB_conv = longbuf(j+1:i-1)
    c_gvb_conv = .true.
   case('skip_uno')
    read(longbuf(j+1:i-1),*) nskip_uno
   case('root')
    read(longbuf(j+1:i-1),*) iroot
    ss_opt = .true. ! State-specific orbital optimization
   case('xmult')
    read(longbuf(j+1:i-1),*) xmult
   case('force')
    force = .true.
   case('charge')
    bgchg = .true.
   case('inherit')
    inherit = .true.
   case('ri')
    RI = .true.
   case('rijk_bas')
    RIJK_bas = longbuf(j+1:i-1)
   case('ric_bas')
    RIC_bas = longbuf(j+1:i-1)
   case('f12')
    F12 = .true.; RI = .true.
   case('f12_cabs')
    F12_cabs = longbuf(j+1:i-1)
   case('dlpno')
    DLPNO = .true.; RI = .true.; FIC = .true.
   case('fic')
    FIC = .true.
   case('otpdf')
    otpdf = longbuf(j+1:i-1)
   case('on_thres')
    read(longbuf(j+1:i-1),*) on_thres
   case('uno_thres')
    read(longbuf(j+1:i-1),*) uno_thres
   case('npair') ! numbers of pairs for non-GVB calculations
    if(npair_wish /= 0) then
     write(6,'(/,A)') 'ERROR in subroutine parse_keyword: npair is specified by&
                      &more than once.'
     write(6,'(A)') 'Please check your input file.'
     stop
    end if
    read(longbuf(j+1:i-1),*) npair_wish
   case('nmr')
    nmr = .true.
   case('icss')
    ICSS = .true.; nmr = .true.
   case('excludexh')
    excludeXH = .true.
   case('onlyxh')
    excludeXH = .true.; onlyXH = .true.
   case('mixed_spin')
    mixed_spin = .true.
   case('nstates')
    read(longbuf(j+1:i-1),*) nstate ! SA-CASSCF (ground state not included)
    gvb = .false.; casscf = .false.; sa_cas = .true.; excited = .true.
    eist = 1
   case('tdhf')
    TDHF = .true.
   case('qd')
    QD = .true.
   case('fcgvb')
    fcgvb = .true.; c_fcgvb = .true.
   case('nofcgvb')
    fcgvb = .false.
   case('hfonly')
    HFonly = .true.
   case('nodmrgno')
    dmrg_no = .false.
   case('locdocc')
    LocDocc = .true.
   case default
    write(6,'(/,A)') "ERROR in subroutine parse_keyword: keyword '"//longbuf(1:j-1)&
                    //"' not recognized in {}."
    stop
   end select

   ! delete the keyword which has been specified
   longbuf(1:i) = ' '
   longbuf = ADJUSTL(longbuf)
   if(LEN_TRIM(longbuf) == 0) exit
  end do ! for while

  dkh2_or_x2c = (DKH2 .or. X2C)
  if(nstate > 999) then
   write(6,'(/,A)') 'ERROR in subroutine parse_keyword: too large nstates.'
   write(6,'(A)') 'The paramter Nstates must be <=999.'
   stop
  end if

  if(readrhf .or. readuhf .or. readno) then
   if(frag_guess) then
    write(6,'(/,A)') 'ERROR in subroutine parse_keyword: frag_guess can only be&
                    & used when none'
    write(6,'(A)') 'of readrhf/readuhf/readno is used.'
    stop
   end if
   skiphf = .true.
   call check_sanity_of_provided_fch(DKH2, X2C, hf_fch)
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
  case(6) ! GVB/STO-6G GVB -> GVB/target basis set
   gvb = .true.; skiphf = .true.
  case default
   write(6,'(/,A)') "ERROR in subroutine parse_keyword: the parameter 'ist' is &
                    & out of range."
   write(6,'(A)') 'Only 1~6 are allowed. See MOKIT manual $4.4.4 for details of&
                  & ist.'
   stop
  end select

  if(.not. mcpdft) otpdf = 'NONE'

  dyn_corr = (caspt2 .or. nevpt2 .or. mrmp2 .or. mrcisd .or. mrcisdt .or. &
              mcpdft .or. caspt3 .or. nevpt3)
  if(RI) call determine_auxbas(basis,RIJK_bas, dyn_corr,RIC_bas, F12,F12_cabs)
  call prt_strategy()
 end subroutine parse_keyword

subroutine prt_strategy()
 implicit none
 write(6,'(/,A,I0)') 'No. Strategy = ', ist

 write(6,'(5(A,L1,3X))') 'readRHF = ', readrhf, 'readUHF = ', readuhf,&
      'readNO  = ', readno, 'skipHF  = ', skiphf, 'Cart    = ', cart

 write(6,'(5(A,L1,3X))') 'Vir_Proj= ',vir_proj, 'UNO     = ', uno    ,&
      'GVB     = ', gvb   , 'CASCI   = ',  casci, 'CASSCF  = ', casscf

 write(6,'(5(A,L1,3X))') 'DMRGCI  = ',  dmrgci, 'DMRGSCF = ', dmrgscf,&
      'CASPT2  = ', caspt2, 'CASPT2K = ',caspt2k, 'CASPT3  = ', caspt3

 write(6,'(5(A,L1,3X))') 'MRMP2   = ',   mrmp2, 'OVBMP2  = ',  ovbmp2,&
      'SDSPT2  = ', sdspt2, 'MRCISD  = ', mrcisd, 'MRCISDT = ', MRCISDT

 write(6,'(5(A,L1,3X))') 'NEVPT2  = ',  nevpt2, 'NEVPT3  = ',  nevpt3,&
      'MCPDFT  = ', mcpdft, 'MRCC    = ',   mrcc, 'CIonly  = ', CIonly

 write(6,'(5(A,L1,3X))') 'dyn_corr= ',dyn_corr, 'DKH2    = ',    DKH2,&
      'X2C     = ',    X2C, 'RI      = ',     RI, 'FIC     = ', FIC

 write(6,'(5(A,L1,3X))') 'DLPNO   = ',   DLPNO, 'F12     = ',     F12,&
      'HardWFN = ',hardwfn, 'CrazyWFN= ',crazywfn, 'OnlyXH  = ', onlyXH

 write(6,'(5(A,L1,3X))') 'BgCharge= ',   bgchg, 'Ana_Grad= ',  force,&
      'Pop     = ',    pop, 'NMR     = ',    nmr, 'ICSS    = ', ICSS

 write(6,'(5(A,L1,3X))') 'TDHF    = ',    TDHF, 'SA_CAS  = ',  sa_cas,&
      'Excited = ',excited, 'QD      = ',     QD, 'SOC     = ', SOC

 write(6,'(A,L1,2X,3(A,L1,3X),A)') 'Mixed_Spin=',Mixed_Spin, 'RigidScan=',&
      rigid_scan, 'RelaxScan=',relaxed_scan, 'Inherit = ',inherit, &
      'GVB_conv= '//TRIM(GVB_conv)

 write(6,'(A,I2,3X,2(A,I1,3X),A,I5,3X,A,L1)') 'Skip_UNO=', nskip_uno, &
      'CtrType = ', CtrType, 'MRCC_type=',mrcc_type, 'MaxM =', maxM,&
      'excludeXH=', excludeXH

 write(6,'(A,F7.5,1X,A,F7.5)') 'LocalM  = '//TRIM(localm)//'  ON_thres= ',&
      on_thres, 'OtPDF='//TRIM(otpdf)//'  UNO_thres= ', uno_thres

 write(6,'(A)',advance='no') 'RIJK_bas='//TRIM(RIJK_bas)//' RIC_bas='//&
      TRIM(RIC_bas)//'  F12_cabs='//TRIM(F12_cabs)//' HF_fch='

 if(skiphf) then
  write(6,'(A)') TRIM(hf_fch)
 else
  write(6,'(A)') 'NONE'
 end if
end subroutine prt_strategy

subroutine check_kywd_compatible()
 implicit none
 integer :: i
 logical :: alive(3)
 character(len=10) :: cas_prog
 character(len=43), parameter :: error_warn = 'ERROR in subroutine check_kywd_compatible: '

 write(6,'(/,A)') 'Check if the keywords are compatible with each other...'

 if(readrhf .or. readuhf .or. readno) then
  if(ist == 6) then
   write(6,'(/,A)') 'ERROR in subroutine check_kywd_compatible: ist=6 is not &
                    &compatible with any'
   write(6,'(A)') 'keyword of readrhf/readuhf/readno.'
   stop
  end if
  call check_cart_in_fch(hf_fch, cart)
 end if

 if(on_thres<0d0 .or. on_thres>1d0) then
  write(6,'(/,A)') error_warn//'ON_thres must be in [0.0,1.0].'
  write(6,'(A,E12.5)') 'Your input ON_thres=', on_thres
  stop
 end if

 if(uno_thres<0d0 .or. uno_thres>1d0) then
  write(6,'(/,A)') error_warn//'uno_thres must be in [0.0,1.0].'
  write(6,'(A,E12.5)') 'Your input uno_thres=', uno_thres
  stop
 end if

 if(iroot < 0) then
  write(6,'(/,A)') error_warn//'iroot must be non-negative.'
  write(6,'(A,I0)') 'Your input iroot=', iroot
  stop
 else if(iroot > 0) then ! State-Specific CASSCF
  if(nstate > 0) then
   write(6,'(/,A)') error_warn//'you can specify only one of Nstates and Root.'
   stop
  end if
  if(casci) then
   write(6,'(/,A)') error_warn//'Root is expected to be used in CASSCF, not CASCI.'
   stop
  end if
 end if

 if(nmr) then
  if((.not.casscf) .and. iroot==0) then
   write(6,'(/,A)') error_warn//'NMR is supposed to be used with the'
   write(6,'(A)') 'CASSCF method. But neither CASSCF nor Root=N is specified.'
   stop
  end if
  if(TRIM(dalton_path) == 'NOT FOUND') then
   write(6,'(/,A)') error_warn//'it seems Dalton is not installed.'
   stop
  else
   call check_exe_exist(dalton_path)
  end if
 end if

 if(icss .and. dalton_mpi) then
  write(6,'(/,A)') error_warn//'you are using MPI version of Dalton.'
  write(6,'(A)') 'Please use MKL version of Dalton for CASSCF ICSS computations&
                 &, which is'
  write(6,'(A)') 'supposed to be faster.'
  stop
 end if

 if(DKH2 .and. X2C) then
  write(6,'(/,A)') error_warn//"'DKH2' and 'X2C' cannot both be activated."
  stop
 end if

 if(DKH2 .and. TRIM(hf_prog)=='pyscf') then
  write(6,'(/,A)') 'Warning: DKH2 not supported in PySCF. HF in PySCF will use&
                   & X2C Hamiltonian'
  write(6,'(A)') 'instead. MOs obtained by these two Hamiltonians are usually&
                & very similar.'
  write(6,'(A)') 'If you want to use DKH2 during HF calculations, you can spec&
                 &ify HF_prog=PSI4/ORCA.'
 end if

 if(X2C .and. .not.(TRIM(hf_prog)=='pyscf' .or. TRIM(hf_prog)=='psi4') .and. &
    (.not.skiphf)) then
  write(6,'(/,A)') 'Warning: X2C is activated but currently HF_prog is not PySC&
                   &F/PSI4. Switching'
  write(6,'(A)') 'to HF_prog=PySCF automatically.'
  hf_prog = 'pyscf'
 end if

 if(casci .or. casscf) then
  if(casci) then
   cas_prog = casci_prog
  else
   cas_prog = casscf_prog
  end if
 end if

 if(RI) then
  if(DKH2 .or. X2C) then
   write(6,'(/,A)') error_warn//'currently RI cannot be applied in DKH2/X2C co&
                   &mputations.'
   stop
  end if

  if(.not. (casci .or. casscf)) then
   write(6,'(/,A)') error_warn//'RI activated. But neither CASCI nor CASSCF is&
                   & invoked.'
   stop
  end if
  select case(cas_prog)
  case('pyscf','orca','openmolcas','psi4','molpro')
  case default
   write(6,'(/,A)') error_warn//'CASCI/CASSCF with RI-JK is not supported'
   write(6,'(A)') 'for CASCI_prog or CASSCF_prog='//TRIM(cas_prog)
   write(6,'(A)') 'You should specify CASCI_prog or CASSCF_prog=PySCF/&
                  &ORCA/OpenMolcas/Molpro/PSI4.'
   stop
  end select
 end if

 if(F12) then
  if(.not. RI) then
   write(6,'(/,A)') error_warn//'F12 must be combined with RI. But RI is set'
   write(6,'(A)') 'to be False. Impossible.'
   stop
  end if
  if(.not. (nevpt2 .or. mrcisd)) then
   write(6,'(/,A)') error_warn//'F12 can only be used in NEVPT2 or MRCISD.'
   write(6,'(A)') 'But neither of NEVPT2/MRCISD is specified.'
   stop
  end if
  if(nevpt2) then
   if(nevpt2_prog /= 'orca') then
    write(6,'(/,A)') error_warn//'NEVPT2-F12 is only supported with ORCA.'
    write(6,'(A)') 'But currently NEVPT2_prog='//TRIM(nevpt2_prog)
    stop
   end if
   if(.not. FIC) then
    write(6,'(/,A)') error_warn//'SC-NEVPT2-F12 is not supported in ORCA.'
    write(6,'(A)') 'Only FIC-NEVPT2-F12 is supported. You need to add&
                    & keyword FIC in mokit{}.'
    stop
   end if
  end if
  if(mrcisd .and. mrcisd_prog/='molpro') then
   write(6,'(/,A)') error_warn//'MRCISD-F12 is only supported with Molpro.'
   write(6,'(A)') 'But currently MRCISD_prog='//TRIM(mrcisd_prog)
   stop
  end if
 end if

 alive = .false. ! remember to initialize
 select case(TRIM(casci_prog))
 case('gaussian','gamess','orca')
  alive(1) = .true.
 end select
 select case(TRIM(casscf_prog))
 case('gaussian','gamess','orca')
  alive(2) = .true.
 end select
 alive(3) = ((casci .and. alive(1)) .or. (casscf .and. alive(2)))
 if(X2C .and. alive(3)) then
  write(6,'(/,A)') error_warn//'CASCI/CASSCF with Gaussian/GAMESS/ORCA is&
                  & incompatible with'
  write(6,'(A)') 'X2C. You can use PySCF/Molpro/OpenMolcas.'
  stop
 end if

 alive = .false. ! remember to initialize
 select case(TRIM(casci_prog))
 case('pyscf','bdf')
  alive(1) = .true.
 end select
 select case(TRIM(casscf_prog))
 case('pyscf','bdf')
  alive(2) = .true.
 end select
 alive(3) = ( ((dmrgci .or. casci) .and. alive(1)) .or. &
              ((dmrgscf .or. casscf) .and. alive(2)) )
 if(DKH2 .and. alive(3)) then
  write(6,'(/,A)') error_warn//'CASCI/CASSCF with DKH2 is not'
  write(6,'(A)') 'supported by PySCF/BDF.'
  write(6,'(A)') 'For CASCI/CASSCF, you can use CASCI_prog/CASSCF_prog=Molpro/&
                 &OpenMolcas/GAMESS'
  write(6,'(A)') 'ORCA/Gaussian.'
  write(6,'(A)') 'For DMRG-CASCI/DMRG-CASSCF, you can use CASCI_prog/CASSCF_pr&
                 &og=OpenMolcas.'
  stop
 end if

 select case(dmrgci_prog)
 case('pyscf', 'openmolcas')
 case default
  write(6,'(/,A)') error_warn//'currently DMRG-CASCI is only supported by PySCF'
  write(6,'(A)') 'or OpenMolcas. Wrong DMRGCI_prog='//TRIM(dmrgci_prog)
  stop
 end select

 select case(dmrgscf_prog)
 case('pyscf', 'openmolcas')
 case default
  write(6,'(/,A)') error_warn//'currently DMRG-CASSCF is only supported by PySCF'
  write(6,'(A)') 'or OpenMolcas. Wrong DMRGSCF_prog='//TRIM(dmrgscf_prog)
  stop
 end select

 alive(1) = (.not.(casci .or. casscf .or. dmrgci .or. dmrgscf) .and. gvb)
 if(alive(1)) then
  ! if (post-)GVB calculation is requested, set to GAMESS default
  ! if (post-)CASSCF calculation is requested, use 5d-4 and FcGVB=.T.
  if(.not. c_gvb_conv) GVB_conv = '1d-5'
  if(.not. c_fcgvb) fcgvb = .false.
  if(X2C) then
   write(6,'(/,A)') error_warn//'GVB with GAMESS is incompatible with X2C.'
   stop
  end if
 end if

 if(hardwfn .and. crazywfn) then
  write(6,'(/,A)') error_warn//"'hardwfn' or 'crazywfn' cannot both be activated."
  stop
 end if

 if(.not. (TRIM(localm)=='pm' .or. TRIM(localm)=='boys')) then
  write(6,'(/,A)') error_warn//"only 'PM' or 'Boys' localization is supported."
  write(6,'(A)') 'Wrong LocalM='//TRIM(localm)
  stop
 end if

 alive = [readrhf, readuhf, readno]
 i = COUNT(alive .eqv. .true.)
 if(i > 1) then
  write(6,'(/,A)') error_warn//"more than one of 'readrhf',"
  write(6,'(A)') "'readuhf', and 'readno' are specified. These three keywords &
                 &are mutually exclusive."
  stop
 end if

 if(readrhf .and. .not.(ist==3 .or. ist==4)) then
  write(6,'(/,A)') error_warn//"'readrhf' is only compatible with ist=3 or 4."
  stop
 end if

 if(.not.readrhf .and. (ist==3 .or. ist==4)) then
  write(6,'(/,A)') error_warn//"ist=3 or 4 specified, it must be used combined &
                  &with 'readrhf'."
  stop
 end if

 if(readuhf .and. .not.(ist==1 .or. ist==2)) then
  write(6,'(/,A)') error_warn//"'readuhf' is only compatible with ist=1 or 2."
  stop
 end if

 if(.not.readuhf .and. (ist==1 .or. ist==2)) then
  write(6,'(/,A)') error_warn//"ist=1 or 2 specified, it must be used combined &
                  &with 'readuhf'."
  stop
 end if

 if(readno .and. ist/=5) then
  write(6,'(/,A)') error_warn//"'readno' is only compatible with ist=5."
  stop
 end if

 if(.not.readno .and. ist==5) then
  write(6,'(/,A)') error_warn//"ist=5 specified, it must be used combined with &
                  &'readno'."
  stop
 end if

 if(CIonly .and. (.not.caspt2) .and. (.not.nevpt2) .and. (.not.mrcisd) .and. &
    (.not. mcpdft) .and. (.not.caspt3) .and. (.not.mrcc)) then
  write(6,'(/,A)') error_warn//"keyword 'CIonly' can only be used in"
  write(6,'(A)') 'CASPT2/CASPT3/NEVPT2/MRCISD/MC-PDFT/MRCC computations. But&
                 & none of them is specified.'
  stop
 end if

 if(CIonly .and. TRIM(nevpt2_prog)=='bdf') then
  write(6,'(/,A)') error_warn//'currently CASCI-NEVPT2 is not sopported in BDF &
                  &program.'
  write(6,'(A)') 'You may use NEVPT2_prog=PySCF, Molpro, ORCA or OpenMolcas.'
  stop
 end if

 select case(TRIM(gvb_prog))
 case('gamess','gaussian','qchem')
 case default
  write(6,'(/,A)') error_warn//'only GAMESS/Gaussian/QChem is supported for the'
  write(6,'(A)') 'GVB computation. User specified GVB program cannot be ident&
                 &ified: '//TRIM(gvb_prog)
  stop
 end select

 if(mcpdft .and. TRIM(mcpdft_prog)=='gamess' .and. bgchg) then
  write(6,'(/,A)') error_warn
  write(6,'(A)') 'Currently MC-PDFT with point charges is incompatible with GAMESS.'
  write(6,'(A)') 'You can use PySCF/OpenMolcas.'
  stop
 end if

 select case(TRIM(mcpdft_prog))
 case('pyscf','openmolcas','gamess')
 case default
  write(6,'(/,A)') error_warn
  write(6,'(A)') 'User specified MC-PDFT program cannot be identified: '&
                //TRIM(mcpdft_prog)
 end select

 select case(TRIM(casci_prog))
 case('gaussian','gamess','openmolcas','pyscf','orca','molpro','bdf','psi4','dalton')
 case default
  write(6,'(/,A)') error_warn
  write(6,'(A)') 'User specified CASCI program cannot be identified: '//TRIM(casci_prog)
  stop
 end select

 select case(TRIM(casscf_prog))
 case('gaussian','gamess','openmolcas','pyscf','orca','molpro','bdf','psi4','dalton')
 case default
  write(6,'(/,A)') error_warn
  write(6,'(A)') 'User specified CASSCF program cannot be identified: '//TRIM(casscf_prog)
  stop
 end select

 select case(TRIM(mrcisd_prog))
 case('gaussian', 'orca', 'openmolcas', 'molpro','psi4','dalton','gamess')
 case default
  write(6,'(/,A)') error_warn
  write(6,'(A)') 'User specified MRCISD program is not supported: '//&
                  TRIM(mrcisd_prog)
  stop
 end select

 select case(TRIM(mrcisdt_prog))
 case('openmolcas','dalton','psi4','gamess')
 case default
  write(6,'(/,A)') error_warn
  write(6,'(A)') 'User specified MRCISDT program is not supported: '//&
                  TRIM(mrcisdt_prog)
  stop
 end select

 if(mrcisdt) then
  CtrType = 1
  write(6,'(/,A)') 'Only uncontracted MRCISDT is supported. Automatically &
                   &setting CtrType=1.'
 end if

 if(mrcisd) then
  select case(CtrType)
  case(1) ! uncontracted MRCISD
   if(mrcisd_prog == 'molpro') then
    write(6,'(/,A)') error_warn
    write(6,'(A)') 'Currently (uc-)MRCISD cannot be done with Molpro.'
    stop
   end if
   select case(TRIM(mrcisd_prog))
   case('gaussian','orca','gamess')
    if(X2C) then
     write(6,'(/,A)') error_warn
     write(6,'(A)') 'MRCISD in Gaussian/ORCA/GAMESS is incompatible with X2C.'
     stop
    end if
   end select
  case(2) ! ic-MRCISD
   select case(TRIM(mrcisd_prog))
   case('openmolcas', 'molpro')
   case default
    write(6,'(/,A)') error_warn
    write(6,'(A)') 'The ic-MRCISD are only supported by OpenMolcas and Molpro.&
                   & But you specify'
    write(6,'(A)') 'MRCISD_prog='//TRIM(mrcisd_prog)
    stop
   end select
  case(3) ! FIC-MRCISD
   if(nacto_wish > 15) then
    mrcisd_prog = 'pyscf'
    write(6,'(/,A)') 'Large size of active space. MRCISD_prog=PySCF is automati&
                     &cally set. PySCF+Block2'
    write(6,'(A)') 'would be called to perform the DMRG-FIC-MRCISD calculation.'
   else if(nacto_wish > 0) then
    mrcisd_prog = 'orca'
    write(6,'(/,A)') 'Small/Medium size of active space. MRCISD_prog=ORCA is au&
                     &tomatically set.'
   end if
   ! It is possible that nacto_wish=0 here, so the following 'select' is needed.
   select case(TRIM(mrcisd_prog))
   case('orca')
    if(X2C) then
     write(6,'(/,A)') error_warn
     write(6,'(A)') 'FIC-MRCISD with ORCA is incompatible with X2C currently.'
     stop
    end if
   case('pyscf')
    ! the check of DKH2 has been done in previous CASCI/CASSCF
   case default
    write(6,'(/,A)') error_warn
    write(6,'(A)') 'MRCISD_prog for FIC-MRCISD/DMRG-FIC-MRCISD method is ORCA/&
                   &PySCF, respectively.'
    write(6,'(A)') 'But current MRCISD_prog='//TRIM(mrcisd_prog)
    stop
   end select
  case default
   write(6,'(/,A)') error_warn//'invalid CtrType.'
   write(6,'(/,A)') 'MRCISD has many variants, please read Section 4.4.17 MRCI&
                    &SD_prog in MOKIT manual.'
   write(6,'(A)') 'You need to specify CtrType=1/2/3 for uncontracted/ic-/FIC-&
                  &MRCISD, respectively.'
   stop
  end select

  select case(TRIM(mrcisd_prog))
  case('gaussian','psi4','dalton','gamess')
   if(CtrType /= 1) then
    write(6,'(/,A)') error_warn
    write(6,'(A)') 'Gaussian/PSI4/Dalton/GAMESS only supports uncontracted MRCISD,'
    write(6,'(A,I0)') 'i.e. CtrType=1. But you specify CtrType=',CtrType
    stop
   end if
  end select

  if(TRIM(mrcisd_prog)=='orca' .and. cart) then
   write(6,'(/,A)') error_warn//'conflict settings.'
   write(6,'(A)') 'ORCA is set as the MRCI_prog, and it only supports spherical&
                  & harmonic functions,'
   write(6,'(A)') "but Cart = True, you should delete the keyword 'cart', or pr&
                  &ovide a .fch file"
   write(6,'(A)') 'with spherical harmonic functions.'
   stop
  end if
 end if

 if((TRIM(casci_prog)=='orca' .or. TRIM(casscf_prog)=='orca') .and. cart) then
  write(6,'(A)') error_warn
  write(6,'(A)') 'ORCA is set as CASCI_prog/CASSCF_prog, and it only supports &
                 &spherical harmonic'
  write(6,'(A)') 'functions, but Cart=True. Please use another program or prov&
                 &ide a .fch file'
  write(6,'(A)') 'with spherical harmonic functions.'
  stop
 end if

 select case(TRIM(caspt2_prog))
 case('openmolcas', 'molpro','orca')
 case default
  write(6,'(/,A)') error_warn
  write(6,'(A)') 'Supported CASPT2_prog=OpenMolcas/Molpro/ORCA.'
  write(6,'(A)') 'User specified CASPT2 program cannot be identified: '//TRIM(caspt2_prog)
  stop
 end select

 select case(TRIM(nevpt2_prog))
 case('pyscf','molpro','openmolcas','orca','bdf')
 case default
  write(6,'(/,A)') error_warn
  write(6,'(A)') 'Supported NEVPT2_prog=PySCF/OpenMolcas/Molpro/ORCA/BDF.'
  write(6,'(A)') 'User specified NEVPT2 program cannot be identified: '//TRIM(nevpt2_prog)
  stop
 end select

 if(mrmp2_prog /= 'gamess') then
  write(6,'(/,A)') error_warn
  write(6,'(A)') 'Only MRMP2_prog=GAMESS is supported.'
  write(6,'(A)') 'User specified MRMP2 program cannot be identified: '//TRIM(mrmp2_prog)
  stop
 end if

 select case(TRIM(mrcc_prog))
 case('orca','nwchem')
 case default
  write(6,'(A)') error_warn
  write(6,'(A)') 'Currently the MRCC method is only supported by ORCA/NWChem. &
                 & But got MRCC_prog='//TRIM(mrcc_prog)
  stop
 end select

 if(force) then
  if(casscf .and. (.not.dyn_corr)) casscf_force = .true.
  if(caspt2) caspt2_force = .true.
  if(nevpt2) nevpt2_force = .true.
  if(mcpdft) mcpdft_force = .true.
  if(TRIM(casscf_prog) == 'psi4') then
   write(6,'(/,A)') error_warn
   write(6,'(A)') 'CASSCF analytical gradients are not supported in PSI4. Pleas&
                  &e use another CASSCF_prog.'
   stop
  end if
 end if

 if(casscf_force .and. cart .and. TRIM(casscf_prog)=='pyscf') then
  write(6,'(/,A)') error_warn//"current version of PySCF can only compute force"
  write(6,'(A)') 'using spherical harmonic basis fucntions.'
  stop
 end if

 if(mrmp2 .and. X2C) then
  write(6,'(/,A)') error_warn//'MRMP2 with GAMESS is incompatible with X2C.'
  write(6,'(A)') 'You can use the DKH2 Hamiltonian instead.'
  stop
 end if

 if(nevpt2) then
  if(DKH2 .and. (TRIM(nevpt2_prog)=='pyscf' .or. TRIM(nevpt2_prog)=='bdf')) then
   write(6,'(/,A)') error_warn//'NEVPT2 with DKH2 is not supported by PySCF or BDF.'
   write(6,'(A)') 'You can use NEVPT2_prog=Molpro or ORCA.'
   stop
  else if(X2C .and. TRIM(nevpt2_prog)=='orca') then
   write(6,'(/,A)') error_warn//'NEVPT2 with X2C is not supported by ORCA.'
   write(6,'(A)') 'You can use NEVPT2_prog=Molpro, OpenMolcas, ORCA or BDF.'
   stop
  end if
  if(TRIM(nevpt2_prog)=='bdf' .and. bgchg) then
   write(6,'(/,A)') error_warn//'NEVPT2 with BDF program is incompatible with'
   write(6,'(A)') 'background point charges. You can use NEVPT2_prog=Molpro or ORCA.'
   stop
  end if
  if(FIC .and. TRIM(nevpt2_prog)=='pyscf') then
   write(6,'(/,A)') error_warn//'FIC-NEVPT2 is not supported by PySCF.'
   write(6,'(A)') 'You can use NEVPT2_prog=Molpro,BDF,ORCA,OpenMolcas.'
   stop
  end if
  if(RI .and. TRIM(nevpt2_prog)=='openmolcas') then
   write(6,'(/,A)') error_warn//'RI not supported in DMRG-NEVPT2 using OpenMolcas.'
   stop
  end if
 end if

 if((sdspt2.or.nevpt3) .and. bgchg) then
  write(6,'(/,A)') error_warn//'SDSPT2 or NEVPT3 with BDF program is incompati-'
  write(6,'(A)') 'ble with background point charges.'
  stop
 end if

 if((DKH2 .or. X2C) .and. cart) then
  write(6,'(/,A)') error_warn//'relativistic calculations using Cartesian'
  write(6,'(A)') 'functions may cause severe numerical instability. Please use&
                 & spherical'
  write(6,'(A)') 'harmonic type basis set.'
  stop
 end if

 if(excludeXH .and. TRIM(gvb_prog)/='gamess') then
  write(6,'(/,A)') 'ERROR in subroutine check_kywd_compatible: the keyword excl&
                 &udeXH currently'
  write(6,'(A)') 'is supported only for GAMESS. But got GVB_prog='//TRIM(gvb_prog)
  stop
 end if

 write(6,'(A)') 'Check done. All keywords are compatible.'
 write(6,'(/,A)') REPEAT('-',79)
 write(6,'(A)') "Note: in any following output which starts with '$' symbol, &
                &it is the Shell co-"
 write(6,'(A)') '      mmand used to submit the corresponding job.'
 write(6,'(A)') REPEAT('-',79)
end subroutine check_kywd_compatible

! turn letters in buf into lower case, except those in symbol ''
subroutine lower(buf)
 implicit none
 integer :: i, j, k, k1, k2
 character(len=240), intent(inout) :: buf

 k = LEN_TRIM(buf)
 k1 = INDEX(buf,"'")
 k2 = INDEX(buf,"'",back=.true.)

 do i = 1, k, 1
  if(i>k1 .and. i<k2) cycle

  j = IACHAR(buf(i:i))
  if(j>64 .and. j<91) buf(i:i) = ACHAR(j+32)
 end do ! for i
end subroutine lower

 ! read background point charge(s) from .gjf file
 subroutine read_bgchg_from_gjf(no_coor)
  use mol, only: natom, nbgchg, bgcharge, ptchg_e, nuc_pt_e, nuc, coor
  implicit none
  integer :: i, fid, nblank, nblank0
  character(len=240) :: buf
  character(len=41), parameter :: error_warn='ERROR in subroutine read_bgchg_from_gjf: '
  logical, intent(in) :: no_coor

  nblank = 0
  if(no_coor) then ! no Cartesian Coordinates
   if(rigid_scan .or. relaxed_scan) then
    nblank0 = 3
   else
    nblank0 = 2
   end if
  else             ! there exists Cartesian Coordinates
   if(rigid_scan .or. relaxed_scan) then
    nblank0 = 4
   else
    nblank0 = 3
   end if
  end if

  open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(LEN_TRIM(buf) == 0) nblank = nblank + 1
   if(nblank == nblank0) exit
  end do ! for while

  if(i /= 0) then
   write(66,'(A)') error_warn//'wrong format of background point charges.'
   close(fid)
   stop
  end if

  nbgchg = 0

  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i/=0 .or. LEN_TRIM(buf)==0) exit
   nbgchg = nbgchg + 1
  end do ! for while

  if(nbgchg == 0) then
   write(6,'(A)') error_warn//'no background point charge(s) found'
   write(6,'(A)') 'in file '//TRIM(gjfname)
   close(fid)
   stop
  end if

  write(6,'(A,I0)') 'Background point charge specified: nbgchg = ', nbgchg
  allocate(bgcharge(4,nbgchg), source=0d0)

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
  write(6,'(A,F18.8,A)') 'Self energy of the charges =', ptchg_e, ' a.u.'

  i = INDEX(gjfname, '.gjf', back=.true.)
  chgname = gjfname(1:i-1)//'.chg'
  call write_charge_into_chg(nbgchg, bgcharge, chgname)

  call calc_nuc_pt_e(nbgchg, bgcharge, natom, nuc, coor, nuc_pt_e)
  write(6,'(A,F18.8,A)') 'Nuclei-charges interaction =',nuc_pt_e,' a.u.'
  write(6,'(A)') 'Note: these two energies are included in all electronic energ&
                 &ies below.'
 end subroutine read_bgchg_from_gjf

 ! check any readrhf/readuhf/readno in the given .gjf file
 function check_readfch(gjfname) result(has_readfch)
  implicit none
  integer :: fid
  character(len=240) :: buf
  character(len=240), intent(in) :: gjfname
  logical :: has_readfch
 
  has_readfch = .false.
  open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
 
  do while(.true.)
   read(fid,'(A)') buf
   if(LEN_TRIM(buf) == 0) exit
  end do ! for while
 
  read(fid,'(A)') buf
  close(fid)
 
  if(LEN_TRIM(buf) == 0) return
  call lower(buf)
 
  if(index(buf,'readrhf')>0 .or. index(buf,'readuhf')>0 .or. &
     index(buf,'readno')>0) has_readfch = .true.
 end function check_readfch

end module mr_keyword

! read mem, nproc and Route Section from an opened .gjf file
subroutine read_mem_nproc_route(fid, mem, nproc, buf)
 implicit none
 integer :: i, ifail
 integer, intent(in) :: fid
 integer, intent(out) :: mem, nproc
 character(len=240), intent(out) :: buf
 character(len=42), parameter :: warn_str = 'ERROR in subroutine read_mem_nproc&
                                             &_route: '
 buf = ' '
 do while(.true.)
  read(fid,'(A)',iostat=ifail) buf
  if(ifail /= 0) exit
  call lower(buf)
  if(buf(1:1) == '#') exit
  if(buf(1:5) == '%mem=') then
   i = INDEX(buf,'gb')
   if(i == 0) then
    write(6,'(/,A)') warn_str//' memory unit only GB is accepted.'
    close(fid)
    stop
   end if
   read(buf(6:i-1),*) mem
  else if(buf(1:6) == '%nproc') then
   i = INDEX(buf,'=')
   read(buf(i+1:),*) nproc
  else
   write(6,'(/,A)') warn_str//' only %mem and %nproc are allowed.'
   write(6,'(A)') 'Bot now got '//TRIM(buf)
   close(fid)
   stop
  end if
 end do ! for while

 if(ifail /= 0) then
  write(6,'(/,A)') warn_str//' end-of-file detected.'
  write(6,'(A)') 'The input file is problematic.'
  close(fid)
  stop
 end if
end subroutine read_mem_nproc_route

! If RIJK_bas, RIC_bas and F12_cabs are given, check validity
! If not given, automatically determine them
subroutine determine_auxbas(basis, RIJK_bas, dyn, RIC_bas, F12, F12_cabs)
 implicit none
 integer :: i, j
 character(len=21) :: basis1
 character(len=21), intent(in) :: basis
 character(len=21), intent(inout) :: RIJK_bas, RIC_bas, F12_cabs
 logical, intent(in) :: dyn ! dynamic correlation
 logical, intent(in) :: F12 ! F12

 if(TRIM(RIJK_bas) /= 'NONE') then
  call lower(RIJK_bas)
  if(INDEX(RIJK_bas, '/jk') == 0) then
   write(6,'(/,A)') REPEAT('-',79)
   write(6,'(A)') 'Warning in subroutine determine_auxbas: RI-JK auxiliary basi&
                  &s set does not'
   write(6,'(A)') "contain key '/JK'. Did you specify wrong auxiliary basis set&
                  &? Caution!"
   write(6,'(A)') REPEAT('-',79)
  end if
 end if

 if(dyn) then
  if(TRIM(RIC_bas) /= 'NONE') then
   call lower(RIC_bas)
   if(INDEX(RIC_bas, '/c') == 0) then
    write(6,'(/,A)') REPEAT('-',79)
    write(6,'(A)') 'Warning in subroutine determine_auxbas: dynamic correlation&
                   & computations'
    write(6,'(A)') "activated. But your provided RIC_bas does not include key '&
                   &/C'. Caution!"
    write(6,'(A)') REPEAT('-',79)
   end if
   select case(RIC_bas(1:5)) ! auto-switch Gaussian to ORCA syntax
   case('def2s','def2t','def2q')
    RIC_bas = RIC_bas(1:4)//'-'//TRIM(RIC_bas(5:))
   end select
  end if
 else   ! no dynamic correlation computation
  if(TRIM(RIC_bas) /= 'NONE') then
   write(6,'(/,A)') 'ERROR in subroutine determine_auxbas: no dynamic correlati&
                    &on computation is'
   write(6,'(A)') 'activated. But you provide RIC_bas='//TRIM(RIC_bas)//'.'
   stop
  end if
 end if

 if(F12) then
  if(TRIM(F12_cabs) /= 'NONE') then
   call lower(F12_cabs)
   if(INDEX(F12_cabs,'-f12') == 0) then
    write(6,'(/,A)') REPEAT('-',79)
    write(6,'(A)') 'Warning in subroutine determine_auxbas: F12 computation act&
                   &ivated. But your'
    write(6,'(A)') "provided basis set does not contain key '-F12'. Caution!"
    write(6,'(A)') REPEAT('-',79)
   end if
   if(INDEX(F12_cabs,'-CABS') == 0) then
    write(6,'(/,A)') REPEAT('-',79)
    write(6,'(A)') 'Warning in subroutine determine_auxbas: F12 computation act&
                   &ivated. But your'
    write(6,'(A)') "provided F12_cabs does not contain key '-cabs'. Caution!"
    write(6,'(A)') REPEAT('-',79)
   end if
  end if
 else
  if(TRIM(F12_cabs) /= 'NONE') then
   write(6,'(/,A)') 'ERROR in subroutine determine_auxbas: F12 not activated, b&
                    &ut you provide'
   write(6,'(A)') 'F12_cabs='//TRIM(F12_cabs)//'.'
   stop
  end if
 end if

 ! RI-J, RIJCOSX are not supported in AutoMR. Only RIJK for CASSCF and RI for
 ! MRPT2/MRCISD is supported.
 select case(TRIM(basis))
 case('def2sv(p)','def2svp','def2svpd','def2tzvp(-f)','def2tzvp','def2tzvpp',&
      'def2tzvpd','def2tzvppd','def2qzvpp','def2qzvpd','def2qzvppd')
  if(TRIM(RIJK_bas) == 'NONE') RIJK_bas = 'def2/JK'
  if(dyn .and. TRIM(RIC_bas)=='NONE') then
   if(basis(5:5) == 's') then
    RIC_bas = 'def2-SVP/C'
   else if(basis(5:5) == 'q') then
    RIC_bas = 'def2-QZVPP/C'
   else if(basis(1:9) == 'def2tzvpp') then
    RIC_bas = 'def2-TZVPP/C'
   else
    RIC_bas = 'def2-TZVP/C'
   end if
  end if

 case('ma-def2sv(p)','ma-def2svp','ma-def2tzvp(-f)','ma-def2tzvp','ma-def2tzvpp',&
      'ma-def2qzvpp')
  if(TRIM(RIJK_bas) == 'NONE') RIJK_bas = 'def2/JK'
  if(dyn .and. TRIM(RIC_bas)=='NONE') then
   if(basis(8:8) == 's') then
    RIC_bas = 'def2-SVPD/C'
   else if(basis(8:8) == 'q') then
    RIC_bas = 'def2-QZVPPD/C'
   else if(basis(1:12) == 'ma-def2tzvpp') then
    RIC_bas = 'def2-TZVPPD/C'
   else
    RIC_bas = 'def2-TZVPD/C'
   end if
  end if

 case('cc-pvdz','cc-pvtz','cc-pvqz','cc-pv5z','aug-cc-pvdz','aug-cc-pvtz',&
      'aug-cc-pvqz','aug-cc-pv5z','mar-cc-pv5z','apr-cc-pvqz','apr-cc-pv5z',&
      'may-cc-pvtz','may-cc-pvqz','may-cc-pv5z','jun-cc-pvdz','jun-cc-pvtz',&
      'jun-cc-pvqz','jun-cc-pv5z')
  if(TRIM(RIJK_bas) == 'NONE') RIJK_bas = TRIM(basis)//'/JK'
  if(dyn .and. TRIM(RIC_bas)=='NONE') RIC_bas = TRIM(basis)//'/C'
  if(F12) then
   if(basis(1:4) == 'aug-') then
    F12_cabs = TRIM(basis(5:))//'-F12-CABS'
   else
    F12_cabs = TRIM(basis)//'-F12-CABS'
   end if
  end if

 case('cc-pwcvdz','cc-pwcvtz','cc-pwcvqz','cc-pwcv5z','aug-cc-pwcvdz',&
      'aug-cc-pwcvtz','aug-cc-pwcvqz','aug-cc-pwcv5z')
  i = INDEX(basis, 'wc')
  basis1(1:i-1) = basis(1:i-1)
  basis1(i:) = basis(i+2:)
  if(TRIM(RIJK_bas) == 'NONE') RIJK_bas = TRIM(basis1)//'/JK'
  if(dyn .and. TRIM(RIC_bas)=='NONE') RIC_bas = TRIM(basis)//'/C'
  if(F12) then
   if(basis1(1:4) == 'aug-') then
    F12_cabs = TRIM(basis1(5:))//'-F12-CABS'
   else
    F12_cabs = TRIM(basis1)//'-F12-CABS'
   end if
  end if

 case('cc-pvdz-pp','cc-pvtz-pp','cc-pvqz-pp','aug-cc-pvdz-pp','aug-cc-pvtz-pp',&
      'aug-cc-pvqz-pp')
  i = INDEX(basis, '-pp')
  basis1 = basis(1:i-1)
  if(TRIM(RIJK_bas) == 'NONE') RIJK_bas = TRIM(basis1)//'/JK'
  if(dyn .and. TRIM(RIC_bas)=='NONE') RIC_bas = TRIM(basis)//'/C'

 case('cc-pwcvdz-pp','cc-pwcvtz-pp','cc-pwcvqz-pp','aug-cc-pwcvdz-pp',&
      'aug-cc-pwcvtz-pp','aug-cc-pwcvqz-pp')
  i = INDEX(basis, 'wc'); j = INDEX(basis, '-pp')
  basis1(1:i-1) = basis(1:i-1)
  basis1(i:) = basis(i+2:j-1)
  if(TRIM(RIJK_bas) == 'NONE') RIJK_bas = TRIM(basis1)//'/JK'
  if(dyn .and. TRIM(RIC_bas)=='NONE') RIC_bas = TRIM(basis)//'/C'

 case('cc-pvdz-f12','cc-pvtz-f12','cc-pvqz-f12')
  if(TRIM(RIJK_bas) == 'NONE') then
   i = INDEX(basis, '-f12', back=.true.)
   RIJK_bas = basis(1:i-1)//'/JK'
  end if
  if(dyn .and. TRIM(RIC_bas)=='NONE') then
   select case(basis)
   case('cc-pvdz-f12')
    RIC_bas = 'cc-pVDZ-F12-MP2Fit'
   case('cc-pvtz-f12')
    RIC_bas = 'cc-pVTZ-F12-MP2Fit'
   case('cc-pvqz-f12')
    RIC_bas = 'cc-pVQZ-F12-MP2Fit'
   end select
  end if
  if(F12) F12_cabs = TRIM(basis)//'-CABS'

 case('cc-pcvdz-f12','cc-pcvtz-f12','cc-pcvqz-f12')
  i = INDEX(basis, 'c'); j = INDEX(basis, '-f12')
  basis1(1:i-1) = basis(1:i-1)
  basis1(i:) = basis(i+1:j-1)
  if(TRIM(RIJK_bas) == 'NONE') RIJK_bas = TRIM(basis1)//'/JK'
  if(dyn .and. TRIM(RIC_bas)=='NONE') then
   select case(basis)
   case('cc-pcvdz-f12')
    RIC_bas = 'cc-pCVDZ-F12-MP2Fit'
   case('cc-pcvtz-f12')
    RIC_bas = 'cc-pCVTZ-F12-MP2Fit'
   case('cc-pcvqz-f12')
    RIC_bas = 'cc-pCVQZ-F12-MP2Fit'
   end select
  end if
  if(F12) F12_cabs = TRIM(basis1)//'-F12-CABS'

 case('cc-pvdz-pp-f12','cc-pvtz-pp-f12','cc-pvqz-pp-f12')
  i = INDEX(basis, '-pp')
  j = INDEX(basis, '-f12')
  if(TRIM(RIJK_bas) == 'NONE') RIJK_bas = basis(1:i-1)//'/JK'
  if(dyn .and. TRIM(RIC_bas)=='NONE') RIC_bas = basis(1:j-1)//'/C'
  if(F12) F12_cabs = basis(1:i-1)//'-F12-CABS'

 case('STO-3G','STO-6G','3-21G','6-31G','6-31G(d)','6-31G*','6-31G(d,p)','6-31G**')
  write(6,'(/,A)') REPEAT('-',79)
  write(6,'(A)') 'Warning in subroutine determine_auxbas: are you sure that you&
                 & want to use'
  write(6,'(A)') 'Pople-type basis set? This is less recommended. But automr wi&
                 &ll continue.'
  write(6,'(A)') REPEAT('-',79)
  RIJK_bas = 'def2/JK'
  if(dyn .and. TRIM(RIC_bas)=='NONE') RIC_bas = 'def2-SVP/C'

 case('6-311G','6-311G(d)','6-311G*','6-311G(d,p)','6-311G**')
  write(6,'(/,A)') REPEAT('-',79)
  write(6,'(A)') 'Warning in subroutine determine_auxbas: are you sure that you&
                 & want to use'
  write(6,'(A)') 'Pople-type basis set? This is less recommended. But automr wi&
                 &ll continue.'
  write(6,'(A)') REPEAT('-',79)
  RIJK_bas = 'def2/JK'
  if(dyn .and. TRIM(RIC_bas)=='NONE') RIC_bas = 'def2-TZVP/C'

 case default
  if(TRIM(RIJK_bas) == 'NONE') then
   write(6,'(/,A)') REPEAT('-',79)
   write(6,'(A)') 'Warning from subroutine determine_auxbas: RI-JK auxbas canno&
                  &t be determined'
   write(6,'(A)') 'from the orbital basis set. RIJK_bas=def2/JK will be used.'
   write(6,'(A)') REPEAT('-',79)
   RIJK_bas = 'def2/JK'
  end if
  if(dyn .and. TRIM(RIC_bas)=='NONE') then
   write(6,'(/,A)') REPEAT('-',79)
   write(6,'(A)') 'Warning from subroutine determine_auxbas: RIC auxbas cannot &
                  &be determined'
   write(6,'(A)') 'from the orbital basis set. RIC_bas=def2-TZVP/C will be used.'
   write(6,'(A)') REPEAT('-',79)
   RIC_bas = 'def2-TZVP/C'
  end if
 end select

 if(F12 .and. TRIM(F12_cabs)=='NONE') then
  write(6,'(/,A)') 'ERROR in subroutine determine_auxbas: near-complete auxilia&
                   &ry basis set for'
  write(6,'(A)') "F12 calculations '-CABS' cannot be automatically determined. &
                 &You should provide"
  write(6,'(A)') 'them in mokit{}.'
  stop
 end if
end subroutine determine_auxbas

! convert the name of auxiliary basis set from ORCA format to other format
subroutine auxbas_convert(inbas, outbas, itype)
 implicit none
 integer :: i
 integer, intent(in) :: itype ! 1/2 for PySCF/Molpro
 character(len=21) :: inbas1
 character(len=21), intent(in) :: inbas
 character(len=21), intent(out) :: outbas

 outbas = ' '
 inbas1 = inbas
 call lower(inbas1)

 select case(inbas1)
 case('cc-pvdz/jk','cc-pvtz/jk','cc-pvqz/jk','cc-pv5z/jk','aug-cc-pvdz/jk',&
      'aug-cc-pvtz/jk','aug-cc-pvqz/jk','aug-cc-pv5z/jk','mar-cc-pv5z/jk',&
      'apr-cc-pvqz/jk','apr-cc-pv5z/jk','may-cc-pvtz/jk','may-cc-pvqz/jk',&
      'may-cc-pv5z/jk','jun-cc-pvdz/jk','jun-cc-pvtz/jk','jun-cc-pvqz/jk',&
      'jun-cc-pv5z/jk')
  select case(itype)
  case(1) ! PySCF
   i = INDEX(inbas1, '/jk')
   outbas = inbas(1:i-1)//'-jkfit'
  case(2) ! Molpro
   outbas = TRIM(inbas)//'fit'
  case default
   write(6,'(/,A)') 'ERROR in subroutine auxbas_convert: invalid itype.'
   write(6,'(A,I0)') 'inbas='//TRIM(inbas)//', itype=', itype
   stop
  end select
 case('def2/jk')
  select case(itype)
  case(1) ! PySCF
   outbas = 'def2-universal-jkfit'
  case(2) ! Molpro
   outbas = 'qzvpp/jkfit'
   ! In the future(maybe Molpro 2021), we better use 'univjkfit/jkfit'.
   ! But so far there is a tiny bug for definition of def2-universal-JKFIT in
   ! Molpro, thus we have to use 'qzvpp/jkfit'. And these two are identical in fact.
  case default
   write(6,'(/,A)') 'ERROR in subroutine auxbas_convert: invalid itype.'
   write(6,'(A,I0)') 'inbas='//TRIM(inbas)//', itype=', itype
   stop
  end select
 case default
  write(6,'(/,A)') 'ERROR in subroutine auxbas_convert: inbas out of range.'
  write(6,'(A)') 'inbas1='//TRIM(inbas1)
  stop
 end select
end subroutine auxbas_convert

! calculate the Coulomb interaction energy of point charges
subroutine calc_Coulomb_energy_of_charges(n, charge, e)
 use phys_cons, only: Bohr_const
 implicit none
 integer :: i, j
 integer, intent(in) :: n
 real(kind=8) :: q1, q2, rtmp1(3), rtmp2(3)
 real(kind=8), intent(in) :: charge(4,n)
 ! charge(1:3,i) is the Cartesian coordinates of the i-th point charge
 ! charge(4,i) is the electronic charge of the i-th point charge
 real(kind=8), intent(out) :: e
 real(kind=8), parameter :: zero1 = 1d-2, zero2 = 1d-3
 real(kind=8), allocatable :: r(:,:)

 e = 0d0
 allocate(r(n,n))

 do i = 1, n-1, 1
 rtmp1 = charge(1:3,i)

  do j = i+1, n, 1
   rtmp2 = (rtmp1 - charge(1:3,j))/Bohr_const

   r(j,i) = DSQRT(DOT_PRODUCT(rtmp2, rtmp2))

   if(r(j,i) < zero2) then
    write(6,'(A)') 'ERROR in subroutine calc_Coulomb_energy_of_charges:'
    write(6,'(2(A,I0))') 'There exists two point charges too close: i=',i,',j=',j
    write(6,'(A,F15.8)') 'r(j,i) = ', r(j,i)
    stop
   else if(r(j,i) < zero1) then
    write(6,'(A)') 'Warning in subroutine calc_Coulomb_energy_of_charges:'
    write(6,'(2(A,I0))') 'There exists two point charges very close: i=',i,',j=',j
    write(6,'(A,F15.8)') 'r(j,i) = ', r(j,i)
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
end subroutine write_charge_into_chg

! calculate nuclear-point_charge interaction energy
subroutine calc_nuc_pt_e(nbgchg, bgcharge, natom, nuc, coor, nuc_pt_e)
 use phys_cons, only: Bohr_const
 implicit none
 integer :: i, j
 integer, intent(in) :: nbgchg, natom
 integer, intent(in) :: nuc(natom)
 real(kind=8) :: rtmp1(3), rtmp2(3), dis, pt_e
 real(kind=8), intent(in) :: bgcharge(4,nbgchg), coor(3,natom)
 real(kind=8), intent(out) :: nuc_pt_e

 nuc_pt_e = 0d0

 do i = 1, nbgchg, 1
  pt_e = bgcharge(4,i)
  rtmp1 = bgcharge(1:3,i)

  do j = 1, natom, 1
   rtmp2 = (rtmp1 - coor(:,j))/Bohr_const
   dis = DSQRT(DOT_PRODUCT(rtmp2, rtmp2))

   nuc_pt_e = nuc_pt_e + pt_e*DBLE(nuc(j))/dis
  end do ! for j
 end do ! for i

end subroutine calc_nuc_pt_e

! check whether a given binary file exists
subroutine check_exe_exist(path)
 implicit none
 character(len=240), intent(in) :: path
 logical :: alive

 inquire(file=TRIM(path),exist=alive)
 if(.not. alive) then
  write(6,'(A)') 'ERROR in subroutine check_exe_exist: the given binary file d&
                 &oes not exist.'
  write(6,'(A)') 'path='//TRIM(path)
  stop
 end if
end subroutine check_exe_exist

subroutine get_molpro_path(molpro_path)
 implicit none
 integer :: i, fid, system
 character(len=240), intent(out) :: molpro_path

 molpro_path = ' '
 i = SYSTEM("which molpro >mokit.molpro 2>&1")
 if(i /= 0) then
  molpro_path = 'mokit.molpro'
  call delete_file(TRIM(molpro_path))
  molpro_path = 'NOT FOUND'
  return
 end if

 open(newunit=fid,file='mokit.molpro',status='old',position='rewind')
 read(fid,'(A)',iostat=i) molpro_path
 close(fid,status='delete')

 if(i /= 0) then
  molpro_path = 'NOT FOUND'
 else
  if(LEN_TRIM(molpro_path) == 0) then
   molpro_path = 'NOT FOUND'
  else if(index(molpro_path,'no molpro') > 0) then
   molpro_path = 'NOT FOUND'
  end if
 end if
end subroutine get_molpro_path

subroutine get_psi4_path(psi4_path)
 implicit none
 integer :: i, fid, system
 character(len=240), intent(out) :: psi4_path

 psi4_path = ' '
 call getenv('PSI4', psi4_path)
 if(LEN_TRIM(psi4_path) > 0) return
 ! $PSI4 has higher priority than 'which psi4'

 i = SYSTEM("which psi4 >mokit.psi4 2>&1")
 open(newunit=fid,file='mokit.psi4',status='old',position='rewind')

 read(fid,'(A)',iostat=i) psi4_path
 close(fid,status='delete')

 if(i /= 0) then
  psi4_path = 'NOT FOUND'
 else
  if(LEN_TRIM(psi4_path) == 0) then
   psi4_path = 'NOT FOUND'
  else if(index(psi4_path,'no psi4') > 0) then
   psi4_path = 'NOT FOUND'
  end if
 end if
end subroutine get_psi4_path

subroutine get_dalton_path(dalton_path)
 implicit none
 integer :: i, fid, system
 character(len=240), intent(out) :: dalton_path

 dalton_path = ' '
 i = SYSTEM("which dalton >mokit.dalton 2>&1")

 open(newunit=fid,file='mokit.dalton',status='old',position='rewind')
 read(fid,'(A)',iostat=i) dalton_path
 close(fid,status='delete')

 if(i /= 0) then
  dalton_path = 'NOT FOUND'
 else
  if(LEN_TRIM(dalton_path) == 0) then
   dalton_path = 'NOT FOUND'
  else if(index(dalton_path,'no dalton') > 0) then
   dalton_path = 'NOT FOUND'
  end if
 end if
end subroutine get_dalton_path

! set memory and nproc in the module mr_keyword
subroutine set_mem_and_np_in_mr_keyword(mem_in, np_in)
 use mr_keyword, only: mem, nproc
 integer, intent(in) :: mem_in, np_in

 mem = mem_in
 nproc = np_in
end subroutine set_mem_and_np_in_mr_keyword

! check the sanity of the user-provided .fch(k) file (when readrhf/readuhf/readno
!  is specified)
subroutine check_sanity_of_provided_fch(DKH2, X2C, hf_fch)
 implicit none
 integer :: i
 character(len=240) :: buf
 character(len=240), intent(inout) :: hf_fch
 logical, intent(in) :: DKH2, X2C

 call require_file_exist(hf_fch)
 call read_basis_name_from_fch(hf_fch, buf)
 write(6,'(A79)') REPEAT('-',79)
 write(6,'(A)') 'Note: read user-specified .fch file. The basis set and ECP(i&
                &f any) will be im-'
 write(6,'(A)') 'ported from .fch file. The basis set you wrote in Route Sect&
                &ion(#p ...) will not'
 write(6,'(A)') "be used. The basis set name copied from specified .fch file &
                &is: '"//TRIM(buf)//"'"
 write(6,'(A79)') REPEAT('-',79)

 i = INDEX(hf_fch, '.fchk', back=.true.)
 if(i /= 0) then
  buf = hf_fch(1:i-1)//'.fch'
  call sys_copy_file(hf_fch, buf, .false.)
  hf_fch = buf
 end if

 if(DKH2) then
  call add_DKH2_into_fch(hf_fch)
 else if(X2C) then
  call add_X2C_into_fch(hf_fch)
 else
  call find_irel_in_fch(hf_fch, i)
  if(i == -3) then
   write(6,'(/,A)') REPEAT('-',55)
   write(6,'(A)') "Warning in subroutine parse_keyword: 'X2C' keyword detecte&
                  &d in file"
   write(6,'(A)') TRIM(hf_fch)//". But no 'X2C' keyword found in mokit{}. If &
                 &you do"
   write(6,'(A)') 'not want to perform X2C computations, please kill this job&
                  & immediately and'
   write(6,'(A)') "delete 'X2C' in .fch."
   write(6,'(A)') REPEAT('-',55)
  else if(i > 0) then
   write(6,'(/,A)') REPEAT('-',55)
   write(6,'(A)') 'Warning in subroutine parse_keyword: DKH related keywords &
                  &detected in file'
   write(6,'(A)') TRIM(hf_fch)//". But no 'DKH2' keyword found in mokit{}. If&
                  & you do"
   write(6,'(A)') 'not want to perform DKH2 computations, please kill this jo&
                  &b immediately'
   write(6,'(A)') 'and delete DKH related keywords in .fch.'
   write(6,'(A)') REPEAT('-',55)
  end if
 end if
end subroutine check_sanity_of_provided_fch

