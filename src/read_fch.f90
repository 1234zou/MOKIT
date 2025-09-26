! moved by jxzou at 20200322 from fch2mkl.f90 to this new file
! updated by jxzou at 20210215: add detection of 'Number of primitive s'

! The 'Shell types' array in Gaussian .fch file:
!
!   Spherical     |     Cartesian
! -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5
!  H  G  F  D  L  S  P  D  F  G  H
!
! 'L' is 'SP' in Pople-type basis sets

module fch_content
 use phys_cons, only: Bohr_const
 implicit none
 logical :: is_uhf = .false. ! is UHF or not
 integer :: irel = -1        ! -3/-2/-1/0/2/4 for sfX2C/RESC/None/DKH0/DKH2/DKH4
 ! relativistic Hamiltonian. Note that irel=0 means DKH0, not None
 integer :: nbf, nif         ! number of basis functions/MOs
 integer :: na, nb, nopen    ! number of alpha/beta/open shell electrons
 integer :: ncontr, nprim    ! Number of contracted/primitive shells
 integer :: charge = 0       ! total charge
 integer :: mult  = 1        ! spin multiplicity
 integer :: natom = 0        ! number of atoms
 integer :: LenNCZ = 0       ! ECP-LenNCZ
 integer :: nmode = 0        ! Number of Normal Modes
 integer, parameter :: period_nelem = 118    ! 118 elements, H-Og
 integer, parameter :: shltyp2nbf(-5:5) = (/11,9,7,5,4,1,3,6,10,15,21/)
 integer, allocatable :: ielem(:)            ! elements, 6 for 'C', etc
 integer, allocatable :: iatom_type(:)       ! Int Atom Types
 integer, allocatable :: shell_type(:)       ! Shell types
 integer, allocatable :: prim_per_shell(:)   ! Number of primitives per shell
 integer, allocatable :: shell2atom_map(:)   ! Shell to atom map
 integer, allocatable :: KFirst(:,:)         ! ECP-KFirst, size (natom,10)
 integer, allocatable :: KLast(:,:)          ! ECP-KLast, size (natom,10)
 integer, allocatable :: Lmax(:), LPSkip(:)  ! ECP-LMax, ECP-LPSkip, size natom
 integer, allocatable :: NLP(:)              ! ECP-NLP, size LenNCZ
 real(kind=8) :: virial = 2d0 ! Virial Ratio
 real(kind=8) :: tot_e = 0d0 ! Total Energy
 real(kind=8), allocatable :: coor(:,:)         ! Cartesian coordinates
 real(kind=8), allocatable :: prim_exp(:)       ! Primitive exponents
 real(kind=8), allocatable :: contr_coeff(:)    ! Contraction coefficients
 real(kind=8), allocatable :: contr_coeff_sp(:) ! (S=P) Contraction coefficients
 real(kind=8), allocatable :: eigen_e_a(:)      ! Alpha energy levels
 real(kind=8), allocatable :: eigen_e_b(:)      ! Beta energy levels
 real(kind=8), allocatable :: alpha_coeff(:,:)  ! Alpha MOs
 real(kind=8), allocatable :: beta_coeff(:,:)   ! Beta MOs
 real(kind=8), allocatable :: RNFroz(:)         ! ECP-RNFroz(core e), size natom
 real(kind=8), allocatable :: CLP(:), ZLP(:)    ! ECP-CLP1, ECP-ZLP, size LenNCZ
 real(kind=8), allocatable :: CLP2(:)           ! ECP-CLP2, size LenNCZ, for GHF
 real(kind=8), allocatable :: tot_dm(:,:)       ! Total SCF Density
 real(kind=8), allocatable :: spin_dm(:,:)      ! Spin SCF Density
 real(kind=8), allocatable :: mull_char(:)      ! Mulliken Charges
 real(kind=8), allocatable :: force_const(:,:)  ! Cartesian Force Constants
 real(kind=8), allocatable :: vibe2(:)          ! Vib-E2, size 14*nmode
 real(kind=8), allocatable :: norm_mode(:,:)    ! Vib-Modes
 real(kind=8), allocatable :: rnuc(:) ! nuclear charges with real(kind=8) type
 character(len=2), allocatable :: elem(:)       ! elements ('H ', 'C ', etc)

 ! Frozen core orbitals used in RHF_proj and dynamic correlation computations.
 ! This table is copied from the figure in ORCA 5.0.1 manual 9.11 Frozen Core
 ! Options.
 ! Note: frozen electrons are two times of the frozen core orbitals
 integer, parameter :: core_orb(0:period_nelem) = (/0,&           ! Bq
     0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  5,  5,  5, & ! H~P
     5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5, & ! S~Zn
     9,  9,  9,  9,  9,  9,  9,  9, 14, 14, 14, 14, 14, 14, 14, & ! Ga~Rh
    14, 14, 14, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, & ! Pd~Nd
    18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 23, 23, 23, 23, 23, & ! Pm~Re
    23, 23, 23, 23, 23, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, & ! Os~Th
    34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 50, 50, & ! Pa~Db
    50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50/)          ! Sg~Cn

 ! Taken from
 ! https://www.ciaaw.org/atomic-weights.htm#m
 ! https://lynceans.org/tag/extended-periodic-table
 ! https://www.nist.gov/sites/default/files/images/2019/07/26/nist_periodictable_july2019_final_front.jpg
! real(kind=8), parameter :: ram(0:period_nelem) = (/0d0, 1.008d0,4.003d0,     &! Bq, H, He
!  6.938d0, 9.012d0, 10.806d0, 12.01d0, 14.006d0, 16d0, 18.998d0, 20.18d0,     &! Li~Ne
!  22.99d0, 24.304d0, 26.982d0, 28.084d0, 30.974d0, 32.059d0, 35.446d0,        &! Na~Cl
!  39.792d0, 39.098d0, 40.078d0, 44.956d0, 47.867d0, 50.942d0, 51.996d0,       &! Ar~Cr
!  54.938d0, 55.845d0, 58.933d0, 58.693d0, 63.546d0, 65.382d0, 69.723d0,       &! Mn~Ga
!  72.631d0, 74.922d0, 78.972d0, 79.901d0, 83.798d0, 85.468d0, 87.621d0,       &! Ge~Sr
!  88.906d0, 91.224d0, 92.906d0, 95.951d0, 98.907d0, 101.072d0, 102.905d0,     &! Y~Rh
!  106.421d0, 107.868d0, 112.414d0, 114.818d0, 118.711d0, 121.76d0, 127.603d0, &! Pd~Te
!  126.904d0, 131.294d0, 132.905d0, 137.328d0, 138.905d0, 140.116d0, 140.908d0,&! I~Pr
!  144.242d0, 144.913d0, 150.362d0, 151.964d0, 157.253d0, 158.925d0, 162.5d0,  &! Nd~Dy
!  164.93d0, 167.259d0, 168.934d0, 173.045d0, 174.967d0, 178.487d0, 180.948d0, &! Ho~Ta
!  183.841d0, 186.207d0, 190.233d0, 192.217d0, 195.085d0, 196.967d0, 200.592d0,&! W~Hg
!  204.382d0, 206.14d0, 208.98d0, 208.982d0, 209.987d0, 222.018d0, 223.02d0,   &! Tl~Fr
!  226.025d0, 227.028d0, 232.038d0, 231.036d0, 238.029d0, 237.048d0, 244.064d0,&! Ra~Pu
!  243.061d0, 247.07d0, 247.07d0, 251.08d0, 254d0, 257.095d0, 258.1d0,         &! Am~Md
!  259.101d0, 266d0, 267d0, 268d0, 269d0, 270d0, 269d0, 278d0, 281d0, 282d0,   &! No~Rg
!  285d0, 286d0, 289d0, 289d0, 293d0, 294d0, 294d0/)                            ! Cn~Og

 ! H~Lr taken from the G16 output file.
 real(kind=8), parameter :: ram(0:period_nelem)=(/0d0,1.00782504d0,4.0026032d0,&! Bq,H,He
   7.0160045d0,   9.0121825d0,  11.0093053d0,  12.0000000d0,  14.0030740d0,  &! Li~N
  15.9949146d0,  18.9984032d0,  19.9924391d0,  22.9897697d0,  23.9850450d0,  &! O~Mg
  26.9815413d0,  27.9769284d0,  30.9737634d0,  31.9720718d0,  34.9688527d0,  &! Al~Cl
  39.9623831d0,  38.9637079d0,  39.9625907d0,  44.9559136d0,  47.9479467d0,  &! Ar~Ti
  50.9439625d0,  51.9405097d0,  54.9380463d0,  55.9349393d0,  58.9331978d0,  &! V~Co
  57.9353471d0,  62.9295992d0,  63.9291454d0,  68.9255809d0,  73.9211788d0,  &! Ni~Ge
  74.9215955d0,  79.9165205d0,  78.9183361d0,  83.9115064d0,  84.9117000d0,  &! As~Rb
  87.9056000d0,  88.9054000d0,  89.9043000d0,  92.9060000d0,  97.9055000d0,  &! Sr~Mo
  98.9063000d0, 101.9037000d0, 102.9048000d0, 105.9032000d0, 106.9050900d0,  &! Tc~Ag
 113.9036000d0, 114.9041000d0, 117.9018000d0, 120.9038000d0, 129.9067000d0,  &! Cd~Te
 126.9004000d0, 131.9042000d0, 132.9054290d0, 137.9050000d0, 138.9061000d0,  &! I~La
 139.9053000d0, 140.9074000d0, 141.9075000d0, 144.9127000d0, 151.9195000d0,  &! Ce~Sm
 152.9209000d0, 157.9241000d0, 158.9250000d0, 163.9288000d0, 164.9303000d0,  &! Eu~Ho
 165.9304000d0, 168.9344000d0, 173.9390000d0, 174.9409000d0, 179.9468000d0,  &! Er~Hf
 180.9480000d0, 183.9510000d0, 186.9560000d0, 189.9586000d0, 192.9633000d0,  &! Ta~Ir
 194.9648000d0, 196.9666000d0, 201.9706000d0, 204.9745000d0, 207.9766000d0,  &! Pt~Pb
 208.9804000d0, 208.9825000d0, 210.9875000d0, 222.0175000d0, 223.0198000d0,  &! Bi~Fr
 226.0254000d0, 227.0278000d0, 232.0382000d0, 231.0359000d0, 238.0508000d0,  &! Ra~U
 237.0480000d0, 242.0587000d0, 243.0614000d0, 246.0674000d0, 247.0702000d0,  &! Np~Bk
 249.0748000d0, 252.0829000d0, 252.0827000d0, 255.0906000d0, 259.1010000d0,  &! Cf~No
 262.1097d0, 267d0, 268d0, 269d0, 270d0, 269d0, 278d0, 281d0, 282d0, 285d0,  &! Lr~Cn
 286d0, 289d0, 289d0, 293d0, 294d0, 294d0/)                                   ! Nh~Og

 character(len=2), parameter :: period_elem(0:period_nelem) = (/'Bq',&
   'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
   'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca', &
   'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', &
   'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr', &
   'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', &
   'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', &
   'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', &
   'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', &
   'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', &
   'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', &
   'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', &
   'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'/)

contains

! read the number of electrons from a given .fch file
subroutine read_ne_from_fch(fchname, ne)
 implicit none
 integer :: i, fid
 integer, intent(out) :: ne
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname

 ne = 0
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:14) == 'Number of elec') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_ne_from_fch: no 'Number of elec' f&
                   &ound in file "//TRIM(fchname)
  stop
 end if

 read(buf(50:),*) ne
end subroutine read_ne_from_fch

! check whether GHF-type MOs are hold in a given .fch(k) file
subroutine check_ghf_in_fch(fchname, ghf)
 implicit none
 integer :: i, nif, norb, fid
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
 logical, intent(out) :: ghf

 nif = 0; norb = 0; ghf = .false.
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:13) == 'Number of ind') exit
 end do ! for while

 if(i /= 0) then
  close(fid)
  write(6,'(/,A)') "ERROR in subroutine check_ghf_in_fch: no 'Number of ind' fo&
                   &und in file "//TRIM(fchname)
  stop
 end if

 read(buf(50:),*) nif

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i < 0) exit ! end-of-file
  if(i > 0) then
   write(6,'(/,A)') 'ERROR in subroutine check_ghf_in_fch: failed to read file &
                   &'//TRIM(fchname)
   close(fid)
   stop
  end if

  if(buf(1:7) == 'Alpha O') then
   read(buf(50:),*) norb
   exit
  end if
 end do ! for while

 close(fid)
 if(norb == 2*nif) ghf = .true.
end subroutine check_ghf_in_fch

! read geometry, basis sets, ECP(if any) and MOs from .fch file (Gaussian)
subroutine read_fch(fchname, uhf)
 implicit none
 integer :: i, j, fid
 integer :: ncoeff ! = nbf*nif
 integer, allocatable :: itmp(:)
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
 real(kind=8) :: sum_res
 real(kind=8), allocatable :: coor0(:) ! Cartesian coordinates
 real(kind=8), allocatable :: coeff(:) ! MOs, nbf*nif
 logical, intent(in) :: uhf
 logical :: ecp, has_rnuc, int_atom_typ

 ecp = .false.; has_rnuc = .false.; buf = ' '

 call find_specified_suffix(fchname, '.fch', i)
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 ! find variables: charge, mult, na, nb, nbf, nif
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:6) == 'Charge') exit
 end do
 BACKSPACE(fid)
 read(fid,'(A49,2X,I10)') buf, charge
 read(fid,'(A49,2X,I10)') buf, mult
 read(fid,'(A)') buf   ! skip 'Number of electrons'
 read(fid,'(A49,2X,I10)') buf, na
 read(fid,'(A49,2X,I10)') buf, nb
 read(fid,'(A49,2X,I10)') buf, nbf
 read(fid,'(A49,2X,I10)') buf, nif
 ncoeff = nbf*nif

 nopen = na - nb
 if(nopen < 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_fch: The number of alpha electrons&
                   & is less than'
  write(6,'(A)') 'that of beta electrons!'
  write(6,'(A)') 'fchname='//TRIM(fchname)
  write(6,'(2(A,I0))') 'na=', na, ', nb=', nb
  close(fid)
  stop
 end if

 ! find and read natom, elem
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:14) == 'Atomic numbers') exit
 end do ! for while

 BACKSPACE(fid)
 read(fid,'(A49,2X,I10)') buf, natom
 allocate(ielem(natom), source=0)
 read(fid,'(6(1X,I11))') (ielem(i),i=1,natom)
 allocate(elem(natom))
 elem = ' '
 forall(i = 1:natom) elem(i) = nuc2elem(ielem(i))

 ! find 'Nuclear charges'
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:12) == 'Nuclear char') then
   has_rnuc = .true.
   exit
  end if
  if(buf(1:12) == 'Current cart') then
   BACKSPACE(fid)
   exit
  end if
 end do ! for while

 if(has_rnuc) then
  allocate(rnuc(natom))
  read(fid,'(5(1X,ES15.8))') (rnuc(i),i=1,natom)
  allocate(itmp(natom), source=NINT(rnuc))
  ! if it is ab inito all-electron basis set, this array is useless
  if(SUM(IABS(itmp-ielem)) == 0) deallocate(rnuc)
  deallocate(itmp)
  ! If the array rnuc still exists, the case should be either ab inito calculat-
  ! ion with ECP/PP, or semi-empirical calculation with only valence electrons.
 end if

 ! find and read coordinates
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:12) == 'Current cart') exit
 end do ! for while

 allocate(coor0(3*natom), coor(3,natom))
 read(fid,'(5(1X,ES15.8))') (coor0(i),i=1,3*natom)
 coor = RESHAPE(coor0,(/3,natom/))
 deallocate(coor0)
 coor = coor*Bohr_const ! convert Bohr to Anstrom

 ! Find and read 'Int Atom Types'. The fch(k) file generated by some programs
 ! does not have this section, and all atoms would be viewed as real atoms in
 ! such case.
 int_atom_typ = .false.
 do while(.true.)
  read(fid,'(A)') buf
  select case(buf(1:13))
  case('Int Atom Type')
   int_atom_typ = .true.
   exit
  case('MM charges   ','Real atomic w','Atom fragment')
   exit
  case('Number of con')
   BACKSPACE(fid)
   exit
  end select
 end do ! for while

 allocate(iatom_type(natom), source=0)
 if(int_atom_typ) read(fid,'(6(1X,I11))') iatom_type
 ! 0 for real atoms, 1000 for ghost atoms

 ! find and read variables ncontr, nprim
 i = 0
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:18) == 'Number of contract') then
   i = 1
   exit
  end if
  if(buf(1:21) == 'Number of primitive s') then
   i = 2
   exit
  end if
 end do ! for while

 BACKSPACE(fid)
 select case(i)
 case(1)
  read(fid,'(A49,2X,I10)') buf, ncontr
  read(fid,'(A49,2X,I10)') buf, nprim
 case(2)
  read(fid,'(A49,2X,I10)') buf, nprim
  read(fid,'(A49,2X,I10)') buf, ncontr
 case default
  write(6,'(/,A,I0)') 'ERROR in subroutine read_fch: invalid i=', i
  write(6,'(A)') 'Unsupported file format, or this file is incomplete.'
  write(6,'(A)') 'Filename='//TRIM(fchname)
  close(fid)
  stop
 end select

 ! find and read Shell types, Number of primitives per shell, Shell to atom map,
 !  Primitive exponents, Contraction coefficients
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:14) == 'Shell types') exit
 end do ! for while
 allocate(shell_type(ncontr), source=0)
 read(fid,'(6(1X,I11))') (shell_type(i), i=1,ncontr)

 read(fid,'(A)') buf
 allocate(prim_per_shell(ncontr), source=0)
 read(fid,'(6(1X,I11))') (prim_per_shell(i), i=1,ncontr)

 read(fid,'(A)') buf
 allocate(shell2atom_map(ncontr), source=0)
 read(fid,'(6(1X,I11))') (shell2atom_map(i), i=1,ncontr)

 read(fid,'(A)') buf
 allocate(prim_exp(nprim), source=0d0)
 read(fid,'(5(1X,ES15.8))') (prim_exp(i),i=1,nprim)

 read(fid,'(A)') buf
 allocate(contr_coeff(nprim), source=0d0)
 read(fid,'(5(1X,ES15.8))') (contr_coeff(i),i=1,nprim)

 read(fid,'(A)') buf
 if(buf(1:6) == 'P(S=P)') then
  allocate(contr_coeff_sp(nprim), source=0d0)
  read(fid,'(5(1X,ES15.8))') (contr_coeff_sp(i),i=1,nprim)
 end if

 ! check if there is any ECP/PP data
 ! Note that def2TZVP for the 4th row will cause LenNCZ=0
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:10) == 'ECP-LenNCZ') then
   read(buf(52:),*) LenNCZ
   if(LenNCZ > 0) then
    ecp = .true.
    BACKSPACE(fid)
   end if
   exit
  end if

  if(buf(1:12) == 'Virial Ratio') then
   BACKSPACE(fid)
   exit
  end if
  if(buf(1:13) == 'Alpha Orbital') then
   BACKSPACE(fid)
   exit
  end if
 end do ! for while

 ! if ecp is true, read ECP/PP data
 if(ecp) then
  read(fid,'(A49,2X,I10)') buf, LenNCZ
  allocate(KFirst(natom,10), source=0)
  allocate(KLast(natom,10), source=0)
  allocate(LPSkip(natom), source=0)
  allocate(RNFroz(natom), source=0d0)
  allocate(NLP(LenNCZ), source=0)
  allocate(CLP(LenNCZ), source=0d0)
  allocate(CLP2(LenNCZ), source=0d0)
  allocate(ZLP(LenNCZ), source=0d0)

  allocate(Lmax(10*natom), source=0)
  read(fid,'(A)') buf
  read(fid,'(6(1X,I11))') (Lmax(i), i=1,10*natom)
  KFirst = RESHAPE(Lmax, (/natom,10/))

  read(fid,'(A)') buf
  read(fid,'(6(1X,I11))') (Lmax(i), i=1,10*natom)
  KLast  = RESHAPE(Lmax, (/natom,10/))
  deallocate(Lmax)

  allocate(Lmax(natom), source=0)
  read(fid,'(A)') buf
  read(fid,'(6(1X,I11))') (Lmax(i), i=1,natom)
  read(fid,'(A)') buf
  read(fid,'(6(1X,I11))') (LPSkip(i), i=1,natom)
  read(fid,'(A)') buf
  read(fid,'(5(1X,ES15.8))') (RNFroz(i), i=1,natom)
  read(fid,'(A)') buf
  read(fid,'(6(1X,I11))') (NLP(i), i=1,LenNCZ)
  read(fid,'(A)') buf
  read(fid,'(5(1X,ES15.8))') (CLP(i), i=1,LenNCZ)
  read(fid,'(A)') buf
  read(fid,'(5(1X,ES15.8))') (CLP2(i), i=1,LenNCZ)
  read(fid,'(A)') buf
  read(fid,'(5(1X,ES15.8))') (ZLP(i), i=1,LenNCZ)
 end if

 ! read Virial Ratio and Total Energy
 do while(.true.)
  read(fid,'(A)') buf
  select case(buf(1:12))
  case('Virial Ratio')
   read(buf(45:),*) virial
  case('Total Energy')
   read(buf(45:),*) tot_e
  case('Alpha Orbita')
   exit
  case default ! do nothing
  end select
 end do ! for while

 ! read Alpha Orbital Energies
 allocate(eigen_e_a(nif), source=0d0)
 read(fid,'(5(1X,ES15.8))') (eigen_e_a(i),i=1,nif)

 ! if '-uhf' specified, read Beta Orbital Energies
 if(uhf) then
  read(fid,'(A)') buf
  if(buf(1:12) /= 'Beta Orbital') then
   write(6,'(/,A)') "ERROR in subroutine read_fch: no 'Beta Orbital' found in t&
                    &his .fch(k) file,"
   write(6,'(A)') "but you specify '-uhf'. fchname="//TRIM(fchname)
   close(fid)
   stop
  end if
  allocate(eigen_e_b(nif), source=0d0)
  read(fid,'(5(1X,ES15.8))') (eigen_e_b(i),i=1,nif)
 end if

 ! read Alpha MO
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:8) == 'Alpha MO') exit
 end do
 allocate(coeff(ncoeff), source=0d0)
 !read(fid,'(5(1X,ES15.8))') (coeff(i),i=1,ncoeff)
 read(fid,*) coeff
 allocate(alpha_coeff(nbf,nif), source=RESHAPE(coeff,[nbf,nif]))
 deallocate(coeff)

 ! if '-uhf' specified, read Beta MO
 if(uhf) then
  allocate(coeff(ncoeff), source=0d0)
  read(fid,'(A)') buf   ! skip the Beta MO line
  !read(fid,'(5(1X,ES15.8))') (coeff(i),i=1,ncoeff)
  read(fid,*) coeff
  allocate(beta_coeff(nbf,nif), source=RESHAPE(coeff,[nbf,nif]))
  deallocate(coeff)
 end if

 do while(.true.)   ! read Total SCF Density
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) then
   close(fid)
   return
  end if
  if(buf(1:9) == 'Total SCF') exit
 end do ! for while

 allocate(tot_dm(nbf,nbf))
 read(fid,'(5(1X,ES15.8))') ((tot_dm(j,i),j=1,i),i=1,nbf)
 forall(i=1:nbf,j=1:nbf,i>j) tot_dm(i,j) = tot_dm(j,i)

 if(uhf) then   ! read Spin SCF Density
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) then
    close(fid)
    return
   end if
   if(buf(1:9) == 'Spin SCF ') exit
  end do ! for while
  allocate(spin_dm(nbf,nbf))
  read(fid,'(5(1X,ES15.8))') ((spin_dm(j,i),j=1,i),i=1,nbf)
  forall(i=1:nbf,j=1:nbf,i>j) spin_dm(i,j) = spin_dm(j,i)
 end if

 close(fid)

 if(allocated(rnuc) .and. ecp) then
  sum_res = 0d0
  do i = 1, natom, 1
   if(iatom_type(i) == 1000) cycle ! skip ghost atoms
   sum_res = sum_res + DABS(DBLE(ielem(i))- RNFroz(i) - rnuc(i))
  end do ! for i
  if(sum_res < 1d-6) then
   deallocate(rnuc)
  else
   write(6,'(/,A)') 'ERROR in subroutine read_fch: unknown basis type.'
   write(6,'(A)') 'This seems not an all-electron basis set, not with ECP/PP, n&
                  &ot semi-empirical.'
   stop
  end if
 end if
end subroutine read_fch

! map a nuclear charge to an element (e.g. 6->'C')
pure function nuc2elem(i) result(s)
 implicit none
 integer, intent(in) :: i
 character(len=2) :: s

 s = period_elem(i)
end function nuc2elem

! convert atomic numbers to atomic symbols if the user uses atomic numbers
subroutine possible_nuc2elem(natom, elem)
 implicit none
 integer :: i, j
 integer, intent(in) :: natom
 character(len=2), intent(inout) :: elem(natom)

 do i = 1, natom, 1
  j = IACHAR(elem(i)(1:1))
  if(j>47 .and. j<58) then
   read(elem(i),*) j
   elem(i) = period_elem(j)
  end if
 end do ! for i
end subroutine possible_nuc2elem

! map an element to a nuclear charge (e.g. 'C'->6)
! This is a pure function, `write` or `stop` cannot be used here. It is the
! user's responsibility to make sure the input element is valid.
pure function elem2nuc(s) result(i)
 implicit none
 integer :: i
 character(len=2), intent(in) :: s

 do i = 0, period_nelem, 1
  if(period_elem(i) == s) return
 end do ! for i
end function elem2nuc

! deallocate arays in module fch_content
subroutine free_arrays_in_fch_content()
 implicit none

 deallocate(ielem, iatom_type, shell_type, prim_per_shell, shell2atom_map, &
            coor, prim_exp, contr_coeff, eigen_e_a, alpha_coeff)
 if(is_uhf) deallocate(eigen_e_b, beta_coeff)
 if(LenNCZ > 0) deallocate(KFirst, KLast, Lmax, LPSkip, NLP, RNFroz, CLP, CLP2,&
                           ZLP)
 if(allocated(rnuc)) deallocate(rnuc)
 if(allocated(elem)) deallocate(elem)
 if(allocated(contr_coeff_sp)) deallocate(contr_coeff_sp)
 if(allocated(tot_dm)) deallocate(tot_dm)
 if(allocated(spin_dm)) deallocate(spin_dm)
 if(allocated(mull_char)) deallocate(mull_char)
end subroutine free_arrays_in_fch_content

end module fch_content

! expansion coefficients matrices of spherical harmonic -> Cartesian functions
module parameter_sph2cart
 implicit none
 real(kind=8), parameter :: r1 =0.5d0, r2=0.5d0*DSQRT(3d0), r3=1.5d0/DSQRT(5d0)
 real(kind=8), parameter :: r4 =DSQRT(0.375d0), r5=DSQRT(3d0/40d0), r6=DSQRT(1.2d0)
 real(kind=8), parameter :: r7 =DSQRT(0.625d0), r8=3d0/DSQRT(8d0), r9=DSQRT(27d0/35d0)
 real(kind=8), parameter :: r10=3d0/8d0, r11=0.25d0*r9, r12=DSQRT(10d0/7d0)
 real(kind=8), parameter :: r13=DSQRT(9d0/56d0), r14=DSQRT(45d0/56d0), r15=DSQRT(27d0/28d0)
 real(kind=8), parameter :: r16=0.25*DSQRT(5d0), r17=DSQRT(9d0/7d0), r18=DSQRT(5d0/28d0)
 real(kind=8), parameter :: r19=0.125d0*DSQRT(35d0), r20=0.75d0*DSQRT(3d0), r21=2d0*r16
 real(kind=8), parameter :: r22=DSQRT(25d0/21d0), r23=0.625d0, r24=DSQRT(15d0/112d0)
 real(kind=8), parameter :: r25=DSQRT(5d0/3d0), r26=DSQRT(9d0/28d0), r27=0.125d0*r25
 real(kind=8), parameter :: r28=DSQRT(45d0/28d0), r29=DSQRT(5d0/112d0), r30=0.125d0*DSQRT(15d0)
 real(kind=8), parameter :: r31=DSQRT(35d0/48d0), r32=DSQRT(5d0/12d0), r33=DSQRT(1.5d0)
 real(kind=8), parameter :: r34=DSQRT(35d0/128d0), r35=DSQRT(5d0/6d0), r36=0.25d0*r35
 real(kind=8), parameter :: r37=DSQRT(175d0/128d0), r38=1.25d0*r33, r39=DSQRT(63d0/128d0)

! The transformation table (rd,rf,rg,rh) is originally copied from http://sobereva.com/97
 real(kind=8), parameter :: rd(6,5) = RESHAPE([-r1, -r1, 1d0, 0d0, 0d0, 0d0,&
  0d0, 0d0, 0d0, 0d0, 1d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 1d0, &
   r2, -r2, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 1d0, 0d0, 0d0], [6,5])
 real(kind=8), parameter :: &
  rf(10,7) = RESHAPE([0d0, 0d0, 1d0, 0d0, 0d0, -r3, 0d0, 0d0, -r3, 0d0,&
                      -r4, 0d0, 0d0, -r5, 0d0, 0d0,  r6, 0d0, 0d0, 0d0,&
                      0d0, -r4, 0d0, 0d0, -r5, 0d0, 0d0,  r6, 0d0, 0d0,&
                      0d0, 0d0, 0d0, 0d0, 0d0,  r2, 0d0, 0d0, -r2, 0d0,&
                      0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 1d0,&
                       r7, 0d0, 0d0, -r8, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,&
                      0d0, -r7, 0d0, 0d0,  r8, 0d0, 0d0, 0d0, 0d0, 0d0],[10,7])
 real(kind=8), parameter :: &
  rg(15,9) = RESHAPE([1d0, 0d0, -r9, 0d0, r10, 0d0, 0d0, 0d0, 0d0, -r9, 0d0, r11, 0d0, 0d0, r10,&
                      0d0, 0d0, 0d0, 0d0, 0d0, r12, 0d0,-r13, 0d0, 0d0, 0d0, 0d0,-r14, 0d0, 0d0,&
                      0d0, r12, 0d0,-r14, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r13, 0d0, 0d0, 0d0, 0d0,&
                      0d0, 0d0,-r15, 0d0, r16, 0d0, 0d0, 0d0, 0d0, r15, 0d0, 0d0, 0d0, 0d0,-r16,&
                      0d0, 0d0, 0d0, 0d0, 0d0, 0d0, r17, 0d0,-r18, 0d0, 0d0, 0d0, 0d0,-r18, 0d0,&
                      0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, -r8, 0d0, 0d0, 0d0, 0d0,  r7, 0d0, 0d0,&
                      0d0, 0d0, 0d0, -r7, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,  r8, 0d0, 0d0, 0d0, 0d0,&
                      0d0, 0d0, 0d0, 0d0, r19, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r20, 0d0, 0d0, r19,&
                      0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r21, 0d0, 0d0, 0d0, 0d0, r21, 0d0],[15,9])
 real(kind=8), parameter :: &
rh(21,11)=RESHAPE([1d0, 0d0,-r22, 0d0, r23, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r22, 0d0, r24, 0d0, 0d0, 0d0, 0d0, r23, 0d0, 0d0,&
                   0d0, 0d0, 0d0, 0d0, 0d0, 0d0, r25, 0d0,-r26, 0d0, r27, 0d0, 0d0, 0d0, 0d0,-r28, 0d0, r29, 0d0, 0d0, r30,&
                   0d0, r25, 0d0,-r28, 0d0, r30, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r26, 0d0, r29, 0d0, 0d0, 0d0, 0d0, r27, 0d0,&
                   0d0, 0d0,-r21, 0d0, r31, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, r21, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r31, 0d0, 0d0,&
                   0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, r25, 0d0,-r32, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r32, 0d0, 0d0, 0d0, 0d0,&
                   0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r33, 0d0, r34, 0d0, 0d0, 0d0, 0d0, r35, 0d0, r36, 0d0, 0d0,-r34,&
                   0d0, 0d0, 0d0,-r35, 0d0, r34, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, r33, 0d0,-r36, 0d0, 0d0, 0d0, 0d0,-r34, 0d0,&
                   0d0, 0d0, 0d0, 0d0, r19, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r20, 0d0, 0d0, 0d0, 0d0, r19, 0d0, 0d0,&
                   0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r21, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, r21, 0d0, 0d0, 0d0, 0d0,&
                   0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, r37, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r38, 0d0, 0d0, r39,&
                   0d0, 0d0, 0d0, 0d0, 0d0, r39, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r38, 0d0, 0d0, 0d0, 0d0, r37, 0d0],&
                  [21,11])
end module parameter_sph2cart

! check whether there exists DKH/X2C keywords in a specified .fch(k) file
subroutine find_irel_in_fch(fchname, irel)
 implicit none
 integer :: i, fid
 integer, intent(out) :: irel
!f2py intent(out) :: irel
! -3: X2C, i.e. sf-x2c1e | -2: RESC
! -1: no relativity      |  0: DKH 0th-order
!  2: DKH2               |  4: DKH4 with SO
 character(len=61) :: buf
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname
 character(len=610) :: longbuf
 logical :: alive(7)

 irel = -1   ! default: no relativity
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:5)=='Route' .or. buf(1:6)=='Charge') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine find_irel_in_fch: neither 'Route' nor '&
                   &Charge' is found"
  write(6,'(A)') 'in file '//TRIM(fchname)//'.'
  write(6,'(A)') 'This .fch(k) file maybe incomplete.'
  close(fid)
  stop
 end if

 if(buf(1:6) == 'Charge') then ! no 'Route' found
  close(fid)
  return
 end if

 longbuf = ' '
 do i = 1, 10
  read(fid,'(A)') buf
  if(buf(1:6) == 'Charge') exit
  longbuf = TRIM(longbuf)//TRIM(buf)
 end do ! for i
 close(fid)

 call upper(longbuf)

 alive = [(index(longbuf,'DKH0')>0), (index(longbuf,'DKHSO')>0), &
          (index(longbuf,'DKH')>0 .and. index(longbuf,'NODKH')==0), &
          (index(longbuf,'DKH2')>0), (index(longbuf,'DOUGLASKROLLHESS')>0), &
          (index(longbuf,'RESC')>0), (index(longbuf,'X2C')>0)]

 if(alive(1)) then      ! DKH0
  irel = 0
 else if(alive(2)) then ! DKH4, DKHSO
  irel = 4
 else if(alive(6)) then ! RESC
  irel = -2
 else if(alive(7)) then ! X2C
  irel = -3
 else if(alive(3) .or. alive(4) .or. alive(5)) then ! DKH2
  irel = 2
 else
  irel = -1             ! no relativity
 end if
end subroutine find_irel_in_fch

! set the type of relativistic Hamiltonian in a specified .fch(k) file
subroutine set_irel_in_fch(fchname, irel)
 implicit none
 integer :: i, j, fid, fid1, RENAME
 integer, intent(in) :: irel
!f2py intent(in) :: irel
 character(len=7) :: ham(-3:4) = [',X2C)  ',',RESC) ',')      ',',DKH0) ',&
                                  ',DKH1) ',',DKH2) ',',DKH3) ',',DKHSO)']
 character(len=240) :: buf, fchname1
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname
 character(len=960) :: route

 call find_specified_suffix(fchname, '.fch', i)
 fchname1 = fchname(1:i-1)//'.t'
 call read_route_from_fch(fchname, route)

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(fchname1),status='replace')

 if(LEN_TRIM(route) == 0) then ! no Route in .fch(k)
  do while(.true.)
   read(fid,'(A)') buf
   if(buf(1:6) == 'Charge') exit
   write(fid1,'(A)') TRIM(buf)
  end do ! for while

  if(irel == -1) then
   write(fid1,'(A,38X,A)') 'Route','C   N=           4'
   write(fid1,'(A)') '#p HF chkbasis nosymm int=nobasistransform'
  else
   write(fid1,'(A,38X,A)') 'Route','C   N=           5'
   write(fid1,'(A)') '#p HF chkbasis nosymm int(nobasistransform'//TRIM(ham(irel))
  end if
  write(fid1,'(A)') TRIM(buf)

 else                          ! there is Route in .fch(k)
  do while(.true.)
   read(fid,'(A)') buf
   write(fid1,'(A)') TRIM(buf)
   if(buf(1:5) == 'Route') exit
  end do ! for while

  i = INDEX(route,' int(')
  if(i == 0) i = INDEX(route,' int=')
  if(i == 0) then
   write(6,'(/,A)') "ERROR in subroutine set_irel_in_fch: neither 'int(' nor 'i&
                    &nt=' is found"
   write(6,'(A)') 'in file '//TRIM(fchname)
   close(fid)
   close(fid1,status='delete')
   stop
  end if
  j = INDEX(route(i+6:),' ')
  route = route(1:i)//'int(nobasistransform'//TRIM(ham(irel))//TRIM(route(i+j+5:))
  write(fid1,'(A)') TRIM(route)

  do while(.true.)
   read(fid,'(A)') buf
   if(buf(1:6) == 'Charge') exit
  end do ! for while
  write(fid1,'(A)') TRIM(buf)
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(fchname1), TRIM(fchname))
end subroutine set_irel_in_fch

! check whether 'nobasistransform' exists in a given .fch(k) file
! if not, print warning
subroutine check_nobasistransform_in_fch(fchname)
 implicit none
 integer:: i, fid
 character(len=240) :: buf
 character(len=1200) :: longbuf
 character(len=240), intent(in) :: fchname

 buf = ' '; longbuf = ' '
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:5) == 'Route') exit

  ! no Route Section in fch, maybe produced by g09 or older version
  ! do not print warning in this case
  if(buf(1:5) == 'Charg') then
   close(fid)
   return
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine check_nobasistransform_in_fch: problema&
                   &tic file '//TRIM(fchname)
  close(fid)
  stop
 end if

 do i = 1, 5
  read(fid,'(A)') buf
  if(buf(1:5) == 'Charg') exit

  if(i == 1) then
   longbuf = buf
  else
   longbuf = TRIM(longbuf)//TRIM(buf)
  end if
 end do ! for i

 close(fid)
 call lower(longbuf)

 if(index(longbuf,'nobasistransform') == 0) then
  write(6,'(/,A)') REPEAT('-',56)
  write(6,'(A)') "Warning in subroutine check_nobasistransform_in_fch: keyword &
                 &'int=nobasistransform'"
  write(6,'(A)') 'not detected in file '//TRIM(fchname)
  write(6,'(A)') 'It might be dangerous to transfer orbitals if you did not spe&
                 &cify this keyword'
  write(6,'(A)') 'in .gjf file.'
  write(6,'(A)') REPEAT('-',56)
 end if
end subroutine check_nobasistransform_in_fch

! check whether 'nosymm' exists in a given .fch(k) file
! if not, print warning
subroutine check_nosymm_in_fch(fchname)
 implicit none
 integer:: i, fid
 character(len=240) :: buf
 character(len=1200) :: longbuf
 character(len=240), intent(in) :: fchname

 buf = ' '; longbuf = ' '
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:5) == 'Route') exit

  ! no Route Section in fch, maybe produced by g09 or older version
  ! do not print warning in this case
  if(buf(1:5) == 'Charg') then
   close(fid)
   return
  end if
 end do ! for while

 if(i /= 0) then
  close(fid)
  write(6,'(/,A)') 'ERROR in subroutine check_nosymm_in_fch: problematic file '&
                   //TRIM(fchname)
  stop
 end if

 do i = 1, 5
  read(fid,'(A)') buf
  if(buf(1:5) == 'Charg') exit

  if(i == 1) then
   longbuf = buf
  else
   longbuf = TRIM(longbuf)//TRIM(buf)
  end if
 end do ! for i

 close(fid)
 call lower(longbuf)

 if(index(longbuf,'nosymm') == 0) then
  write(6,'(/,A)') REPEAT('-',56)
  write(6,'(A)') "Warning in subroutine check_nosymm_in_fch: keyword 'nosymm' n&
                 &ot detected in"
  write(6,'(A)') 'file '//TRIM(fchname)
  write(6,'(A)') 'It might be dangerous to transfer orbitals if you did not spe&
                 &cify this keyword'
  write(6,'(A)') 'in .gjf file.'
  write(6,'(A)') REPEAT('-',56)
 end if
end subroutine check_nosymm_in_fch

!! expand the AO overlap integral matrix from spherical harmonic functions to
!! pure Cartesian functions
!subroutine ovlp_sph2cart(ncontr, shltyp, nbf0, nbf1, sph_ovlp, cart_ovlp)
! use parameter_sph2cart, only: rd, rf, rg, rh
! implicit none
! integer :: i, j, k, m, nbf
! integer, intent(in) :: ncontr, nbf0, nbf1
! integer, intent(in) :: shltyp(ncontr)
! real(kind=8), parameter :: r1=10d0/7d0, r2=243d0/140d0, r3=9d0/7d0, &
!  r4=5d0/3d0, r5=205d0/84d0, r6=51d0/28d0
! real(kind=8), parameter :: norm_g(2:14) = [r1,r2,r1,1d0,r1,r3,r3,r1,r2,r3,r2,r1,r1]
! real(kind=8), parameter :: norm_h(2:20) = [r4,r5,r5,r4,1d0,r4,r4,r6,r4,r4,r5,&
!  r6,r6,r5,r5,r4,r5,r4,r4]
! real(kind=8), intent(in) :: sph_ovlp(nbf0,nbf0)
! real(kind=8), intent(out) :: cart_ovlp(nbf1,nbf1)
! real(kind=8), allocatable :: rdt(:,:), rft(:,:), rgt(:,:), rht(:,:)
!
! allocate(rdt(5,6), source=TRANSPOSE(rd))
! allocate(rft(7,10), source=TRANSPOSE(rf))
! rft(:,4:9) = rft(:,4:9)/DSQRT(1.2d0)
!
! allocate(rgt(9,15), source=TRANSPOSE(rg))
! do i = 2, 14
!  rgt(:,i) = rgt(:,i)/DSQRT(norm_g(i))
! end do ! for i
!
! allocate(rht(11,21), source=TRANSPOSE(rh))
! do i = 2, 20
!  rht(:,i) = rht(:,i)/DSQRT(norm_h(i))
! end do ! for i
!
! nbf = 0; j = 0; cart_ovlp = 0d0
!
! do i = 1, ncontr, 1
!  k = nbf + 1; m = j + 1
!
!  select case(shltyp(i))
!  case( 0) ! S
!   cart_ovlp(k,k) = sph_ovlp(m,m)
!   nbf = nbf + 1; j = j + 1
!  case( 1) ! P
!   cart_ovlp(k:k+2,k:k+2) = sph_ovlp(m:m+2,m:m+2)
!   nbf = nbf + 3; j = j + 3
!  case(-1) ! L, SP
!   cart_ovlp(k:k+3,k:k+3) = sph_ovlp(m:m+3,m:m+3)
!   nbf = nbf + 4; j = j + 4
!  case(-2) ! 5D
!   call calc_CTSC(5, 6, rdt, sph_ovlp(m:m+4,m:m+4), cart_ovlp(k:k+5,k:k+5))
!   nbf = nbf + 6; j = j + 5
!  case(-3) ! 7F
!   call calc_CTSC(7, 10, rft, sph_ovlp(m:m+6,m:m+6), cart_ovlp(k:k+9,k:k+9))
!   nbf = nbf + 10; j = j + 7
!  case(-4) ! 9G
!   call calc_CTSC(9, 15, rgt, sph_ovlp(m:m+8,m:m+8), cart_ovlp(k:k+14,k:k+14))
!   nbf = nbf + 15; j = j + 9
!  case(-5) ! 11H
!   call calc_CTSC(11,21,rht, sph_ovlp(m:m+10,m:m+10), cart_ovlp(k:k+20,k:k+20))
!   nbf = nbf + 21; j = j + 11
!  case default
!   write(6,'(/,A)') 'ERROR in subroutine ovlp_sph2cart: shltyp(i) out of range!'
!   write(6,'(A,3I6)') 'i, ncontr, shltyp(i)=', i, ncontr, shltyp(i)
!   stop
!  end select
! end do ! for i
!
! deallocate(rdt, rft, rgt, rht)
! if(nbf/=nbf1 .or. j/=nbf0) then
!  write(6,'(/,A)') 'ERROR in subroutine ovlp_sph2cart: nbf/=nbf1 or j/=nbf0.'
!  write(6,'(A,4I6)') 'j, nbf, nbf0, nbf1=', j, nbf, nbf0, nbf1
!  stop
! end if
!end subroutine ovlp_sph2cart

! expand MO coefficients from spherical harmonic functions to pure Cartesian
! functions
subroutine mo_sph2cart(ncontr, shltyp, nbf0, nbf1, nmo, sph_coeff, cart_coeff)
 use parameter_sph2cart, only: rd, rf, rg, rh
 implicit none
 integer :: i, j, nbf
 integer, intent(in) :: ncontr, nbf0, nbf1, nmo
 integer, intent(inout) :: shltyp(ncontr)
 real(kind=8), intent(in) :: sph_coeff(nbf0,nmo)
 real(kind=8), intent(out) :: cart_coeff(nbf1,nmo)
 
 nbf = 0; j = 0; cart_coeff = 0d0

 do i = 1, ncontr, 1
  select case(shltyp(i))
  case( 0) ! S
   cart_coeff(nbf+1,:) = sph_coeff(j+1,:)
   nbf = nbf + 1; j = j + 1
  case( 1) ! P
   cart_coeff(nbf+1:nbf+3,:) = sph_coeff(j+1:j+3,:)
   nbf = nbf + 3; j = j + 3
  case(-1) ! L, SP
   cart_coeff(nbf+1:nbf+4,:) = sph_coeff(j+1:j+4,:)
   nbf = nbf + 4; j = j + 4
  case(-2) ! 5D
   cart_coeff(nbf+1:nbf+6,:) = MATMUL(rd, sph_coeff(j+1:j+5,:))
   nbf = nbf + 6; j = j + 5
   shltyp(i) = 2
  case(-3) ! 7F
   cart_coeff(nbf+1:nbf+10,:) = MATMUL(rf, sph_coeff(j+1:j+7,:))
   nbf = nbf + 10; j = j + 7
   shltyp(i) = 3
  case(-4) ! 9G
   cart_coeff(nbf+1:nbf+15,:) = MATMUL(rg, sph_coeff(j+1:j+9,:))
   nbf = nbf + 15; j = j + 9
   shltyp(i) = 4
  case(-5) ! 11H
   cart_coeff(nbf+1:nbf+21,:) = MATMUL(rh, sph_coeff(j+1:j+11,:))
   nbf = nbf + 21; j = j + 11
   shltyp(i) = 5
  case default
   write(6,'(/,A)') 'ERROR in subroutine mo_sph2cart: shltyp(i) out of range!'
   write(6,'(A,3I6)') 'i, ncontr, shltyp(i)=', i, ncontr, shltyp(i)
   stop
  end select
 end do ! for i

 if(nbf/=nbf1 .or. j/=nbf0) then
  write(6,'(/,A)') 'ERROR in subroutine mo_sph2cart: nbf/=nbf1 or j/=nbf0.'
  write(6,'(A,4I6)') 'j, nbf, nbf0, nbf1=', j, nbf, nbf0, nbf1
  stop
 end if
end subroutine mo_sph2cart

! create/generate a new .fch file using arrays stored in memory
subroutine write_fch(fchname)
 use fch_content
 implicit none
 integer :: i, j, k, ncoeff, nhigh, ncontr_l, len_dm, len_fc, fid
 integer, allocatable :: ilsw(:)
 real(kind=8), allocatable :: shell_coor(:), coor1(:)
 character(len=4) :: hf_str
 character(len=240), intent(in) :: fchname
 logical :: uhf

! if(allocated(eigen_e_b)) then
!  uhf = .true.
! else
!  uhf= .false.
! end if
 uhf = is_uhf
 ncoeff = nif*nbf

 i = MAXVAL(shell_type)
 if(i < 0) i = -i
 j = MINVAL(shell_type)
 if(j < 0) j = -j
 nhigh = MAX(i,j)
 ncontr_l = MAXVAL(prim_per_shell)

 ! Note: for semi-empirical methods, the array rnuc may be allocated outside of
 ! this subroutine, and used to keep the effective nuclear charges.
 ! For semi-empirical methods, there is usually no ECP/PP for core electrons,
 ! yet with GTO for valence electrons.
 if(.not. allocated(rnuc)) then
  allocate(rnuc(natom))
  forall(i = 1:natom) rnuc(i) = DBLE(ielem(i))
  ! If ECP/PP is used, the pseudo core electrons need to be substracted
  if(allocated(RNFroz)) rnuc = rnuc - RNFroz
 end if

 allocate(coor1(3*natom))
 forall(i = 1:natom) coor1(3*i-2:3*i) = coor(:,i)
 coor1 = coor1/Bohr_const

 allocate(shell_coor(3*ncontr))
 call calc_shell_coor(ncontr, natom, shell2atom_map, coor, shell_coor)

 allocate(ilsw(100), source=0)
 if(uhf) ilsw(1) = 1
 ilsw(5) = 2
 ilsw(12) = -1
 ilsw(13) = 5
 ilsw(25:26) = 1
 ilsw(33) = 100000
 ilsw(35) = -1
 ilsw(46) = 1
 ilsw(50) = -2000000000
 ilsw(51) = 1
 ilsw(57) = 4
 ilsw(58) = 52
 ilsw(70) = natom
 ilsw(72) = 1

 hf_str = 'RHF '
 if(uhf) then
  hf_str = 'UHF '
 else
  if(mult /= 1) hf_str = 'ROHF'
 end if

 open(newunit=fid,file=TRIM(fchname),status='replace')
 write(fid,'(A)') 'Generated by MOKIT'
 write(fid,'(A2,8X,A4,56X,A5)') 'SP', hf_str, 'MYBAS'
 write(fid,'(A,I17)') 'Number of atoms                            I',natom

 ! The Route Section is not necessary for a fch file. But when the DKH/X2C key-
 ! words are written into the fch file, fch2xxx utilities will recognize it and
 ! write related relativistic keywords into generated files
 if(irel /= -1) then
  write(fid,'(A,38X,A,3X,A,I12)') 'Route','C','N=',5
  write(fid,'(A)',advance='no') '#p '//TRIM(hf_str)//' chkbasis nosymm int(noba&
                                &sistransform,'
  select case(irel)
  case(-3)
   write(fid,'(A)') 'X2C)'
  case(-2)
   write(fid,'(A)') 'RESC)'
  case(0)
   write(fid,'(A)') 'DKH0)'
  case(2)
   write(fid,'(A)') 'DKH2)'
  case(4)
   write(fid,'(A)') 'DKHSO)'
  case default
   write(6,'(/,A)') 'ERROR in subroutine write_fch: irel out of range.'
   write(6,'(A,I0)') 'irel=', irel
   stop
  end select
 end if

 write(fid,'(A,I17)') 'Charge                                     I',charge
 write(fid,'(A,I17)') 'Multiplicity                               I',mult
 write(fid,'(A,I17)') 'Number of electrons                        I',na+nb
 write(fid,'(A,I17)') 'Number of alpha electrons                  I',na
 write(fid,'(A,I17)') 'Number of beta electrons                   I',nb
 write(fid,'(A,I17)') 'Number of basis functions                  I',nbf
 write(fid,'(A,I17)') 'Number of independent functions            I',nif
 write(fid,'(A,I12)') 'Atomic numbers                             I   N=',natom
 write(fid,'(6I12)') ielem
 write(fid,'(A,I12)') 'Nuclear charges                            R   N=',natom
 write(fid,'(5(1X,ES15.8))') rnuc
 deallocate(rnuc)
 write(fid,'(A,I12)') 'Current cartesian coordinates              R   N=',natom*3
 write(fid,'(5(1X,ES15.8))') coor1
 deallocate(coor1)
 if(.not. allocated(iatom_type)) allocate(iatom_type(natom), source=0)
 write(fid,'(A,I12)') 'Int Atom Types                             I   N=',natom
 write(fid,'(6I12)') iatom_type
 write(fid,'(A,I12)') 'Real atomic weights                        R   N=',natom
 write(fid,'(5(1X,ES15.8))') (ram(ielem(i)),i=1,natom)
 write(fid,'(A,I17)') 'Number of contracted shells                I',ncontr
 write(fid,'(A,I17)') 'Number of primitive shells                 I',nprim
 write(fid,'(A,I17)') 'Highest angular momentum                   I',nhigh
 write(fid,'(A,I17)') 'Largest degree of contraction              I',ncontr_l
 write(fid,'(A,I12)') 'Shell types                                I   N=',ncontr
 write(fid,'(6I12)') shell_type
 write(fid,'(A,I12)') 'Number of primitives per shell             I   N=',ncontr
 write(fid,'(6I12)') prim_per_shell
 write(fid,'(A,I12)') 'Shell to atom map                          I   N=',ncontr
 write(fid,'(6I12)') shell2atom_map
 write(fid,'(A,I12)') 'Primitive exponents                        R   N=',nprim
 write(fid,'(5(1X,ES15.8))') prim_exp
 write(fid,'(A,I12)') 'Contraction coefficients                   R   N=',nprim
 write(fid,'(5(1X,ES15.8))') contr_coeff
 if(allocated(contr_coeff_sp)) then
  write(fid,'(A,I12)') 'P(S=P) Contraction coefficients            R   N=',nprim
  write(fid,'(5(1X,ES15.8))') contr_coeff_sp
 end if
 write(fid,'(A,I12)') 'Coordinates of each shell                  R   N=',3*ncontr
 write(fid,'(5(1X,ES15.8))') shell_coor
 deallocate(shell_coor)
 write(fid,'(A)') 'Num ILSW                                   I              100'
 write(fid,'(A)') 'ILSW                                       I   N=         100'
 write(fid,'(6I12)') ilsw
 deallocate(ilsw)
 if(LenNCZ > 0) then
  write(fid,'(A)') 'ECP-MxAtEC                                 I           250000'
  write(fid,'(A)') 'ECP-MaxLECP                                I               10'
  write(fid,'(A)') 'ECP-MaxAtL                                 I          2250000'
  write(fid,'(A)') 'ECP-MxTECP                                 I         11250000'
  write(fid,'(A,I17)') 'ECP-LenNCZ                                 I', LenNCZ
  write(fid,'(A,I12)') 'ECP-KFirst                                 I   N=',10*natom
  write(fid,'(6I12)') KFirst
  write(fid,'(A,I12)') 'ECP-KLast                                  I   N=',10*natom
  write(fid,'(6I12)') KLast
  write(fid,'(A,I12)') 'ECP-LMax                                   I   N=',natom
  write(fid,'(6I12)') Lmax
  write(fid,'(A,I12)') 'ECP-LPSkip                                 I   N=',natom
  write(fid,'(6I12)') LPSkip
  write(fid,'(A,I12)') 'ECP-RNFroz                                 R   N=',natom
  write(fid,'(5(1X,ES15.8))') RNFroz
  write(fid,'(A,I12)') 'ECP-NLP                                    I   N=',LenNCZ
  write(fid,'(6I12)') NLP
  write(fid,'(A,I12)') 'ECP-CLP1                                   R   N=',LenNCZ
  write(fid,'(5(1X,ES15.8))') CLP
  write(fid,'(A,I12)') 'ECP-CLP2                                   R   N=',LenNCZ
  if(.not. allocated(CLP2)) allocate(CLP2(LenNCZ), source=0d0)
  ! CLP2 is used for SO-ECP, not supported currently, set to zero
  write(fid,'(5(1X,ES15.8))') CLP2
  write(fid,'(A,I12)') 'ECP-ZLP                                    R   N=',LenNCZ
  write(fid,'(5(1X,ES15.8))') ZLP
 end if
 write(fid,'(A,E27.15)') 'Virial Ratio                               R',virial
 write(fid,'(A,E27.15)') 'Total Energy                               R',tot_e
 write(fid,'(A,I12)') 'Alpha Orbital Energies                     R   N=',nif
 write(fid,'(5(1X,ES15.8))') eigen_e_a
 if(uhf) then
  write(fid,'(A,22X,A,I12)') 'Beta Orbital Energies','R   N=', nif
  write(fid,'(5(1X,ES15.8))') eigen_e_b
 end if

 write(fid,'(A,I12)') 'Alpha MO coefficients                      R   N=',ncoeff
 write(fid,'(5(1X,ES15.8))') alpha_coeff
 if(uhf) then
  write(fid,'(A,23X,A,I12)') 'Beta MO coefficients','R   N=', ncoeff
  write(fid,'(5(1X,ES15.8))') beta_coeff
 end if

 len_dm = nbf*(nbf+1)/2
 write(fid,'(A,I12)') 'Total SCF Density                          R   N=',len_dm
 if(allocated(tot_dm)) then
  write(fid,'(5(1X,ES15.8))') ((tot_dm(j,i),j=1,i),i=1,nbf)
 else
  write(fid,'(5(1X,ES15.8))') ((0d0,j=1,i),i=1,nbf) ! ugly workaround for EOF
 end if

 if(allocated(spin_dm)) then
  write(fid,'(A,I12)') 'Spin SCF Density                           R   N=',len_dm
  write(fid,'(5(1X,ES15.8))') ((spin_dm(j,i),j=1,i),i=1,nbf)
 end if

 if(allocated(mull_char)) then
  write(fid,'(A,I12)') 'Mulliken Charges                           R   N=',natom
  write(fid,'(5(1X,ES15.8))') mull_char
 end if

 if(allocated(force_const)) then
  k = 3*natom
  len_fc = k*(k+1)/2
  write(fid,'(A,I12)') 'Cartesian Force Constants                  R   N=',len_fc
  write(fid,'(5(1X,ES15.8))') ((force_const(j,i),j=1,i),i=1,k)
 end if

 if(nmode > 0) then
  write(fid,'(A,9X,I8)') 'Number of Normal Modes                     I',nmode
  write(fid,'(A)') 'Vib-NDFDPol                                I                0'
  write(fid,'(A)') 'Vib-NumROA                                 I                0'
  write(fid,'(A,9X,I8)') 'Vib-NDim                                   I',nmode
  write(fid,'(A,9X,I8)') 'Vib-NDim0                                  I',nmode
  write(fid,'(A,4X,I8)') 'Vib-AtMass                                 R   N=',natom
  write(fid,'(5(1X,ES15.8))') (ram(ielem(i)),i=1,natom)
  write(fid,'(A,4X,I8)') 'Vib-E2                                     R   N=',nmode*14
  write(fid,'(5(1X,ES15.8))') vibe2
  j = nmode*3*natom
  write(fid,'(A,4X,I8)') 'Vib-Modes                                  R   N=',j
  write(fid,'(5(1X,ES15.8))') norm_mode
 end if

 close(fid)
end subroutine write_fch

! calculate/determine the number of core orbitals
subroutine calc_ncore(fchname, chem_core, ecp_core)
 use fch_content, only: core_orb, RNFroz
 implicit none
 integer :: i, natom, fid
 integer, intent(out) :: chem_core, ecp_core
 integer, allocatable :: nuc(:)
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
 logical :: ecp

 chem_core = 0; buf = ' ' ! initialization

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:14) == 'Number of atom') exit
 end do ! for i
 read(buf(51:),*) natom

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:14) == 'Atomic numbers') exit
 end do ! for while
 allocate(nuc(natom), source=0)
 read(fid,'(6(1X,I11))') (nuc(i), i=1,natom)

 ecp = .false.
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:10) == 'ECP-RNFroz') then
   ecp = .true.
   exit
  end if
  if(buf(1:7) == 'Alpha O') then
   close(fid)
   exit
  end if
 end do ! for while

 if(ecp) then
  if(allocated(RNFroz)) deallocate(RNFroz)
  allocate(RNFroz(natom), source=0d0)
  read(fid,'(5(1X,ES15.8))') (RNFroz(i), i=1,natom)
  close(fid)
  ecp_core = NINT(0.5d0*SUM(RNFroz)) ! half of core electrons
  ! Sometimes the user would use large-core ECP/PP, so here the larger one
  ! between core_orb(nuc(i)) and RNFroz(i) should be taken. And chem_core
  ! may be large than sum_i core_orb(nuc(i))
  do i = 1, natom, 1
   chem_core = chem_core + MAX(core_orb(nuc(i)), NINT(0.5d0*RNFroz(i)))
  end do ! for i
  deallocate(RNFroz)
 else ! all-electron basis set
  do i = 1, natom, 1
   chem_core = chem_core + core_orb(nuc(i))
  end do ! for i
 end if

 deallocate(nuc)
end subroutine calc_ncore

! read the position marks of (P,SP)/5D/6D/7F/10F/9G/15G/11H/21H from the array
! shell_type
subroutine read_pdfgh_mark_from_shltyp(ncontr, shltyp, no_mark, mark)
 implicit none
 integer :: i, nbf
 integer, intent(in) :: ncontr
 integer, intent(in) :: shltyp(ncontr)
 integer, intent(out) :: no_mark(9) ! No. of marks per angular momentum
 integer, intent(out) :: mark(ncontr,9)
 ! 9: P, 5D, 6D, 7F, 10F, 9G, 15G, 11H, 21H

 nbf = 0; no_mark = 0; mark = 0

 do i = 1, ncontr, 1
  select case(shltyp(i))
  case( 0)   ! S
   nbf = nbf + 1
  case( 1)   ! 3P
   no_mark(1) = no_mark(1) + 1
   mark(no_mark(1),1) = nbf + 1
   nbf = nbf + 3
  case(-1)   ! SP or L
   no_mark(1) = no_mark(1) + 1
   mark(no_mark(1),1) = nbf + 2
   nbf = nbf + 4
  case(-2)   ! 5D
   no_mark(2) = no_mark(2) + 1
   mark(no_mark(2),2) = nbf + 1
   nbf = nbf + 5
  case( 2)   ! 6D
   no_mark(3) = no_mark(3) + 1
   mark(no_mark(3),3) = nbf + 1
   nbf = nbf + 6
  case(-3)   ! 7F
   no_mark(4) = no_mark(4) + 1
   mark(no_mark(4),4) = nbf + 1
   nbf = nbf + 7
  case( 3)   ! 10F
   no_mark(5) = no_mark(5) + 1
   mark(no_mark(5),5) = nbf + 1
   nbf = nbf + 10
  case(-4)   ! 9G
   no_mark(6) = no_mark(6) + 1
   mark(no_mark(6),6) = nbf + 1
   nbf = nbf + 9
  case( 4)   ! 15G
   no_mark(7) = no_mark(7) + 1
   mark(no_mark(7),7) = nbf + 1
   nbf = nbf + 15
  case(-5)   ! 11H
   no_mark(8) = no_mark(8) + 1
   mark(no_mark(8),8) = nbf + 1
   nbf = nbf + 11
  case( 5)   ! 21H
   no_mark(9) = no_mark(9) + 1
   mark(no_mark(9),9) = nbf + 1
   nbf = nbf + 21
  case default
   write(6,'(/,A)') 'ERROR in subroutine read_pdfgh_mark_from_shltyp:'
   write(6,'(A,I0)') 'Invalid shltyp(i)=', shltyp(i)
   stop
  end select
 end do ! for i
end subroutine read_pdfgh_mark_from_shltyp

! Find the number of kinds of elements in the array nuc. E.g. {H,C,O} -> 3
subroutine find_nelem_in_nuc(natom, nuc, nelem)
 use fch_content, only: period_nelem
 implicit none
 integer :: i
 integer, intent(in) :: natom
 integer, intent(in) :: nuc(natom)
 integer, intent(out) :: nelem

 nelem = 0
 do i = 1, period_nelem, 1
  if(COUNT(nuc == i) > 0) nelem = nelem + 1
 end do ! for i
end subroutine find_nelem_in_nuc

! Find the number of atoms for each kind of element. Note: the input nuc should
! have beend sorted in ascending order.
subroutine find_natom_per_elem(natom, nuc, nelem, natom_per_elem)
 implicit none
 integer :: i, j
 integer, intent(in) :: natom, nelem
 integer, intent(in) :: nuc(natom)
 integer, intent(out) :: natom_per_elem(nelem)

 natom_per_elem = 0; natom_per_elem(1) = 1; j = 1

 do i = 2, natom, 1
  if(nuc(i) /= nuc(i-1)) j = j + 1
  natom_per_elem(j) = natom_per_elem(j) + 1
 end do ! for i
end subroutine find_natom_per_elem

