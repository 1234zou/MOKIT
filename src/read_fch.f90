! moved by jxzou at 20200322 from fch2mkl.f90 to this new file

! The 'Shell types' array in Gaussian .fch file:
!
!   Spherical     |     Cartesian
! -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5
!  H  G  F  D  L  S  P  D  F  G  H
!
! 'L' is 'SP' in Pople-type basis sets

module fch_content
 implicit none
 integer :: nbf, nif         ! number of basis functions/MOs
 integer :: na, nb, nopen    ! number of alpha/beta/open shell electrons
 integer :: ncontr, nprim    ! Number of contracted/primitive shells
 integer :: charge, mult     ! charge and spin multiplicity
 integer :: natom            ! number of atoms
 integer :: LenNCZ           ! ECP-LenNCZ
 integer, parameter :: iout = 6              ! print id, default on screen
 integer, allocatable :: ielem(:)            ! elements, 6 for 'C', etc
 integer, allocatable :: shell_type(:)       ! Shell types
 integer, allocatable :: prim_per_shell(:)   ! Number of primitives per shell
 integer, allocatable :: shell2atom_map(:)   ! Shell to atom map
 integer, allocatable :: KFirst(:,:)         ! ECP-KFirst, size (natom,10)
 integer, allocatable :: KLast(:,:)          ! ECP-KLast, size (natom,10)
 integer, allocatable :: Lmax(:), LPSkip(:)  ! ECP-LMax, ECP-LPSkip, size natom
 integer, allocatable :: NLP(:)              ! ECP-NLP, size LenNCZ
 real(kind=8), parameter :: Bohr_const = 0.52917721092d0
 real(kind=8), allocatable :: coor(:,:)                  ! Cartesian coordinates
 real(kind=8), allocatable :: prim_exp(:)                ! Primitive exponents
 real(kind=8), allocatable :: contr_coeff(:)             ! Contraction coefficients
 real(kind=8), allocatable :: contr_coeff_sp(:)          ! (S=P) Contraction coefficients
 real(kind=8), allocatable :: eigen_e_a(:), eigen_e_b(:) ! Alpha/Beta energy levels
 real(kind=8), allocatable :: alpha_coeff(:,:), beta_coeff(:,:) ! Alpha/Beta MOs
 real(kind=8), allocatable :: RNFroz(:)                  ! ECP-RNFroz(core e), size natom
 real(kind=8), allocatable :: CLP(:), ZLP(:)             ! ECP-CLP1, ECP-ZLP, size LenNCZ
 character(len=2), allocatable :: elem(:)                ! elements ('H ', 'C ', etc)

 integer, parameter :: period_nelem = 112
 character(len=2), parameter :: period_elem(period_nelem) = &
  ['H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
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
   'Rg', 'Cn' ]

contains

! read geometry, basis sets, ECP(if any) and MOs from .fch file (Gaussian)
subroutine read_fch(fchname, uhf)
 implicit none
 integer :: i, fid
 integer :: ncoeff ! = nbf*nif
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
 real(kind=8), allocatable :: coor0(:) ! Cartesian coordinates
 real(kind=8), allocatable :: coeff(:) ! MOs, nbf*nif
 logical, intent(in) :: uhf
 logical :: ecp

 i = INDEX(fchname,'.fch',back=.true.)
 if(i == 0) then
  write(iout,'(A)') "ERROR in subroutine read_fch: input filename does not&
                   & contain '.fch' suffix!"
  write(iout,'(A)') 'fchname='//TRIM(fchname)
  stop
 end if

 buf = ' '
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
  write(iout,'(A)') 'ERROR in subroutine read_fch: The number of alpha electrons&
                   & is less than that of beta electrons!'
  write(iout,'(A)') 'fchname='//TRIM(fchname)
  write(iout,'(2(A,I0))') 'na=', na, ', nb=', nb
  stop
 end if

 ! find and read natom, elem
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:14) == 'Atomic numbers') exit
 end do
 BACKSPACE(fid)
 read(fid,'(A49,2X,I10)') buf, natom
 allocate(ielem(natom), source=0)
 read(fid,'(6(1X,I11))') (ielem(i),i=1,natom)
 allocate(elem(natom))
 elem = ' '
 forall(i = 1:natom) elem(i) = nuc2elem(ielem(i))

 ! find and read coordinates
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:12) == 'Current cart') exit
 end do
 allocate(coor0(3*natom), source=0.0d0)
 allocate(coor(3,natom), source=0.0d0)
 read(fid,'(5(1X,ES15.8))') (coor0(i),i=1,3*natom)
 coor = RESHAPE(coor0,(/3,natom/))
 deallocate(coor0)
 coor = coor*Bohr_const ! convert Bohr to Anstrom

 ! find and read variables ncontr, nprim
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:18) == 'Number of contract') exit
 end do
 BACKSPACE(fid)
 read(fid,'(A49,2X,I10)') buf, ncontr
 read(fid,'(A49,2X,I10)') buf, nprim

 ! find and read Shell types, Number of primitives per shell, Shell to atom map,
 !  Primitive exponents, Contraction coefficients
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:14) == 'Shell types') exit
 end do
 allocate(shell_type(ncontr), source=0)
 read(fid,'(6(1X,I11))') (shell_type(i), i=1,ncontr)

 read(fid,'(A)') buf
 allocate(prim_per_shell(ncontr), source=0)
 read(fid,'(6(1X,I11))') (prim_per_shell(i), i=1,ncontr)

 read(fid,'(A)') buf
 allocate(shell2atom_map(ncontr), source=0)
 read(fid,'(6(1X,I11))') (shell2atom_map(i), i=1,ncontr)

 read(fid,'(A)') buf
 allocate(prim_exp(nprim), source=0.0d0)
 read(fid,'(5(1X,ES15.8))') (prim_exp(i),i=1,nprim)

 read(fid,'(A)') buf
 allocate(contr_coeff(nprim), source=0.0d0)
 read(fid,'(5(1X,ES15.8))') (contr_coeff(i),i=1,nprim)

 read(fid,'(A)') buf
 if(buf(1:6) == 'P(S=P)') then
  allocate(contr_coeff_sp(nprim), source=0.0d0)
  read(fid,'(5(1X,ES15.8))') (contr_coeff_sp(i),i=1,nprim)
 end if

 ! check if there is any ECP/PP data
 ! Note that def2TZVP for the 4th row will cause LenNCZ=0
 ecp = .false.
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
  allocate(RNFroz(natom), source=0.0d0)
  allocate(NLP(LenNCZ), source=0)
  allocate(CLP(LenNCZ), source=0.0d0)
  allocate(ZLP(LenNCZ), source=0.0d0)

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
  do while(.true.)
   read(fid,'(A)') buf
   if(buf(1:7) == 'ECP-ZLP') exit
  end do
  read(fid,'(5(1X,ES15.8))') (ZLP(i), i=1,LenNCZ)
 end if

 ! find and read Alpha Orbital Energies
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:13) == 'Alpha Orbital') exit
 end do
 allocate(eigen_e_a(nif), source=0.0d0)
 read(fid,'(5(1X,ES15.8))') (eigen_e_a(i),i=1,nif)

 ! if '-uhf' specified, read Beta Orbital Energies
 if(uhf) then
  read(fid,'(A)') buf
  if(buf(1:12) /= 'Beta Orbital') then
   write(iout,'(A)') "ERROR in subroutine read_fch: no 'Beta Orbital' found in&
                    & this .fch(k) file, but you specify '-uhf'."
   write(iout,'(A)') 'fchname='//TRIM(fchname)
   stop
  end if
  allocate(eigen_e_b(nif), source=0.0d0)
  read(fid,'(5(1X,ES15.8))') (eigen_e_b(i),i=1,nif)
 end if

 ! read Alpha MO
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:8) == 'Alpha MO') exit
 end do
 allocate(coeff(ncoeff), source=0.0d0)
 read(fid,'(5(1X,ES15.8))') (coeff(i),i=1,ncoeff)
 allocate(alpha_coeff(nbf,nif), source=0.0d0)
 alpha_coeff = RESHAPE(coeff,(/nbf,nif/))
 deallocate(coeff)

 ! if '-uhf' specified, read Beta MO
 if(uhf) then
  allocate(coeff(ncoeff), source=0.0d0)
  read(fid,'(A)') buf   ! skip the Beta MO line
  read(fid,'(5(1X,ES15.8))') (coeff(i),i=1,ncoeff)
  allocate(beta_coeff(nbf,nif), source=0.0d0)
  beta_coeff = RESHAPE(coeff,(/nbf,nif/))
  deallocate(coeff)
 end if

 close(fid)
 return
end subroutine read_fch

 ! map a nuclear charge to an element (e.g. 6->'C')
 ! Note: only 1-112 elements are supported!
 pure function nuc2elem(i) result(s)
  implicit none
  integer, intent(in) :: i
  character(len=2) :: s

  s = period_elem(i)
  return
 end function nuc2elem

 ! map an element to a nuclear charge (e.g. 'C'->6)
 ! Note: only 1-112 elements are supported!
 pure function elem2nuc(s) result(i)
  implicit none
  integer :: i
  character(len=2), intent(in) :: s

  do i = 1, period_nelem, 1
   if(period_elem(i) == s) return
  end do
  return
 end function elem2nuc

end module fch_content

