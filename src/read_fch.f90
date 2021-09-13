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
 implicit none
 integer :: nbf, nif         ! number of basis functions/MOs
 integer :: na, nb, nopen    ! number of alpha/beta/open shell electrons
 integer :: ncontr, nprim    ! Number of contracted/primitive shells
 integer :: charge, mult     ! charge and spin multiplicity
 integer :: natom            ! number of atoms
 integer :: LenNCZ           ! ECP-LenNCZ
 integer, parameter :: iout = 6              ! print id, default on screen
 integer, parameter :: period_nelem = 112    ! 112 elements, H-Cn
 integer, allocatable :: ielem(:)            ! elements, 6 for 'C', etc
 integer, allocatable :: shell_type(:)       ! Shell types
 integer, allocatable :: prim_per_shell(:)   ! Number of primitives per shell
 integer, allocatable :: shell2atom_map(:)   ! Shell to atom map
 integer, allocatable :: KFirst(:,:)         ! ECP-KFirst, size (natom,10)
 integer, allocatable :: KLast(:,:)          ! ECP-KLast, size (natom,10)
 integer, allocatable :: Lmax(:), LPSkip(:)  ! ECP-LMax, ECP-LPSkip, size natom
 integer, allocatable :: NLP(:)              ! ECP-NLP, size LenNCZ
 real(kind=8), parameter :: Bohr_const = 0.52917721092d0
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
 character(len=2), allocatable :: elem(:)       ! elements ('H ', 'C ', etc)

 ! Frozen core orbitals used in RHF_proj and dynamic correlation computations.
 ! This table is copied from the figure in ORCA 5.0.1 manual 9.11 Frozen Core
 ! Options.
 ! Note: frozen electrons are two times of the frozen core orbitals
 integer, parameter :: core_orb(period_nelem) = &
 (/  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  5,  5,  5, &
     5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5, &
     9,  9,  9,  9,  9,  9,  9,  9, 14, 14, 14, 14, 14, 14, 14, &
    14, 14, 14, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, &
    18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 23, 23, 23, 23, 23, &
    23, 23, 23, 23, 23, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, &
    34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 50, 50, &
    50, 50, 50, 50, 50, 50, 50/)

 character(len=2), parameter :: period_elem(period_nelem) = &
 (/'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
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
   'Rg', 'Cn'/)

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
  write(iout,'(A)') "ERROR in subroutine read_ne_from_fch: no 'Number of elec'&
                   & found in file "//TRIM(fchname)
  stop
 end if

 read(buf(50:),*) ne
 return
end subroutine read_ne_from_fch

! check whether UHF-type MOs are hold in a given .fch(k) file
subroutine check_uhf_in_fch(fchname, uhf)
 implicit none
 integer :: i, fid
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
 logical, intent(out) :: uhf

 uhf = .false.
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit

  if(buf(1:7) == 'Beta MO') then
   uhf = .true.
   exit
  end if
 end do ! for while

 close(fid)
 return
end subroutine check_uhf_in_fch

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
 if(i == 1) then
  read(fid,'(A49,2X,I10)') buf, ncontr
  read(fid,'(A49,2X,I10)') buf, nprim
 else if(i == 2) then
  read(fid,'(A49,2X,I10)') buf, nprim
  read(fid,'(A49,2X,I10)') buf, ncontr
 else
  write(iout,'(A,I0)') 'ERROR in subroutine read_fch: invalid i=',i
  write(iout,'(A)') 'Unsupported file format, or file is incomplete. File='//TRIM(fchname)
  close(fid)
  stop
 end if

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
 allocate(eigen_e_a(nif), source=0d0)
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
  allocate(eigen_e_b(nif), source=0d0)
  read(fid,'(5(1X,ES15.8))') (eigen_e_b(i),i=1,nif)
 end if

 ! read Alpha MO
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:8) == 'Alpha MO') exit
 end do
 allocate(coeff(ncoeff), source=0d0)
 read(fid,'(5(1X,ES15.8))') (coeff(i),i=1,ncoeff)
 allocate(alpha_coeff(nbf,nif), source=0d0)
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

! read the position marks of 5D,6D,etc from array shell_type
subroutine read_mark_from_shltyp(sph, ncontr, shltyp, nd, nf, ng, nh, d_mark, f_mark, g_mark, h_mark)
 implicit none
 integer :: i, nbf
 integer, intent(in) :: ncontr
 integer, intent(in) :: shltyp(ncontr)
 integer, intent(out) :: nd, nf, ng, nh
 integer, intent(out) :: d_mark(ncontr), f_mark(ncontr), g_mark(ncontr), h_mark(ncontr)
 logical, intent(in) :: sph

 nbf = 0; nd = 0; nf = 0; ng = 0; nh = 0
 d_mark = 0; f_mark = 0; g_mark = 0; h_mark = 0
!     Spherical      |     Cartesian
! -6,-5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5, 6
!  I  H  G  F  D  L  S  P  D  F  G  H  I

 if(sph) then ! spherical harmonic functions
  do i = 1, ncontr, 1
   select case(shltyp(i))
   case( 0)   ! S
    nbf = nbf + 1
   case( 1)   ! 3P
    nbf = nbf + 3
   case(-1)   ! SP or L
    nbf = nbf + 4
   case(-2)   ! 5D
    nd = nd + 1
    d_mark(nd) = nbf + 1
    nbf = nbf + 5
   case(-3)   ! 7F
    nf = nf + 1
    f_mark(nf) = nbf + 1
    nbf = nbf + 7
   case(-4)   ! 9G
    ng = ng + 1
    g_mark(ng) = nbf + 1
    nbf = nbf + 9
   case(-5)   ! 11H
    nh = nh + 1
    h_mark(nh) = nbf + 1
    nbf = nbf + 11
   end select
  end do ! for i

 else ! Cartesian functions
  do i = 1, ncontr, 1
   select case(shltyp(i))
   case( 0)   ! S
    nbf = nbf + 1
   case( 1)   ! 3P
    nbf = nbf + 3
   case(-1)   ! SP or L
    nbf = nbf + 4
   case( 2)   ! 6D
    nd = nd + 1
    d_mark(nd) = nbf + 1
    nbf = nbf + 6
   case( 3)   ! 10F
    nf = nf + 1
    f_mark(nf) = nbf + 1
    nbf = nbf + 10
   case( 4)   ! 15G
    ng = ng + 1
    g_mark(ng) = nbf + 1
    nbf = nbf + 15
   case( 5)   ! 21H
    nh = nh + 1
    h_mark(nh) = nbf + 1
    nbf = nbf + 21
   end select
  end do ! for i
 end if

 return
end subroutine read_mark_from_shltyp

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

! expansion coefficients matrices of spherical harmonic -> Cartesian functions
module r_5D_2_6D
 implicit none
 real(kind=8), parameter :: r1 = 0.5d0, r2 = 0.5d0*DSQRT(3d0), r3 = 1.5d0/DSQRT(5d0)
 real(kind=8), parameter :: r4 = DSQRT(0.375d0), r5 = DSQRT(3d0/40d0), r6 = DSQRT(1.2d0)
 real(kind=8), parameter :: r7 = DSQRT(0.625d0), r8 = 3d0/DSQRT(8d0), r9 = DSQRT(27d0/35d0)
 real(kind=8), parameter :: r10= 3d0/8d0, r11 = 0.25d0*r9, r12 = DSQRT(10d0/7d0)
 real(kind=8), parameter :: r13= DSQRT(9d0/56d0), r14 = DSQRT(45d0/56d0), r15 = DSQRT(27d0/28d0)
 real(kind=8), parameter :: r16= 0.25*DSQRT(5d0), r17 = DSQRT(9d0/7d0), r18 = DSQRT(5d0/28d0)
 real(kind=8), parameter :: r19= 0.125d0*DSQRT(35d0), r20 = 0.75d0*DSQRT(3d0), r21 = 2d0*r16
 real(kind=8), parameter :: r22= DSQRT(25d0/21d0), r23 = 0.625d0, r24 = DSQRT(15d0/112d0)
 real(kind=8), parameter :: r25= DSQRT(5d0/3d0), r26 = DSQRT(9d0/28d0), r27 = 0.125d0*r25
 real(kind=8), parameter :: r28= DSQRT(45d0/28d0), r29= DSQRT(5d0/112d0), r30 = 0.125d0*DSQRT(15d0)
 real(kind=8), parameter :: r31= DSQRT(35d0/48d0), r32 = DSQRT(5d0/12d0), r33 = DSQRT(1.5d0)
 real(kind=8), parameter :: r34= DSQRT(35d0/128d0), r35 = DSQRT(5d0/6d0), r36 = 0.25d0*r35
 real(kind=8), parameter :: r37= DSQRT(175d0/128d0), r38 = 1.25d0*r33, r39 = DSQRT(63d0/128d0)

! This transformation table is taken from http://sobereva.com/97
 real(kind=8), parameter :: rd(6,5) = RESHAPE([-r1, -r1, 1d0, 0d0, 0d0, 0d0,&
                                               0d0, 0d0, 0d0, 0d0, 1d0, 0d0,&
                                               0d0, 0d0, 0d0, 0d0, 0d0, 1d0,&
                                                r2, -r2, 0d0, 0d0, 0d0, 0d0,&
                                               0d0, 0d0, 0d0, 1d0, 0d0, 0d0],[6,5])
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
rh(21,11) = RESHAPE([1d0, 0d0,-r22, 0d0, r23, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r22, 0d0, r24, 0d0, 0d0, 0d0, 0d0, r23, 0d0, 0d0,&
                     0d0, 0d0, 0d0, 0d0, 0d0, 0d0, r25, 0d0,-r26, 0d0, r27, 0d0, 0d0, 0d0, 0d0, r28, 0d0, r29, 0d0, 0d0, r30,&
                     0d0, r25, 0d0,-r26, 0d0, r30, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r26, 0d0, r29, 0d0, 0d0, 0d0, 0d0, r27, 0d0,&
                     0d0, 0d0,-r21, 0d0, r31, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, r21, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r31, 0d0, 0d0,&
                     0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, r25, 0d0,-r32, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r32, 0d0, 0d0, 0d0, 0d0,&
                     0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r33, 0d0, r34, 0d0, 0d0, 0d0, 0d0, r35, 0d0, r36, 0d0, 0d0,-r34,&
                     0d0, 0d0, 0d0,-r35, 0d0, r34, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, r33, 0d0,-r36, 0d0, 0d0, 0d0, 0d0,-r34, 0d0,&
                     0d0, 0d0, 0d0, 0d0, r19, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r20, 0d0, 0d0, 0d0, 0d0, r19, 0d0, 0d0,&
                     0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r21, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, r21, 0d0, 0d0, 0d0, 0d0,&
                     0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, r37, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r38, 0d0, 0d0, r39,&
                     0d0, 0d0, 0d0, 0d0, 0d0, r39, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,-r38, 0d0, 0d0, 0d0, 0d0, r37, 0d0],&
                    [21,11])
end module r_5D_2_6D

! check whether there exists DKH keywords in a given .fch(k) file
subroutine check_DKH_in_fch(fchname, order)
 implicit none
 integer :: i, fid
 integer, intent(out) :: order
! -2: no DKH
! -1: RESC
!  0: DKH 0th-order
!  2: DKH2
!  4: DKH4 with SO
 integer, parameter :: iout = 6
 character(len=61) :: buf
 character(len=240), intent(in) :: fchname
 character(len=610) :: longbuf
 logical :: alive(6)

 order = -2 ! default value, no DKH

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:5)=='Route' .or. buf(1:6)=='Charge') exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine check_DKH_in_fch: neither 'Route'&
                   & nor 'Charge' is found in file "//TRIM(fchname)//'.'
  write(iout,'(A)') 'This .fch(k) file is incomplete.'
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

 alive = [(index(longbuf,'DKH0')/=0), (index(longbuf,'DKHSO')/=0), &
          (index(longbuf,'DKH')/=0 .and. index(longbuf,'NODKH')==0), &
          (index(longbuf,'DKH2')/=0), &
          (index(longbuf,'DOUGLASKROLLHESS')/=0), (index(longbuf,'RESC')/=0)]

 if(alive(1)) then      ! DKH 0th order
  order = 0
 else if(alive(2)) then ! DKH 4th order
  order = 4
 else if(alive(6)) then ! RESC
  order = -1
 else if(alive(3) .or. alive(4) .or. alive(5)) then ! DKH 2nd order
  order = 2
 else
  order = -2           ! no relativistic
 end if

 return
end subroutine check_DKH_in_fch

! check whether there exists X2C keyword in a given .fch(k) file
! Note: normally the 'X2C' keyword is not possbible to be written in the
!  Gaussian .fch(k) file. This is merely because automr will write X2C into
!  the .fch(k) file.
subroutine check_X2C_in_fch(fchname, alive)
 implicit none
 integer :: i, fid
 integer, parameter :: iout = 6
 character(len=61) :: buf
 character(len=240), intent(in) :: fchname
 character(len=610) :: longbuf
 logical, intent(out) :: alive

 alive = .false. ! default value

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:5)=='Route' .or. buf(1:6)=='Charge') exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine check_X2C_in_fch: neither 'Route'&
                   & not 'Charge' is found in file "//TRIM(fchname)//'.'
  stop
 end if

 if(buf(1:5)=='Route') then
  longbuf = ' '
  do i = 1, 5
   read(fid,'(A)') buf
   if(buf(1:6) == 'Charge') exit
   longbuf = TRIM(longbuf)//TRIM(buf)
  end do ! for i

  call upper(longbuf)
  if(index(longbuf,'X2C') /= 0) alive = .true.
 end if

 close(fid)
 return
end subroutine check_X2C_in_fch

! check whether 'int=nobasistransform' exists in a given .fch(k) file
function nobasistransform_in_fch(fchname) result(notrans)
 implicit none
 integer:: i, fid
 character(len=240) :: buf
 character(len=1200) :: longbuf
 character(len=240), intent(in) :: fchname
 logical :: notrans

 notrans = .true.
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:5) == 'Route') exit

  if(buf(1:5) == 'Charg') then
   close(fid)
   return
  end if
 end do ! for while

 notrans = .false.

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
 if(index(longbuf,'nobasistransform') > 0) notrans = .true.

 return
end function nobasistransform_in_fch

! check whether 'nosymm' exists in a given .fch(k) file
function nosymm_in_fch(fchname) result(nosymm)
 implicit none
 integer:: i, fid
 character(len=240) :: buf
 character(len=1200) :: longbuf
 character(len=240), intent(in) :: fchname
 logical :: nosymm

 nosymm = .true.
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:5) == 'Route') exit

  if(buf(1:5) == 'Charg') then
   close(fid)
   return
  end if
 end do ! for while

 nosymm = .false.

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
 if(index(longbuf,'nosymm') > 0) nosymm = .true.

 return
end function nosymm_in_fch

