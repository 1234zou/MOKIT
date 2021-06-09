! written by jxzou at 20200623

module periodic_table
 implicit none
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
 ! the radii below are read from Gaussian output file
 real(kind=8), parameter :: vdw_radii(period_nelem) = &
  [1.4430d0, 1.1810d0, 1.2255d0, 1.3725d0, 2.0415d0, &
   1.9255d0, 1.8300d0, 1.7500d0, 1.6820d0, 1.6215d0, &
   1.4915d0, 1.5105d0, 2.2495d0, 2.1475d0, 2.0735d0, &
   2.0175d0, 1.9735d0, 1.9340d0, 1.9060d0, 1.6995d0, &
   1.6475d0, 1.5875d0, 1.5720d0, 1.5115d0, 1.4805d0, &
   1.4560d0, 1.4360d0, 1.4170d0, 1.7475d0, 1.3815d0, &
   2.1915d0, 2.1400d0, 2.1150d0, 2.1025d0, 2.0945d0, &
   2.0705d0, 2.0570d0, 1.8205d0, 1.6725d0, 1.5620d0, &
   1.5825d0, 1.5260d0, 1.4990d0, 1.4815d0, 1.4645d0, &
   1.4495d0, 1.5740d0, 1.4240d0, 2.2315d0, 2.1960d0, &
   2.2100d0, 2.2350d0, 2.2500d0, 2.2020d0, 2.2585d0, &
   1.8515d0, 1.7610d0, 1.7780d0, 1.8030d0, 1.7875d0, &
   1.7735d0, 1.7600d0, 1.7465d0, 1.6840d0, 1.7255d0, &
   1.7140d0, 1.7045d0, 1.6955d0, 1.6870d0, 1.6775d0, &
   1.8200d0, 1.5705d0, 1.5850d0, 1.5345d0, 1.4770d0, &
   1.5600d0, 1.4200d0, 1.3770d0, 1.6465d0, 1.3525d0, &
   2.1735d0, 2.1485d0, 2.1850d0, 2.3545d0, 2.3750d0, &
   2.3825d0, 2.4500d0, 1.8385d0, 1.7390d0, 1.6980d0, &
   1.7120d0, 1.6975d0, 1.7120d0, 1.7120d0, 1.6905d0, &
   1.6630d0, 1.6695d0, 1.6565d0, 1.6495d0, 1.6430d0, &
   1.6370d0, 1.6240d0, 1.6180d0, 1.7500d0, 1.7500d0, &
   1.7500d0, 1.7500d0, 1.7500d0, 1.7500d0, 1.7500d0, &
   1.7500d0, 1.7500d0]
contains

 ! map a nuclear charge to an element (e.g. 6->'C')
 ! Note: only 1-112 elements are supported!
 pure function nuc2elem(i) result(s)
  implicit none
  integer, intent(in) :: i
  character(len=2) :: s
 
  s = period_elem(i)
  return
 end function nuc2elem

 pure function elem2vdw_radii(elem) result(radii)
  implicit none
  real(kind=8) :: radii
  character(len=2), intent(in) :: elem

  radii = vdw_radii(elem2nuc(elem))
  return
 end function elem2vdw_radii

 ! map an element to a nuclear charge (e.g. 'C'->6)
 ! Note: only 1-112 elements are supported!
 pure function elem2nuc(s) result(i)
  implicit none
  integer :: i
  character(len=2), intent(in) :: s

  do i = 1, period_nelem, 1
   if(period_elem(i) == s) return
  end do ! for i
  return
 end function elem2nuc

 ! read elements array from a given .gjf array
 subroutine read_elem_from_gjf(gjfname, natom, elem, ghost)
  implicit none
  integer :: i, j, k, nblank, fid
  integer, intent(in) :: natom
  character(len=6) :: str
  character(len=240) :: buf
  character(len=240), intent(in) :: gjfname
  character(len=2), intent(out) :: elem(natom)
  logical, intent(out) :: ghost(natom)

  elem = ' '
  ghost = .false.
  open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
  nblank = 0

  do while(.true.)
   read(fid,'(A)') buf
   if(LEN_TRIM(buf) == 0) nblank = nblank + 1
   if(nblank == 2) exit
  end do ! for while

  read(fid,'(A)') buf ! skip charge and mult
  do i = 1, natom, 1
   read(fid,*) str
   str = ADJUSTL(str)
   k = LEN_TRIM(str)
   if(str(k-2:k) == '-Bq') then
    elem(i) = str(1:k-3)
    ghost(i) = .true.
   else if(str(2:3)=='(f' .or. str(2:3)=='(F') then
    call upper(str(1:1))
    elem(i) = str(1:1)//' '
   else if(str(3:4)=='(f' .or. str(3:4)=='(F') then
    call upper(str(1:1))
    call lower(str(2:2))
    elem(i) = str(1:2)
   else
    k = IACHAR(str(1:1))
    if((k>64 .and. k<91) .or. (k>96 .and. k<123)) then ! A-Z, a-z
     if(k>96 .and. k<123) str(1:1) = ACHAR(k-32)
     j = IACHAR(str(2:2))
     if(j>64 .and. j<91) str(2:2) = ACHAR(j+32)
     elem(i) = TRIM(str)
    else
     read(str,*) j
     elem(i) = nuc2elem(j)
    end if
   end if
  end do ! for natom

  close(fid)
  return
 end subroutine read_elem_from_gjf
end module periodic_table

! find the number of atoms in Gaussian .gjf file
subroutine read_natom_from_gjf(gjfname, natom)
 implicit none
 integer :: i, fid, nblank
 integer, intent(out) :: natom
 integer, parameter :: iout = 6
 character(len=240) :: buf
 character(len=240), intent(in) :: gjfname

 nblank = 0
 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) nblank = nblank + 1
  if(nblank == 2) exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine read_natom_from_gjf: incomplete file '//TRIM(gjfname)
  stop
 end if

 read(fid,'(A)') buf ! skip charge and mult

 natom = 0
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) exit
  natom = natom + 1
 end do ! for while

 close(fid)
 return
end subroutine read_natom_from_gjf

! read the number of atoms from a given .fch file
subroutine read_natom_from_fch(fchname, natom)
 implicit none
 integer :: i, fid
 integer, intent(out) :: natom
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname

 natom = 0
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(buf(1:14) == 'Atomic numbers') exit
 end do ! for while

 BACKSPACE(fid)
 read(fid,'(A49,2X,I10)') buf, natom

 close(fid)
 return
end subroutine read_natom_from_fch

! read the number of atoms from a .xyz file
subroutine read_natom_from_xyz(xyzname, natom)
 implicit none
 integer :: i, fid
 integer, intent(out) :: natom
 integer, parameter :: iout = 6
 character(len=240), intent(in) :: xyzname

 natom = 0
 open(newunit=fid,file=TRIM(xyzname),status='old',position='rewind')
 read(fid,*,iostat=i) natom
 close(fid)

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine read_natom_from_xyz: failed to read&
                   & natom from file '//TRIM(xyzname)
  stop
 end if
 return
end subroutine read_natom_from_xyz

! read the number of atoms from a (Open)Molcas output file
subroutine read_natom_from_molcas_out(outname, natom)
 implicit none
 integer :: i, fid
 integer, intent(out) :: natom
 integer, parameter :: iout = 6
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 natom = 0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:21) == '++    Molecular struc') exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine read_natom_from_molcas_out: keywords&
                   & '++    Molecular struc' not found"
  write(iout,'(A)') 'in file '//TRIM(outname)
  close(fid)
  stop
 end if

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(6:11) == 'Center') exit 
 end do ! for while

 do while(.true.)
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
  natom = natom + 1 
 end do ! for while

 close(fid)
 return
end subroutine read_natom_from_molcas_out

! read the number of atoms from a Molpro output file
subroutine read_natom_from_molpro_out(outname, natom)
 implicit none
 integer :: i, fid
 integer, intent(out) :: natom
 integer, parameter :: iout = 6
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 natom = 0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:12) == 'ATOMIC COOR') exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine read_natom_from_molpro_out: no '&
                   &ATOMIC COOR' found in file "//TRIM(outname)
  close(fid)
  stop
 end if

 do i = 1, 3
  read(fid,'(A)') buf
 end do

 do while(.true.)
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
  natom = natom + 1
 end do ! for while

 close(fid)
 return
end subroutine read_natom_from_molpro_out

! read 3 arrays elem, nuc, coor, and the total charge as well as multiplicity
! from a given .gjf file
subroutine read_elem_and_coor_from_gjf(gjfname, natom, elem, nuc, coor, charge, mult)
 use fch_content, only: elem2nuc
 implicit none
 integer :: i, j, k, fid, nblank, ne
 integer, intent(in) :: natom
 integer, intent(out) :: charge, mult, nuc(natom)
 integer, parameter :: iout = 6
 real(kind=8), intent(out) :: coor(3,natom)
 real(kind=8), parameter :: Bohr_const = 0.52917721092d0
 character(len=2), intent(out) :: elem(natom)
 character(len=240) :: buf
 character(len=240), intent(in) :: gjfname
 logical :: bohr

 charge = 0; mult = 1; coor = 0d0; bohr = .false.
 nblank = 0

 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) nblank = nblank + 1
  if(nblank == 2) exit

  if(buf(1:1) =='#') then
   if(index(buf,'units=au') /= 0) bohr = .true.
  end if
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine read_elem_and_coor_from_gjf: incomplete file '//TRIM(gjfname)
  stop
 end if

 read(fid,*,iostat=k) charge, mult
 if(k /= 0) then
  write(iout,'(A)') 'ERROR in subroutine read_elem_and_coor_from_gjf: failed&
                   & to read charge and mult.'
  write(iout,'(A)') 'There exists syntax error in file '//TRIM(gjfname)
  stop
 end if

 do i = 1, natom, 1
  read(fid,'(A)') buf
  j = index(buf,'('); k = index(buf,')') ! in case for the fragment guess wfn
  if(j*k /= 0) buf(j:k) = ' '

  read(buf,*,iostat=k) elem(i), coor(1:3,i)
  if(k /= 0) then
   write(iout,'(A)') 'ERROR in subroutine read_elem_and_coor_from_gjf: only 4-column&
                    & format is supported.'
   stop
  end if
 end do ! for i

 close(fid)

 if(bohr) coor = coor*Bohr_const ! convert Bohr to Angstrom

 forall(i=1:natom)
  elem(i) = ADJUSTL(elem(i))
  nuc(i) = elem2nuc(elem(i))
 end forall

 ne = SUM(nuc) - charge
 if(MOD(ne,2) /= MOD(mult-1,2)) then
  write(iout,'(/,A)') 'ERROR in subroutine read_elem_and_coor_from_gjf:'
  write(iout,'(2(A,I0),A)') 'The combination of multiplicity ',mult,' and ',&
                             ne,' electrons is impossible.'
  stop
 end if

 return
end subroutine read_elem_and_coor_from_gjf

! read charge, spin multiplicities and atom2frag from a given .gjf file
subroutine read_frag_guess_from_gjf(gjfname, natom, atom2frag, nfrag, frag_char_mult)
 implicit none
 integer :: i, j, k, nblank, charge, mult, fid
 integer, parameter :: iout = 6
 integer, intent(in) :: natom, nfrag
 integer, intent(out) :: atom2frag(natom), frag_char_mult(2,nfrag)
 character(len=240) :: buf
 character(len=240), intent(in) :: gjfname

 atom2frag = 0; frag_char_mult = 0
 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')

 nblank = 0
 do while(.true.)
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) nblank = nblank + 1
  if(nblank == 2) exit
 end do ! for while

 read(fid,*,iostat=i) charge, mult, ((frag_char_mult(j,i),j=1,2),i=1,nfrag)
 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine read_frag_guess_from_gjf: failed to read&
                   & charges and spin multiplicities of fragments.'
  write(iout,'(A)') 'Please check syntax in file '//TRIM(gjfname)
  close(fid)
  stop
 end if

 do i = 1, natom, 1
  read(fid,'(A)') buf
  j = index(buf,'='); k = index(buf,')')

  if(j*k == 0) then
   write(iout,'(A)') 'ERROR in subroutine read_frag_guess_from_gjf: failed to read&
                     & atom2frag.'
   write(iout,'(A)') 'Problematic line: '//TRIM(buf)
   write(iout,'(A)') 'Please check syntax in file '//TRIM(gjfname)
   close(fid)
   stop
  end if

  read(buf(j+1:k-1),*) atom2frag(i)
 end do ! for i

 close(fid)
 return
end subroutine read_frag_guess_from_gjf

! read 3 arrays elem, nuc, coor, and the total charge as well as multiplicity
! from a given .fch file
subroutine read_elem_and_coor_from_fch(fchname, natom, elem, nuc, coor, charge, mult)
 use fch_content, only: nuc2elem
 implicit none
 integer :: i, fid
 integer, intent(in) :: natom
 integer, intent(out) :: charge, mult, nuc(natom)
 integer, parameter :: iout = 6
 real(kind=8), intent(out) :: coor(3,natom)
 real(kind=8), allocatable :: coor0(:)
 real(kind=8), parameter :: Bohr_const = 0.52917721092d0
 character(len=2), intent(out) :: elem(natom)
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 ! find charge and mult
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:6) == 'Charge') exit
 end do
 BACKSPACE(fid)
 read(fid,'(A49,2X,I10)') buf, charge
 read(fid,'(A49,2X,I10)') buf, mult

 ! find atomic numbers/nuclear charges
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:14) == 'Atomic numbers') exit
 end do ! for i
 read(fid,'(6(1X,I11))') (nuc(i),i=1,natom)

 ! find and read coordinates
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:12) == 'Current cart') exit
 end do
 allocate(coor0(3*natom), source=0.0d0)
 read(fid,'(5(1X,ES15.8))') (coor0(i),i=1,3*natom)
 coor = RESHAPE(coor0,[3,natom])

 deallocate(coor0)
 close(fid)

 coor = coor*Bohr_const ! convert Bohr to Angstrom
 forall(i=1:natom) elem(i) = nuc2elem(nuc(i))
 return
end subroutine read_elem_and_coor_from_fch

! read Cartesian gradient from a given PySCF output file
subroutine read_grad_from_pyscf_out(outname, natom, grad)
 implicit none
 integer :: i, k, fid
 integer, intent(in) :: natom
 real(kind=8), intent(out) :: grad(3*natom)
 character(len=3) :: elem
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 grad = 0.0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(index(buf,'gradients') /= 0) exit
 end do ! for while

 read(fid,'(A)') buf

 do i = 1, natom, 1
  read(fid,*) k, elem, grad(3*i-2:3*i)
 end do ! for i

 close(fid)
 return
end subroutine read_grad_from_pyscf_out

! read Cartesian gradient from a given .fch file
subroutine read_grad_from_fch(fchname, natom, grad)
 implicit none
 integer :: i, fid
 integer, intent(in) :: natom
 real(kind=8), intent(out) :: grad(3*natom)
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname

 grad = 0.0d0
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:26) == 'Opt point       1 Gradient') exit
 end do ! for while

 read(fid,'(5(1X,ES15.8))') (grad(i),i=1,3*natom)
 close(fid)
 return
end subroutine read_grad_from_fch

! read Cartesian gradient from a given Gaussian .log file
subroutine read_grad_from_gau_log(logname, natom, grad)
 implicit none
 integer :: i, k1, k2, fid
 integer, intent(in) :: natom
 real(kind=8), intent(out) :: grad(3*natom)
 character(len=240) :: buf
 character(len=240), intent(in) :: logname

 grad = 0.0d0
 open(newunit=fid,file=TRIM(logname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(38:52) == 'Forces (Hartree') exit
 end do ! for while

 read(fid,'(A)') buf
 read(fid,'(A)') buf

 do i = 1, natom, 1
  read(fid,*) k1, k2, grad(3*i-2:3*i)
 end do ! for i

 close(fid)

 grad = -grad !!! VIP
 return
end subroutine read_grad_from_gau_log

! read Cartesian gradient from a given GAMESS .gms file
subroutine read_grad_from_gms_gms(outname, natom, grad)
 implicit none
 integer :: i, k, fid
 integer, intent(in) :: natom
 real(kind=8), intent(out) :: grad(3*natom)
 character(len=2) :: elem = ' '
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 grad = 0.0d0

 open(newunit=fid,file=TRIM(outname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(buf(26:47) == 'GRADIENT OF THE ENERGY') exit
 end do ! for while

 do i = 1, 3
  read(fid,'(A)') buf
 end do ! for i

 do i = 1, natom, 1
  read(fid,*) k, elem, grad(3*i-2:3*i)
 end do ! for i

 close(fid)
 return
end subroutine read_grad_from_gms_gms

! read Cartesian gradient from a given GAMESS .dat file
subroutine read_grad_from_gms_dat(datname, natom, grad)
 implicit none
 integer :: i, fid
 integer, intent(in) :: natom
 real(kind=4) :: r
 real(kind=8), intent(out) :: grad(3*natom)
 character(len=2) :: elem
 character(len=240) :: buf
 character(len=240), intent(in) :: datname

 grad = 0.0d0

 open(newunit=fid,file=TRIM(datname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(buf(2:6) == '$GRAD') exit
 end do ! for while

 read(fid,'(A)') buf
 do i = 1, natom, 1
  read(fid,*) elem, r, grad(3*i-2:3*i)
 end do ! for i

 close(fid)
 return
end subroutine read_grad_from_gms_dat

! read Cartesian gradient from a given MOLCAS/OpenMolcas .out file
subroutine read_grad_from_molcas_out(outname, natom, grad)
 implicit none
 integer :: i, fid
 integer, intent(in) :: natom
 real(kind=8), intent(out) :: grad(3*natom)
 character(len=2) :: elem = ' '
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 grad = 0.0d0

 open(newunit=fid,file=TRIM(outname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(buf(17:35) == 'Molecular gradients') exit
 end do ! for while

 do i = 1, 7, 1
  read(fid,'(A)') buf
 end do ! for i

 do i = 1, natom, 1
  read(fid,*) elem, grad(3*i-2:3*i)
 end do ! for i

 close(fid)
 return
end subroutine read_grad_from_molcas_out

! read Cartesian gradient from a given ORCA .out file
subroutine read_grad_from_orca_out(outname, natom, grad)
 implicit none
 integer :: i, k, fid
 integer, intent(in) :: natom
 real(kind=8), intent(out) :: grad(3*natom)
 character(len=1) :: str = ' '
 character(len=2) :: elem = ' '
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 grad = 0.0d0

 open(newunit=fid,file=TRIM(outname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(buf(1:14) == 'CARTESIAN GRAD') exit
 end do ! for while

 read(fid,'(A)') buf
 read(fid,'(A)') buf

 do i = 1, natom, 1
  read(fid,*) k, elem, str, grad(3*i-2:3*i)
 end do ! for i

 close(fid)
 return
end subroutine read_grad_from_orca_out

! read Cartesian gradient from a given Molpro .out file
subroutine read_grad_from_molpro_out(outname, natom, grad)
 implicit none
 integer :: i, k, fid
 integer, intent(in) :: natom
 real(kind=8), intent(out) :: grad(3*natom)
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 grad = 0.0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(buf(2:8) == 'MC GRAD') exit
 end do ! for while

 do i = 1, 3
  read(fid,'(A)') buf
 end do

 do i = 1, natom, 1
  read(fid,*) k, grad(3*i-2:3*i)
 end do ! for i

 close(fid)
 return
end subroutine read_grad_from_molpro_out

! read Cartesian gradient from a given BDF .out file
subroutine read_grad_from_bdf_out(outname, natom, grad)
 implicit none
 integer :: i, k, fid
 integer, intent(in) :: natom
 real(kind=8), intent(out) :: grad(3*natom)
 character(len=10) :: str
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 grad = 0.0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(index(buf,'Molecular gradient - Mol') /= 0) exit
 end do ! for while

 read(fid,'(A)') buf
 read(fid,'(A)') buf

 do i = 1, natom, 1
  read(fid,*) str, grad(3*i-2:3*i)
 end do ! for i

 close(fid)
 return
end subroutine read_grad_from_bdf_out

! read Cartesian xyz coordinates from a .xyz file
! Note: 1) return array coor(3,natom) are in unit Angstrom
!       2) if 'bohr' key is found in the 2nd line of the xyz file,
!          coor will be 
subroutine read_coor_from_xyz(xyzname, natom, coor)
 implicit none
 integer :: i, k, fid
 integer, intent(in) :: natom
 integer, parameter :: iout = 6
 real(kind=8), intent(out) :: coor(3,natom)
 real(kind=8), parameter :: Bohr_const = 0.52917721092d0
 character(len=3) :: elem
 character(len=240) :: buf
 character(len=240), intent(in) :: xyzname
 logical :: bohr

 coor = 0d0
 open(newunit=fid,file=TRIM(xyzname),status='old',position='rewind')
 read(fid,'(A)') buf
 read(fid,'(A)') buf

 bohr = .false.
 call lower(buf)
 if(index(buf,'bohr') > 0) then
  if(index(buf,'angstrom') > 0) then
   write(iout,'(A)') "ERROR in subroutine read_coor_from_xyz: it's confusing&
                    & because both 'bohr' and 'angstrom'"
   write(iout,'(A)') 'are detected in the 2nd line of file '//TRIM(xyzname)
   close(fid)
   stop
  else
   bohr = .true.
  end if
 end if

 do i = 1, natom, 1
  read(fid,*,iostat=k) elem, coor(1:3,i)
  if(k /= 0) then
   write(iout,'(A)') 'ERROR in subroutine read_coor_from_xyz: insufficient&
                    & number of atoms in file '//TRIM(xyzname)
   write(iout,'(2(A,I0))') 'Input natom=', natom, ', but broken at i=', i
   close(fid)
   stop
  end if
 end do ! for i

 close(fid)
 if(bohr) coor = coor*Bohr_const ! convert Bohr to Angstrom
 return
end subroutine read_coor_from_xyz

! read Cartesian coordinates from a (Open)Molcas output file
subroutine read_coor_from_molcas_out(outname, natom, coor)
 implicit none
 integer :: i, j, k, fid
 integer, intent(in) :: natom
 integer, parameter :: iout = 6
 real(kind=8), intent(out) :: coor(3,natom)
 character(len=6) :: str
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 buf = ' '
 coor = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(buf(4:21) == 'This run of MOLCAS') then
   write(iout,'(A)') "ERROR in subroutine read_coor_from_molcas_out: failed to&
                    & find 'Cartesian coordinates in A'"
   write(iout,'(A)') 'keywords in file '//TRIM(outname)
   close(fid)
   stop
  end if
  if(buf(7:32) == 'Cartesian coordinates in A') exit
 end do ! for while

 do i = 1, 3
  read(fid,'(A)') buf
 end do

 do i = 1, natom, 1
  read(fid,*,iostat=j) k, str, coor(1:3,i)
  if(j /= 0) then
   write(iout,'(A)') 'ERROR in subroutine read_coor_from_molcas_out: insufficient&
                    & number of atoms in file '//TRIM(outname)
   write(iout,'(2(A,I0))') 'natom=', natom, ', but broken at i=', i
   close(fid)
   stop
  end if
 end do ! for i

 close(fid)
 return
end subroutine read_coor_from_molcas_out

! read Cartesian coordinates from a Molpro output file
subroutine read_coor_from_molpro_out(outname, natom, coor)
 implicit none
 integer :: i, fid
 integer, intent(in) :: natom
 integer, parameter :: iout = 6
 real(kind=8), intent(out) :: coor(3,natom)
 character(len=2) :: elem = ' '
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 coor = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(buf(2:13) == 'Current geom') exit

  if(buf(2:13) == 'Primary work') then
   write(iout,'(A)') "ERROR in subroutine read_coor_from_molpro_out: no '&
                   &Current geom' found in file "//TRIM(outname)
   close(fid)
   stop
  end if
 end do ! for while

 do i = 1, 3
  read(fid,'(A)') buf
 end do

 do i = 1, natom, 1
  read(fid,*) elem, coor(1:3,i)
 end do ! for i

 close(fid)
 return
end subroutine read_coor_from_molpro_out

! write/create a .xyz file
subroutine write_xyzfile(natom, coor, elem, xyzname)
 implicit none
 integer :: i, fid
 integer, intent(in) :: natom
 real(kind=8), intent(in) :: coor(3,natom)
 character(len=2), intent(in) :: elem(natom)
 character(len=240), intent(in) :: xyzname

 open(newunit=fid,file=TRIM(xyzname),status='replace')
 write(fid,'(I0)') natom
 write(fid,'(A)') 'xyz format file produced by MOKIT'

 do i = 1, natom, 1
  write(fid,'(A2,3(1X,F18.8))') elem(i), coor(:,i)
 end do ! for i

 close(fid)
 return
end subroutine write_xyzfile

