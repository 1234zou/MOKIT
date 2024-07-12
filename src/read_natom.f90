! moved from rwgeom.f90

subroutine read_natom_from_file(fname, natom)
 implicit none
 integer :: i
 integer, intent(out) :: natom
 character(len=240), intent(in) :: fname

 call require_file_exist(fname)
 i = LEN_TRIM(fname)

 select case(fname(i-3:i))
 case('.gjf')
  call read_natom_from_gjf(fname, natom)
 case('.xyz')
  call read_natom_from_xyz(fname, natom)
 case('.fch')
  call read_natom_from_fch(fname, natom)
 case('.pdb')
  call read_natom_from_pdb(fname, natom)
 case('grad')
  call read_natom_from_engrad(fname, natom)
 case('.gms')
  call read_natom_from_gms(fname, natom)
 case('.mkl')
  call read_natom_from_mkl(fname, natom)
 case('.EIn')
  call read_natom_from_EIn(fname, natom)
 case default
  write(6,'(/,A)') 'ERROR in subroutine read_natom_from_file: file format not s&
                   &upported.'
  write(6,'(A)') "Suffix='"//fname(i-3:i)//"'"
  stop
 end select
end subroutine read_natom_from_file

subroutine read_natom_from_gjf_pbc(gjfname, natom)
 implicit none
 integer, intent(out) :: natom
 character(len=240), intent(in) :: gjfname

 call read_natom_from_gjf(gjfname, natom)
 natom = natom - 3
end subroutine read_natom_from_gjf_pbc

! find the number of atoms in a .gjf file
subroutine read_natom_from_gjf(gjfname, natom)
 implicit none
 integer :: i, fid, nblank
 integer, intent(out) :: natom
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
  write(6,'(/,A)') 'ERROR in subroutine read_natom_from_gjf: incomplete file '&
                  //TRIM(gjfname)
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
end subroutine read_natom_from_gjf

! check whether the system is periodic in a given .gjf file
function check_pbc_in_gjf(gjfname) result(pbc)
 implicit none
 integer :: i, fid, nblank
 character(len=240) :: buf0, buf
 character(len=240), intent(in) :: gjfname
 logical :: pbc

 pbc = .false. ! initialization
 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
 nblank = 0
 buf0 = ' '

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) nblank = nblank + 1
  if(nblank == 3) exit
  buf0 = buf
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in function check_pbc_in_gjf: problematic file '//&
                   TRIM(gjfname)
  stop
 end if

 call upper(buf0(1:1))
 call lower(buf0(2:2))
 if(buf0(1:2) == 'Tv') pbc = .true.
end function check_pbc_in_gjf

! read the number of atoms from a given .fch file
subroutine read_natom_from_fch(fchname, natom)
 implicit none
 integer :: i, fid
 integer, intent(out) :: natom
!f2py intent(out) :: natom
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 natom = 0
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(buf(1:14) == 'Atomic numbers') exit
 end do ! for while

 BACKSPACE(fid)
 read(fid,'(A49,2X,I10)') buf, natom

 close(fid)
end subroutine read_natom_from_fch

! read the number of atoms from a .xyz file
subroutine read_natom_from_xyz(xyzname, natom)
 implicit none
 integer :: i, fid
 integer, intent(out) :: natom
 character(len=240), intent(in) :: xyzname

 natom = 0
 open(newunit=fid,file=TRIM(xyzname),status='old',position='rewind')
 read(fid,*,iostat=i) natom
 close(fid)

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_natom_from_xyz: failed to read nat&
                   &om from file '//TRIM(xyzname)
  stop
 end if
end subroutine read_natom_from_xyz

! read the number of atoms from a .pdb file
! Note: if there exists >1 frames, only the first frame will be dectected
subroutine read_natom_from_pdb(pdbname, natom)
 implicit none
 integer :: i, fid
 integer, intent(out) :: natom
 character(len=13) :: buf
 character(len=240), intent(in) :: pdbname

 natom = 0
 open(newunit=fid,file=TRIM(pdbname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:3) == 'END') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_natom_from_pdb: failed to read nat&
                   &om from file '//TRIM(pdbname)
  close(fid)
  stop
 end if

 BACKSPACE(fid)
 BACKSPACE(fid)
 read(fid,'(A)') buf
 if(buf(1:3) == 'TER') then
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
 end if
 close(fid)

 i = INDEX(buf, ' ')
 read(buf(i+1:),*) natom
end subroutine read_natom_from_pdb

! read the number of atoms from a (Open)Molcas output file
subroutine read_natom_from_molcas_out(outname, natom)
 implicit none
 integer :: i, fid
 integer, intent(out) :: natom
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
  write(6,'(/,A)') "ERROR in subroutine read_natom_from_molcas_out: keywords '+&
                   &+    Molecular struc'"
  write(6,'(A)') 'not found in file '//TRIM(outname)
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
end subroutine read_natom_from_molcas_out

! read the number of atoms from a Molpro output file
subroutine read_natom_from_molpro_out(outname, natom)
 implicit none
 integer :: i, fid
 integer, intent(out) :: natom
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
  write(6,'(/,A)') "ERROR in subroutine read_natom_from_molpro_out: no 'ATOMIC &
                   &COOR' found in"
  write(6,'(A)') 'file '//TRIM(outname)
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
end subroutine read_natom_from_molpro_out

! read the number of atoms from ORCA .engrad file
subroutine read_natom_from_engrad(fname, natom)
 implicit none
 integer :: fid
 integer, intent(out) :: natom
 character(len=240) :: buf
 character(len=240), intent(in) :: fname

 natom = 0
 open(newunit=fid,file=TRIM(fname),status='old',position='rewind')
 read(fid,'(A)') buf
 read(fid,'(A)') buf

 if(buf(1:16) /= '# Number of atom') then
  write(6,'(/,A)') "ERROR in subroutine read_natom_from_engrad: '# Number of at&
                   &om' expected, but"
  write(6,'(A)') "got '"//TRIM(buf)//"'"
  close(fid)
  stop
 end if

 read(fid,'(A)') buf
 read(fid,*) natom
 close(fid)
end subroutine read_natom_from_engrad

! read natom from GAMESS .gms file
subroutine read_natom_from_gms(gmsfile, natom)
 implicit none
 integer :: i, fid
 integer, intent(out) :: natom
 character(len=240) :: buf
 character(len=240), intent(in) :: gmsfile

 buf = ' '
 natom = 0
 open(newunit=fid,file=TRIM(gmsfile),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:22) == 'TOTAL NUMBER OF ATOMS') exit
 end do
 close(fid)

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_natom_from_gms: no 'TOTAL NUMBER O&
                   &F ATOMS' found"
  write(6,'(A)') 'in file '//TRIM(gmsfile)//'.'
  stop
 end if

 i = INDEX(buf,'=')
 read(buf(i+1:),*) natom
end subroutine read_natom_from_gms

! find the number of atoms in GAMESS .inp file
subroutine read_natom_from_gms_inp(inpname, natom)
 implicit none
 integer :: i, fid, nline
 integer, intent(out) :: natom
 character(len=1) :: str
 character(len=240):: buf
 character(len=240), intent(in) :: inpname

 natom = 0
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:2) == '$') then
   call upper(buf(3:6))
   if(buf(3:6) == 'DATA') exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine goto_data_section_in_gms_inp: wrong for&
                   &mat in file '//TRIM(inpname)
  close(fid)
  stop
 end if
 read(fid,'(A)') buf
 read(fid,'(A)') buf

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:2) == '$') then
   call upper(buf(3:5))
   if(buf(3:5) == 'END') exit
  end if

  do while(.true.)
   read(fid,'(A)') buf
   read(buf,*,iostat=i) str, nline
   if(i /= 0) exit

   do i = 1, nline, 1
    read(fid,'(A)') buf
   end do ! for i
  end do ! for while

  natom = natom + 1
 end do ! for while

 close(fid)
 if(natom == 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_natom_from_gms_inp: zero atom foun&
                   &d in file '//TRIM(inpname)
  stop
 end if
end subroutine read_natom_from_gms_inp

! read the number of atoms in .mkl file
subroutine read_natom_from_mkl(mklname, natom)
 implicit none
 integer :: i, fid
 integer, intent(out) :: natom
 character(len=240) :: buf
 character(len=240), intent(in) :: mklname

 natom = 0
 open(newunit=fid,file=TRIM(mklname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:5) == '$COOR') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_natom_from_mkl: no '$COOR' found i&
                   &n file "//TRIM(mklname)
  close(fid)
  stop
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:4) == '$END') exit
  natom = natom + 1
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_natom_from_mkl: no '$END' found in&
                   & file "//TRIM(mklname)
  stop
 end if
end subroutine read_natom_from_mkl

! read the number of atoms from a Gaussian .EIn file
subroutine read_natom_from_EIn(EIn, natom)
 implicit none
 integer :: fid
 integer, intent(out) :: natom
 character(len=240), intent(in) :: EIn

 open(newunit=fid,file=TRIM(EIn),status='old',position='rewind')
 read(fid,*) natom
 close(fid)
end subroutine read_natom_from_EIn

! read the number of atoms from an Amesp .amo file
subroutine read_natom_from_amo(amoname, natom)
 implicit none
 integer :: i, fid
 integer, intent(out) :: natom
 character(len=240) :: buf
 character(len=240), intent(in) :: amoname

 natom = 0
 open(newunit=fid,file=TRIM(amoname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:7) == '[Atoms]') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_natom_from_amo: '[Atoms]' section &
                   &not found in"
  write(6,'(A)') 'file '//TRIM(amoname)
  close(fid)
  stop
 end if

 read(fid,'(A)') buf ! 'Natm='
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:1) == '[') exit
  if(LEN_TRIM(buf) == 0) exit
  natom = natom + 1
 end do ! for while

 close(fid)
end subroutine read_natom_from_amo

! read the number of atoms from a .molden file
subroutine read_natom_from_molden(molden, natom)
 implicit none
 integer :: fid
 integer, intent(out) :: natom
 character(len=240) :: buf
 character(len=240), intent(in) :: molden

 open(newunit=fid,file=TRIM(molden),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:7)=='[Atoms]' .or. buf(2:8)=='[Atoms]' .or. buf(1:7)=='[ATOMS]') exit
 end do ! for while

 natom = 0
 do while(.true.)
  read(fid,'(A)') buf
  if(INDEX(buf(1:2),'[') > 0) exit
  natom = natom + 1
 end do ! for while

 close(fid)
end subroutine read_natom_from_molden

subroutine read_natom_from_gau_log(outname, natom)
 implicit none
 integer :: fid
 integer, intent(out) :: natom
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 natom = 0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:8) == 'NAtoms=') exit
 end do ! for while

 close(fid)
 read(buf(9:),*) natom
end subroutine read_natom_from_gau_log

! Read AtomTypes from a Dalton .mol file. For the fch2dal/mkl2dal generated .mol
! file, this number is equal to No. atoms.
subroutine read_natmtyp_from_dalton_mol(molname, natmtyp)
 implicit none
 integer :: i, fid
 integer, intent(out) :: natmtyp
 character(len=240) :: buf
 character(len=240), intent(in) :: molname

 natmtyp = 0
 open(newunit=fid,file=TRIM(molname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:10) == 'AtomTypes=') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_natmtyp_from_dalton_mol: no '^Atom&
                   &Types=' found in"
  write(6,'(A)') 'file '//TRIM(molname)
  stop
 end if

 read(buf(11:),*) natmtyp
end subroutine read_natmtyp_from_dalton_mol

! read elements and Cartesian coordinates from a Dalton .mol file
subroutine read_elem_and_coor_from_dalton_mol(molname, natom, elem, coor, nline)
 implicit none
 integer :: i, fid
 integer, intent(in) :: natom
 integer, intent(out) :: nline(natom) ! No. of lines of basis set data per atom
 real(kind=8), intent(out) :: coor(3,natom)
 character(len=2), intent(out) :: elem(natom)
 character(len=240) :: buf
 character(len=240), intent(in) :: molname

 elem = ' '; coor = 0d0; nline = 0
 open(newunit=fid,file=TRIM(molname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:10) == 'AtomTypes=') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_elem_and_coor_from_dalton_mol: no &
                   &'AtomTypes=' found"
  write(6,'(A)') 'in file '//TRIM(molname)
  close(fid)
  stop
 end if

 read(buf(11:),*) i
 if(i /= natom) then
  close(fid)
  write(6,'(/,A)') 'ERROR in subroutine read_elem_and_coor_from_dalton_mol: No.&
                   & atoms in .mol is'
  write(6,'(A)') 'not equal to input natom.'
  stop
 end if
 read(fid,'(A)') buf

 do i = 1, natom, 1
  read(fid,*) elem(i), coor(:,i)
  do while(.true.)
   read(fid,'(A)') buf
   if(i == natom) then
    if(LEN_TRIM(buf) == 0) exit
   end if
   if(buf(1:7) == 'Charge=') exit
   nline(i) = nline(i) + 1
  end do ! for while
 end do ! for i

 close(fid)
end subroutine read_elem_and_coor_from_dalton_mol

! read total charge and spin multiplicity from a given .gjf file
subroutine read_charge_and_mult_from_gjf(gjfname, charge, mult)
 implicit none
 integer :: i, nblank, fid
 integer, intent(out) :: charge, mult
 character(len=240) :: buf
 character(len=240), intent(in) :: gjfname

 charge = 0; mult = 1
 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
 nblank = 0

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) nblank = nblank + 1
  if(nblank == 2) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_charge_and_mult_from_gjf: problema&
                   &tic file '//TRIM(gjfname)
  stop
  close(fid)
 end if

 read(fid,*,iostat=i) charge, mult
 close(fid)

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_charge_and_mult_from_gjf: problema&
                   &tic charge or'
  write(6,'(A)') 'spin multiplicity in file '//TRIM(gjfname)
  stop
 end if
end subroutine read_charge_and_mult_from_gjf

