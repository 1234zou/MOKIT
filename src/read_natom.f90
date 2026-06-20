! moved from rwgeom.f90

subroutine read_natom_from_file(fname, natom)
 implicit none
 integer :: i
 integer, intent(out) :: natom
!f2py intent(out) :: natom
 character(len=240), intent(in) :: fname
!f2py intent(in) :: fname

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
!f2py intent(out) :: natom
 character(len=240), intent(in) :: gjfname
!f2py intent(in) :: gjfname

 call read_natom_from_gjf(gjfname, natom)
 natom = natom - 3
end subroutine read_natom_from_gjf_pbc

! find the number of atoms in a .gjf file
subroutine read_natom_from_gjf(gjfname, natom)
 implicit none
 integer :: i, fid, nblank
 integer, intent(out) :: natom
!f2py intent(out) :: natom
 character(len=240) :: buf
 character(len=240), intent(in) :: gjfname
!f2py intent(in) :: gjfname

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
subroutine check_pbc_in_gjf(gjfname, pbc)
 implicit none
 integer :: i, fid, nblank
 character(len=240) :: buf0, buf
 character(len=240), intent(in) :: gjfname
!f2py intent(in) :: gjfname
 logical, intent(out) :: pbc
!f2py intent(out) :: pbc

 nblank = 0; buf = ' '; pbc = .false.
 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) nblank = nblank + 1
  if(nblank == 3) exit
  buf0 = buf
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine check_pbc_in_gjf: problematic file '//&
                   TRIM(gjfname)
  stop
 end if

 buf0 = ADJUSTL(buf0)
 call upper(buf0(1:1))
 call lower(buf0(2:2))
 if(buf0(1:2) == 'Tv') pbc = .true.
end subroutine check_pbc_in_gjf

! check whether the system is periodic in a given .xyz file
subroutine check_pbc_in_xyz(xyzname, pbc)
 implicit none
 integer :: fid
 character(len=240) :: buf
 character(len=240), intent(in) :: xyzname
!f2py intent(in) :: xyzname
 logical, intent(out) :: pbc
!f2py intent(out) :: pbc

 pbc = .false.
 open(newunit=fid,file=TRIM(xyzname),status='old',position='rewind')
 read(fid,'(A)') buf
 read(fid,'(A)') buf
 close(fid)

 buf = ADJUSTL(buf)
 if(LEN_TRIM(buf) > 6) then
  call upper(buf(1:1))
  call lower(buf(2:7))
  if(buf(1:7) == 'Lattice') pbc = .true.
 end if
end subroutine check_pbc_in_xyz

! check whether the system is periodic in a specified .pdb file
subroutine check_pbc_in_pdb(pdbname, pbc)
 implicit none
 integer :: fid
 character(len=240) :: buf
 character(len=240), intent(in) :: pdbname
!f2py intent(in) :: pdbname
 logical, intent(out) :: pbc
!f2py intent(out) :: pbc

 pbc = .false.
 open(newunit=fid,file=TRIM(pdbname),status='old',position='rewind')
 read(fid,'(A)') buf
 read(fid,'(A)') buf
 close(fid)
 if(buf(1:6) == 'CRYST1') pbc = .true.
end subroutine check_pbc_in_pdb

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
!f2py intent(out) :: natom
 character(len=240), intent(in) :: xyzname
!f2py intent(in) :: xyzname

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
! Note: if there exists >1 frames, only the first frame will be detected
subroutine read_natom_from_pdb(pdbname, natom)
 implicit none
 integer :: i, fid
 integer, intent(out) :: natom
!f2py intent(out) :: natom
 character(len=13) :: buf
 character(len=240), intent(in) :: pdbname
!f2py intent(in) :: pdbname

 natom = 0
 open(newunit=fid,file=TRIM(pdbname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:3) == 'END') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_natom_from_pdb: failed to read nat&
                   &om from file'
  write(6,'(A)') TRIM(pdbname)
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
!f2py intent(out) :: natom
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
!f2py intent(in) :: outname

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
!f2py intent(out) :: natom
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
!f2py intent(in) :: outname

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
!f2py intent(out) :: natom
 character(len=240) :: buf
 character(len=240), intent(in) :: fname
!f2py intent(in) :: fname

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
!f2py intent(out) :: natom
 character(len=240) :: buf
 character(len=240), intent(in) :: gmsfile
!f2py intent(in) :: gmsfile

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
!f2py intent(out) :: natom
 character(len=1) :: str
 character(len=240):: buf
 character(len=240), intent(in) :: inpname
!f2py intent(in) :: inpname

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
!f2py intent(out) :: natom
 character(len=240) :: buf
 character(len=240), intent(in) :: mklname
!f2py intent(in) :: mklname

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
!f2py intent(out) :: natom
 character(len=240), intent(in) :: EIn
!f2py intent(in) :: EIn

 open(newunit=fid,file=TRIM(EIn),status='old',position='rewind')
 read(fid,*) natom
 close(fid)
end subroutine read_natom_from_EIn

! read the number of atoms from an Amesp .amo file
subroutine read_natom_from_amo(amoname, natom)
 implicit none
 integer :: i, fid
 integer, intent(out) :: natom
!f2py intent(out) :: natom
 character(len=240) :: buf
 character(len=240), intent(in) :: amoname
!f2py intent(in) :: amoname

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
!f2py intent(out) :: natom
 character(len=240) :: buf
 character(len=240), intent(in) :: molden
!f2py intent(in) :: molden

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
 integer :: i, fid
 integer, intent(out) :: natom
!f2py intent(out) :: natom
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
!f2py intent(in) :: outname

 natom = 0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:8) == 'NAtoms=') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_natom_from_gau_log: "NAtoms=" not &
                   &located in file'
  write(6,'(A)') TRIM(outname)
  stop
 end if

 read(buf(9:),*) natom
end subroutine read_natom_from_gau_log

subroutine read_natom_from_orca_out(outname, natom)
 implicit none
 integer :: i, fid
 integer, intent(out) :: natom
!f2py intent(out) :: natom
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
!f2py intent(in) :: outname

 natom = 0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:15) == 'Number of atoms') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_natom_from_orca_out: "Number of at&
                   &oms" not located in'
  write(6,'(A)') 'file '//TRIM(outname)
  stop
 end if

 i = INDEX(buf, '.', back=.true.)
 read(buf(i+1:),*) natom
end subroutine read_natom_from_orca_out

subroutine read_natom_from_orca_prop_txt(txtname, natom)
 implicit none
 integer :: i, fid
 integer, intent(out) :: natom
!f2py intent(out) :: natom
 character(len=240) :: buf
 character(len=240), intent(in) :: txtname
!f2py intent(in) :: txtname

 natom = 0
 open(newunit=fid,file=TRIM(txtname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(5:10) == 'NAtoms') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_natom_from_orca_prop_txt: 'NAtoms'&
                   & not found in"
  write(6,'(A)') 'file '//TRIM(txtname)
  stop
 end if

 i = INDEX(buf, ']')
 read(buf(i+1:),*) natom
end subroutine read_natom_from_orca_prop_txt

subroutine read_natom_from_cp2k_inp(inpname, natom)
 implicit none
 integer :: i, k, fid
 integer, intent(out) :: natom
!f2py intent(out) :: natom
 character(len=46), parameter :: error_warn = 'ERROR in subroutine read_natom_f&
                                              &rom_cp2k_inp: '
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname
!f2py intent(in) :: inpname

 natom = 0
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  buf = ADJUSTL(buf)
  k = LEN_TRIM(buf)
  call upper(buf(1:k))
  if(buf(1:6) == '&COORD') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'failed to locate "&COORD" in file'
  write(6,'(A)') TRIM(inpname)
  close(fid)
  stop
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  buf = ADJUSTL(buf)
  k = LEN_TRIM(buf)
  call upper(buf(1:k))
  if(buf(1:10) == '&END COORD') exit
  natom = natom + 1
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') error_warn//'failed to locate "&COORD" in file'
  write(6,'(A)') TRIM(inpname)
  stop
 end if
end subroutine read_natom_from_cp2k_inp

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

! read lattice vectors (3,3) from a specified .molden file
subroutine read_lat_vec_from_molden(molden, lat_vec)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: lat_vec(3,3)
 character(len=6) :: str6
 character(len=240) :: buf
 character(len=240), intent(in) :: molden

 lat_vec = 0d0
 open(newunit=fid,file=TRIM(molden),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  buf = ADJUSTL(buf)
  str6 = buf(1:6)
  call upper(str6)
  if(str6 == "[CELL]") exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_lat_vec_from_molden: [CELL] sectio&
                   &n not found'
  write(6,'(A)') 'in file '//TRIM(molden)
  close(fid)
  stop
 end if

 do i = 1, 3
  read(fid,*) lat_vec(:,i)
 end do ! for i
 close(fid)
end subroutine read_lat_vec_from_molden

! read lattice vectors (3,3) from a specified .xyz file
subroutine read_lat_vec_from_xyz(xyzname, lat_vec)
 implicit none
 integer :: i, j, fid
 real(kind=8), intent(out) :: lat_vec(3,3)
 character(len=240) :: buf
 character(len=240), intent(in) :: xyzname

 lat_vec = 0d0
 open(newunit=fid,file=TRIM(xyzname),status='old',position='rewind')
 read(fid,'(A)') buf
 read(fid,'(A)') buf
 close(fid)

 i = INDEX(buf, """")
 j = INDEX(buf(i+1:), """")
 if(i==0 .or. j<=i) then
  write(6,'(/,A)') 'ERROR in subroutine read_lat_vec_from_xyz: wrong lattice ve&
                   &ctors found in'
  write(6,'(A)') 'file '//TRIM(xyzname)
  stop
 end if
 read(buf(i+1:i+j-1),*) lat_vec
end subroutine read_lat_vec_from_xyz

! Read lattice vectors (3,3) from a specified CP2K input file. It is allowed
! that inpname does not end with '.inp'. Only the &CELL section is required in
! this file.
subroutine read_lat_vec_from_cp2k_inp(inpname, lat_vec)
 implicit none
 integer :: i, j, fid
 real(kind=8), intent(out) :: lat_vec(3,3)
 character(len=1) :: str1
 character(len=3) :: str3
 character(len=5) :: str5
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname

 lat_vec = 0d0; str3 = ' '; str5 = ' '
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  buf = ADJUSTL(buf)
  str5 = buf(1:5)
  call upper(str5)
  if(str5 == "&CELL") exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_lat_vec_from_cp2k_inp: `&CELL` sec&
                   &tion not found"
  write(6,'(A)') 'in file '//TRIM(inpname)
  close(fid)
  stop
 end if

 read(fid,'(A)') buf
 read(buf,*) str3
 str3 = ADJUSTL(str3)

 if(str3 == 'ABC') then
  i = INDEX(buf, ']')
  if(i == 0) then
   j = INDEX(buf, 'ABC')
   read(buf(j+3:),*) lat_vec(1,1), lat_vec(2,2), lat_vec(3,3)
  else
   read(buf(i+1:),*) lat_vec(1,1), lat_vec(2,2), lat_vec(3,3)
  end if
 else if(str3(1:1) == 'A') then
  BACKSPACE(fid)
  do i = 1, 3
   read(fid,*) str1, lat_vec(:,i)
  end do ! for i
 end if

 close(fid)
end subroutine read_lat_vec_from_cp2k_inp

! read lattice vectors from xTB .coord file
subroutine read_lat_vec_from_coord(coord, lat_vec)
 use phys_cons, only: Bohr_const
 implicit none
 integer :: fid
 real(kind=8), intent(out) :: lat_vec(3,3)
 character(len=240) :: buf
 character(len=240), intent(in) :: coord

 lat_vec = 0d0
 open(newunit=fid,file=TRIM(coord),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(buf(1:8) == '$lattice') exit
 end do ! for while

 read(fid,*) lat_vec
 close(fid)
 lat_vec = lat_vec*Bohr_const
end subroutine read_lat_vec_from_coord

! read lattice parameters (a,b,c,alpha,beta,gamma) from a specified .pdb file
subroutine read_lat_para_from_pdb(pdbname, cell)
 implicit none
 integer :: fid
 real(kind=8), intent(out) :: cell(6)
 character(len=240) :: buf
 character(len=240), intent(in) :: pdbname

 cell = 0d0
 open(newunit=fid,file=TRIM(pdbname),status='old',position='rewind')
 read(fid,'(A)') buf
 read(fid,'(A)') buf
 close(fid)

 if(buf(1:6) /= 'CRYST1') then
  write(6,'(/,A)') 'ERROR in subroutine read_lat_para_from_pdb: "CRYST1" is not&
                   & found in the second'
  write(6,'(A)') 'line in file '//TRIM(pdbname)
  stop
 end if

 read(buf(7:),*) cell(1:6)
end subroutine read_lat_para_from_pdb

! read lattice vectors from a specified .pdb file
subroutine read_lat_vec_from_pdb(pdbname, lat_vec)
 implicit none
 real(kind=8) :: cell(6)
 real(kind=8), intent(out) :: lat_vec(3,3)
 character(len=240), intent(in) :: pdbname

 call read_lat_para_from_pdb(pdbname, cell)
 call lat_para2lat_vec(cell, lat_vec)
end subroutine read_lat_vec_from_pdb

! Read lattice vectors (3,3) from a specified file, where the filetype would be
! detected automatically.
subroutine read_lat_vec_from_file(fname, lat_vec)
 implicit none
 integer :: i, k
 real(kind=8), intent(out) :: lat_vec(3,3)
!f2py intent(out) :: lat_vec
 character(len=240), intent(in) :: fname
!f2py intent(in) :: fname

 k = LEN_TRIM(fname)
 call find_specified_suffix(fname(1:k), '.', i)

 select case(fname(i:k))
 case('.coord')
  call read_lat_vec_from_coord(fname, lat_vec)
 case('.inp')
  call read_lat_vec_from_cp2k_inp(fname, lat_vec)
 case('.molden')
  call read_lat_vec_from_molden(fname, lat_vec)
 case('.pdb')
  call read_lat_vec_from_pdb(fname, lat_vec)
 case('.xyz')
  call read_lat_vec_from_xyz(fname, lat_vec)
 end select
end subroutine read_lat_vec_from_file

