! written by jxzou at 20210213: replace Cartesian coordinates in a input file
! written by jxzou at 20210217: support Molpro

program main
 implicit none
 integer :: i, itype
 character(len=11) :: str
 character(len=240) :: xyzname, inpname

 i = iargc()
 if(i /= 3) then
  write(6,'(/,A)')' ERROR in subroutine replace_xyz_in_inp: wrong command line &
                   &arguments!'
  write(6,'(A)')  ' Example 1 (OpenMolcas): replace_xyz_in_inp a.xyz a.input -molcas'
  write(6,'(A)')  ' Example 2 (OpenMolcas): replace_xyz_in_inp a.out a.input -molcas'
  write(6,'(A)')  ' Example 3 (Molpro)    : replace_xyz_in_inp a.xyz a.com -molpro'
  write(6,'(A,/)')' Example 4 (Molpro)    : replace_xyz_in_inp a.out a.com -molpro'
  stop
 end if

 str = ' '
 call getarg(1, xyzname)
 call getarg(2, inpname)
 call getarg(3, str)
 call require_file_exist(xyzname)
 call require_file_exist(inpname)

 select case(TRIM(str))
 case('-molcas','-openmolcas')
  itype = 1
 case('-molpro')
  itype = 2
 case('-orca')
  itype = 3
 case default
  itype = 0
  write(6,'(A)') 'ERROR in subroutine replace_xyz_in_inp: wrong command line ar&
                 &gument str='//TRIM(str)
  stop
 end select

 call replace_xyz_in_inp(xyzname, inpname, itype)
end program main

! replace Cartesian coordinates in a input file
subroutine replace_xyz_in_inp(xyzname, inpname, itype)
 implicit none
 integer :: i, natom
 integer, intent(in) :: itype ! 1/2/3 for Molcas/Molpro/ORCA
 real(kind=8), allocatable :: coor(:,:)
 character(len=2), allocatable :: elem(:)
 character(len=240), intent(in) :: xyzname, inpname

 i = LEN_TRIM(xyzname)

 if(xyzname(i-3:i) == '.xyz') then
  call read_natom_from_xyz(xyzname, natom)
  allocate(elem(natom), coor(3,natom))
  call read_elem_and_coor_from_xyz(xyzname, natom, elem, coor)
  deallocate(elem)
 else if(xyzname(i-3:i) == '.EIn') then
  call read_natom_from_EIn(xyzname, natom)
  allocate(coor(3,natom))
  call read_coor_from_EIn(xyzname, natom, coor)
 else
  select case(itype)
  case(1)
   call read_natom_from_molcas_out(xyzname, natom)
   allocate(coor(3,natom))
   call read_coor_from_molcas_out(xyzname, natom, coor)
  case(2)
   call read_natom_from_molpro_out(xyzname, natom)
   allocate(coor(3,natom))
   call read_coor_from_molpro_out(xyzname, natom, coor)
  case default
   write(6,'(/,A)') 'ERROR in subroutine replace_xyz_in_inp: file format not su&
                    &pported.'
   write(6,'(A)') 'Filename='//TRIM(xyzname)
   stop
  end select
 end if

 select case(itype)
 case(1)
  call replace_coor_in_molcas_inp(inpname, natom, coor)
 case(2)
  call replace_coor_in_molpro_inp(inpname, natom, coor)
 case(3)
  call replace_coor_in_orca_inp(inpname, natom, coor)
 case default
  write(6,'(/,A)') 'ERROR in subroutine replace_xyz_in_inp: itype out of range.'
  write(6,'(A)') 'Filename='//TRIM(inpname)
  stop
 end select

 deallocate(coor)
end subroutine replace_xyz_in_inp

! replace Cartesian coordinates in a (Open)Molcas input file
subroutine replace_coor_in_molcas_inp(inpname, natom, coor)
 implicit none
 integer :: i, k, fid, fid1
 integer, intent(in) :: natom
 real(kind=8), intent(in) :: coor(3,natom)
 character(len=6) :: elem
 character(len=16) :: key
 character(len=240) :: buf0, buf, inpname1
 character(len=240), intent(in) :: inpname

 i = INDEX(inpname, '.inp', back=.true.)
 inpname1 = inpname(1:i-1)//'_new.input'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 i = 0
 do while(.true.)
  read(fid,'(A)',iostat=k) buf
  if(k /= 0) exit
  write(fid1,'(A)') TRIM(buf)

  key = buf(1:16)
  call lower(key)
  if(key == 'end of basis set') then
   read(buf0,*) elem
   BACKSPACE(fid1)
   BACKSPACE(fid1)
   i = i + 1
   write(fid1,'(A,3(2X,F15.8),A)') TRIM(elem),coor(1:3,i),'   Angstrom'
   write(fid1,'(A)') 'End of basis set'
   if(i == natom) exit
  end if

  buf0 = buf
 end do ! for while

 if(k /= 0) then
  write(6,'(A)') 'ERROR in subroutine replace_coor_in_molcas_inp: insuffic&
                    &ient number of atoms found in file '//TRIM(inpname)
  close(fid1,status='delete')
  close(fid)
  stop
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid)
 close(fid1)
end subroutine replace_coor_in_molcas_inp

! replace Cartesian coordinates in a Molpro input file
subroutine replace_coor_in_molpro_inp(inpname, natom, coor)
 implicit none
 integer :: i, fid, fid1
 integer, intent(in) :: natom
 real(kind=8), intent(in) :: coor(3,natom)
 character(len=6) :: elem
 character(len=8) :: key
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 i = INDEX(inpname, '.com', back=.true.)
 inpname1 = inpname(1:i-1)//'_new.com'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)

  key = buf(1:8)
  call lower(key)
  if(key == 'geometry') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine replace_coor_in_molcas_inp: insuffic&
                    &ient number of atoms found'
  write(6,'(A)') 'in file '//TRIM(inpname)
  close(fid1,status='delete')
  close(fid)
  stop
 end if

 do i = 1, natom, 1
  elem = ' '
  read(fid,*) elem
  write(fid1,'(A,3(2X,F15.8))') TRIM(elem),coor(1:3,i)
 end do ! for i

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid)
 close(fid1)
end subroutine replace_coor_in_molpro_inp

! replace the Cartesian coordinates in an ORCA input file
subroutine replace_coor_in_orca_inp(inpname, natom, coor)
 implicit none
 integer :: i, k, fid, fid1, RENAME
 integer, intent(in) :: natom
 real(kind=8), intent(in) :: coor(3,natom)
 character(len=10) :: str
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 k = LEN_TRIM(inpname)
 i = INDEX(inpname(1:k), '/', back=.true.)
 if(i > 0) then
  inpname1 = inpname(i+1:k)//'.t'
 else
  inpname1 = TRIM(inpname)//'.t'
 end if

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid,'(A)') buf
  write(fid1,'(A)') TRIM(buf)
  if(buf(1:7) == ' Coords') exit
 end do ! for while

 do i = 1, natom, 1
  read(fid,*) str
  if(TRIM(str) == 'NewECP') then
   write(fid1,'(A)') '  NewECP'
   do while(.true.)
    read(fid,'(A)') buf
    write(fid1,'(A)') TRIM(buf)
    if(buf(3:5) == 'end') exit
   end do ! for while
   read(fid,*) str
  end if
  write(fid1,'(1X,A,3(1X,F17.9))') TRIM(str), coor(:,i)
  do while(.true.)
   read(fid,'(A)') buf
   write(fid1,'(A)') TRIM(buf)
   if(buf(3:5) == 'end') exit
  end do ! for while
 end do ! for i

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid, status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine replace_coor_in_orca_inp

! read Cartesian coordinates from a Gaussian .EIn file
subroutine read_coor_from_EIn(EIn, natom, coor)
 use phys_cons, only: Bohr_const
 implicit none
 integer :: i, nuc, job_type, charge, spin, fid
 integer, intent(in) :: natom
 real(kind=8), intent(out) :: coor(3,natom)
 character(len=240), intent(in) :: EIn

 open(newunit=fid,file=TRIM(EIn),status='old',position='rewind')
 read(fid,*) i, job_type, charge, spin

 if(i /= natom) then
  write(6,'(/,A)') 'ERROR in subroutine read_coor_from_EIn: input natom and nat&
                   &om in the .EIn'
  write(6,'(A)') 'file are inconsistent. EIn: '//TRIM(EIn)
  stop
 end if

 do i = 1, natom, 1
  read(fid,*) nuc, coor(:,i)
 end do ! for i

 close(fid)
 coor = coor*Bohr_const ! convert to Angstrom
end subroutine read_coor_from_EIn

