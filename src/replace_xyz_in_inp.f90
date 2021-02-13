! written by jxzou at 20210213: replace Cartesian coordinates in a input file

program main
 implicit none
 integer :: i, itype
 integer, parameter :: iout = 6
 character(len=11) :: str
 character(len=240) :: xyzname, inpname

 i = iargc()
 if(i /= 3) then
  write(iout,'(/,A)')' ERROR in subroutine replace_xyz_in_inp: wrong command&
                     & line arguments!'
  write(iout,'(A)')  ' Example 1 (OpenMolcas): replace_xyz_in_inp a.xyz a.input -molcas'
  write(iout,'(A,/)')' Example 2 (OpenMolcas): replace_xyz_in_inp a.out a.input -molcas'
  stop
 end if

 call getarg(1, xyzname)
 call getarg(2, inpname)
 call getarg(3, str)
 call require_file_exist(xyzname)
 call require_file_exist(inpname)

 select case(TRIM(str))
 case('-molcas','-openmolcas')
  itype = 1
 case('-gamess', '-gms')
  itype = 2
 case default
  itype = 0
  write(iout,'(A)') 'ERROR in subroutine replace_xyz_in_inp: wrong command&
                   & line argument.'
  write(iout,'(A)') "The 3rd argument can only be -molcas, -openmolcas, -gms,&
                   & -gamess. But got argument="//TRIM(str)
  stop
 end select

 call replace_xyz_in_inp(xyzname, inpname, itype)
 stop
end program main

! replace Cartesian coordinates in a input file
subroutine replace_xyz_in_inp(xyzname, inpname, itype)
 implicit none
 integer :: i, natom
 integer, intent(in) :: itype ! 1/2 for Molcas/GAMESS
 integer, parameter :: iout = 6
 real(kind=8), allocatable :: coor(:,:)
 character(len=240), intent(in) :: xyzname, inpname

 i = LEN_TRIM(xyzname)

 if(xyzname(i-3:i) == '.xyz') then
  call read_natom_from_xyz(xyzname, natom)
  allocate(coor(3,natom))
  call read_coor_from_xyz(xyzname, natom, coor)
 else
  select case(itype)
  case(1)
   call read_natom_from_molcas_out(xyzname, natom)
   allocate(coor(3,natom))
   call read_coor_from_molcas_out(xyzname, natom, coor)
  case default
   write(iout,'(A)') 'ERROR in subroutine replace_xyz_in_inp: file format not&
                     & supported.'
   write(iout,'(A)') 'Filename='//TRIM(xyzname)
   stop
  end select
 end if

 select case(itype)
 case(1)
  call replace_coor_in_molcas_inp(inpname, natom, coor)
 case default
  write(iout,'(A)') 'ERROR in subroutine replace_xyz_in_inp: file format not&
                   & supported.'
  write(iout,'(A)') 'Filename='//TRIM(inpname)
  stop
 end select

 deallocate(coor)
 return
end subroutine replace_xyz_in_inp

! replace Cartesian coordinates in a (Open)Molcas .input file
subroutine replace_coor_in_molcas_inp(inpname, natom, coor)
 implicit none
 integer :: i, k, fid, fid1
 integer, intent(in) :: natom
 integer, parameter :: iout = 6
 real(kind=8), intent(in) :: coor(3,natom)
 character(len=6) :: elem
 character(len=16) :: key
 character(len=240) :: buf0, buf, inpname1
 character(len=240), intent(in) :: inpname

 i = index(inpname, '.inp', back=.true.)
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
  write(iout,'(A)') 'ERROR in subroutine replace_coor_in_molcas_inp: insuffic&
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
 return
end subroutine replace_coor_in_molcas_inp

