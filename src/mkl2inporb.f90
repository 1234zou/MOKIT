! written by jxzou at 20230728: a wrapper of mkl2fch and fch2inporb for ORCA ->
! (Open)Molcas

program mkl2inporb
 use util_wrapper, only: mkl2fch_wrap
 implicit none
 integer :: i, k, SYSTEM, RENAME
 character(len=240) :: mklname, fchname, inpname0, inpname, orbfile0, orbfile

 i = iargc()
 if(.not. (i==1 .or. i==2)) then
  write(6,'(/,A)') 'ERROR in program mkl2inporb: wrong command line argument!'
  write(6,'(A)')   'Example 1: mkl2inporb a.mkl'
  write(6,'(A,/)') 'Example 2: mkl2inporb a.mkl b.input'
  stop
 end if

 mklname = ' '
 call getarg(1, mklname)
 call require_file_exist(mklname)

 inpname = ' '
 if(i == 2) then
  call getarg(2, inpname)
 else
  call find_specified_suffix(mklname, '.mkl', i)
  inpname = mklname(1:i-1)//'.input'
 end if

 fchname = ' '
 call find_specified_suffix(mklname, '.mkl', i)
 call get_a_random_int(k)
 write(fchname,'(A,I0,A)') mklname(1:i-1)//'_',k,'.fch'

 call mkl2fch_wrap(mklname, fchname)

 i = SYSTEM('fch2inporb '//TRIM(fchname))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in program mkl2inporb: failed to call utility fch2inp&
                   &orb.'
  write(6,'(A)') 'mkl2inporb is a wrapper of mkl2fch and fch2inporb, so fch2inp&
                 &orb is required to'
  write(6,'(A)') 'be called successfully.'
  stop
 end if

 call delete_file(fchname)
 call find_specified_suffix(fchname, '.fch', i)
 inpname0 = fchname(1:i-1)//'.input'
 orbfile0 = fchname(1:i-1)//'.INPORB'

 call find_specified_suffix(inpname, '.inp', i)
 orbfile = inpname(1:i-1)//'.INPORB'

 i = RENAME(TRIM(inpname0), TRIM(inpname))
 i = RENAME(TRIM(orbfile0), TRIM(orbfile))
 call modify_orbname_in_molcas_inp(inpname, orbfile)
end program mkl2inporb

! modify the orbital filename in a (Open)Molcas .input file
subroutine modify_orbname_in_molcas_inp(inpname, orbfile)
 implicit none
 integer :: i, j, fid, fid1, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname, orbfile

 call find_specified_suffix(inpname, '.inp', i)
 inpname1 = inpname(1:i-1)//'.t'

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:7) == 'FILEORB') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(i /= 0) then
  close(fid)
  close(fid1,status='delete')
  write(6,'(/,A)') "ERROR in subroutine modify_orbname_in_molcas_inp: keyword '&
                   &FILEORB' not found"
  write(6,'(A)') 'in file '//TRIM(inpname)
  stop
 end if

 j = index(buf, '=')
 buf = buf(1:j)//' '//TRIM(orbfile)
 write(fid1,'(A)') TRIM(buf)

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine modify_orbname_in_molcas_inp

