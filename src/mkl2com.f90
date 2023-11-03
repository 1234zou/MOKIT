! written by jxzou at 20230727: a wrapper of mkl2fch and fch2com for ORCA ->
! Molpro

program mkl2com
 use util_wrapper, only: mkl2fch_wrap
 implicit none
 integer :: i, k, SYSTEM, RENAME
 character(len=240) :: mklname, fchname, comname0, comname
 character(len=240) :: a_file0, b_file0, a_file, b_file
 logical :: alive

 i = iargc()
 if(.not. (i==1 .or. i==2)) then
  write(6,'(/,A)') 'ERROR in program mkl2com: wrong command line argument!'
  write(6,'(A)')   'Example 1: mkl2com a.mkl'
  write(6,'(A,/)') 'Example 2: mkl2com a.mkl b.com'
  stop
 end if

 mklname = ' '
 call getarg(1, mklname)
 call require_file_exist(mklname)

 comname = ' '
 if(i == 2) then
  call getarg(2, comname)
 else
  call find_specified_suffix(mklname, '.mkl', i)
  comname = mklname(1:i-1)//'.com'
 end if

 fchname = ' '
 call find_specified_suffix(mklname, '.mkl', i)
 call get_a_random_int(k)
 write(fchname,'(A,I0,A)') mklname(1:i-1)//'_',k,'.fch'

 call mkl2fch_wrap(mklname, fchname)

 i = SYSTEM('fch2com '//TRIM(fchname))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in program mkl2com: failed to call utility fch2com.'
  write(6,'(A)') 'mkl2com is a wrapper of mkl2fch and fch2com, so fch2com is re&
                 &quired to be'
  write(6,'(A)') 'called successfully.'
  stop
 end if

 call delete_file(fchname)
 call find_specified_suffix(fchname, '.fch', i)
 comname0 = fchname(1:i-1)//'.com'
 a_file0 = fchname(1:i-1)//'.a'
 b_file0 = fchname(1:i-1)//'.b'
 call lower(a_file0)
 call lower(b_file0)

 call find_specified_suffix(comname, '.com', i)
 a_file = comname(1:i-1)//'.a'
 b_file = comname(1:i-1)//'.b'
 call lower(a_file)
 call lower(b_file)

 i = RENAME(TRIM(comname0), TRIM(comname))
 i = RENAME(TRIM(a_file0), TRIM(a_file))
 inquire(file=TRIM(b_file0), exist=alive)
 if(alive) i = RENAME(TRIM(b_file0), TRIM(b_file))

 call modify_orbname_in_com(comname, a_file)
end program mkl2com

! modify the orbital filename in .com file
subroutine modify_orbname_in_com(comname, a_file)
 implicit none
 integer :: i, j, k, fid, fid1, RENAME
 character(len=240) :: comname1, b_file
 character(len=240), intent(in) :: comname, a_file
 character(len=500) :: buf
 ! Note: the length of a line in Molpro .com file might exceed 240 characters

 call find_specified_suffix(comname, '.com', i)
 comname1 = comname(1:i-1)//'.t'

 call find_specified_suffix(a_file, '.a', i)
 b_file = a_file(1:i-1)//'.b'

 open(newunit=fid,file=TRIM(comname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(comname1),status='replace')

 do while(.true.)
  read(fid,'(A)') buf
  write(fid1,'(A)') TRIM(buf)
  if(buf(1:2) == 'wf') exit
 end do ! for while

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit

  i = INDEX(buf, 'file=')
  if(i > 0) then
   j = INDEX(buf, '.a', back=.true.)
   if(j > 0) buf = buf(1:i+4)//TRIM(a_file)//';'
   k = INDEX(buf, '.b', back=.true.)
   if(k > 0) buf = buf(1:i+4)//TRIM(b_file)//';'
  end if
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(comname1), TRIM(comname))
end subroutine modify_orbname_in_com

