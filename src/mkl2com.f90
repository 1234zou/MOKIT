! written by jxzou at 20230727: a wrapper of mkl2fch and fch2com for ORCA ->
! Molpro

program mkl2com
 use util_wrapper, only: mkl2fch_wrap
 implicit none
 integer :: i, k, narg, irel, SYSTEM, RENAME
 character(len=240) :: mklname, fchname, comname0, comname
 character(len=240) :: a_file0, b_file0, a_file, b_file
 character(len=240), allocatable :: str_arg(:)
 logical :: alive

 narg = iargc()
 if(narg<1 .or. narg>3) then
  write(6,'(/,A)') ' ERROR in program mkl2com: wrong command line argument!'
  write(6,'(A)')   ' Example 1: mkl2com a.mkl'
  write(6,'(A)')   ' Example 2: mkl2com a.mkl b.com'
  write(6,'(A)')   ' Example 3: mkl2com a.mkl -dkh2'
  write(6,'(A)')   ' Example 4: mkl2com a.mkl -sfx2c'
  write(6,'(A,/)') ' Example 5: mkl2com a.mkl b.com -sfx2c'
  stop
 end if

 allocate(str_arg(narg))
 do i = 1, narg, 1
  call getarg(i, str_arg(i))
 end do ! for i

 mklname = str_arg(1)
 call require_file_exist(mklname)
 call find_specified_suffix(mklname, '.mkl', i)
 comname = mklname(1:i-1)//'.com'
 k = 2; irel = -1

 if(narg > 1) then
  i = LEN_TRIM(str_arg(2))
  if(str_arg(2)(i-3:i) == '.com') then
   comname = str_arg(2)
   k = 3
  end if
  if(k <= narg) then
   select case(TRIM(str_arg(k)))
   case('-sfx2c')
    irel = -3
   case('-dkh2')
    irel = 2
   case default
    write(6,'(/,A)') 'ERROR in subroutine mkl2com: wrong command line argument!'
    write(6,'(A)')   'Example: mkl2com a.mkl -sfx2c'
    stop
   end select
  end if
 end if

 deallocate(str_arg)
 fchname = ' '
 call find_specified_suffix(mklname, '.mkl', i)
 call get_a_random_int(k)
 write(fchname,'(A,I0,A)') mklname(1:i-1)//'_',k,'.fch'

 call mkl2fch_wrap(mklname=mklname,fchname=fchname,irel=irel)

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

