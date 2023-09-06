! written by jxzou at 20230808: a wrapper of mkl2fch and fch2bdf for ORCA->BDF

program mkl2bdf
 use util_wrapper, only: mkl2fch_wrap
 implicit none
 integer :: i, k, SYSTEM, RENAME
 character(len=240) :: mklname, fchname, inpname0, inpname, orbfile0, orbfile, &
  basfile0, basfile

 i = iargc()
 if(.not. (i==1 .or. i==2)) then
  write(6,'(/,A)') 'ERROR in program mkl2bdf: wrong command line argument!'
  write(6,'(A)')   'Example 1: mkl2bdf a.mkl'
  write(6,'(A,/)') 'Example 2: mkl2bdf a.mkl b.inp'
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
  inpname = mklname(1:i-1)//'.inp'
 end if

 fchname = ' '
 call find_specified_suffix(mklname, '.mkl', i)
 call get_a_random_int(k)
 write(fchname,'(A,I0,A)') mklname(1:i-1)//'_',k,'.fch'

 call mkl2fch_wrap(mklname, fchname)

 i = SYSTEM('fch2bdf '//TRIM(fchname))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in program mkl2bdf: failed to call utility fch2bdf.'
  write(6,'(A)') 'mkl2bdf is a wrapper of mkl2fch and fch2bdf, so fch2bdf is re&
                 &quired to be'
  write(6,'(A)') 'called successfully.'
  stop
 end if

 call delete_file(fchname)
 call find_specified_suffix(fchname, '.fch', i)
 inpname0 = fchname(1:i-1)//'_bdf.inp'
 orbfile0 = fchname(1:i-1)//'_bdf.scforb'
 basfile0 = fchname(1:i-1)//'.BAS'
 call upper(basfile0)

 call find_specified_suffix(inpname, '.inp', i)
 orbfile = inpname(1:i-1)//'.scforb'
 basfile = inpname(1:i-1)//'.BAS'
 call upper(basfile)

 i = RENAME(TRIM(inpname0), TRIM(inpname))
 i = RENAME(TRIM(orbfile0), TRIM(orbfile))
 i = RENAME(TRIM(basfile0), TRIM(basfile))

 call modify_basfile_in_bdf_inp(inpname, basfile)
end program mkl2bdf

! modify the orbital filename in BDF .inp file
subroutine modify_basfile_in_bdf_inp(inpname, basfile)
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname, basfile

 call find_specified_suffix(inpname, '.inp', i)
 inpname1 = inpname(1:i-1)//'.t'

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid,'(A)') buf
  write(fid1,'(A)') TRIM(buf)
  if(buf(1:5) == 'Basis') exit
 end do ! for while

 write(fid1,'(A)') TRIM(basfile)
 read(fid,'(A)') buf ! skip the original basfile name

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine modify_basfile_in_bdf_inp

