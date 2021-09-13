! Written by jxzou at 20171106.
! modified by jxzou at 20180610.

! This program is used to move the last few strings in a long line(> 72 characters) to the
! next line, such that GAMESS actvte.x will not recognize this line as a long(>72) one.

program main
 implicit none
 integer :: i, j, k, fid1, fid2
 integer, parameter :: nfile = 23
 integer, parameter :: len_max = 72
 character(len=30) :: fname1, fname2
 character(len=200) :: buffer
 character(len=10) :: filelist(nfile)
 logical :: alive
 data filelist /'comp','cphf','cprohf','fmoh2c','fmohss','grd1','grd2a','gvb',&
  'guess','hess','hss1c','hss2a','hss2b','hss2c','locpol','mexing','parley',&
  'prppop','qmfm','scflib','statpt','vector','vvos'/

 fname1 = ' '
 fname2 = ' '
 buffer = ' '

 do i = 1, nfile, 1
  fname1 = TRIM(filelist(i))//'.src'
  fname2 = TRIM(filelist(i))//'_new.src'
  if(TRIM(fname1) == 'vector.src') then
   inquire(file='vector.src', exist=alive)
   if(.not. alive) cycle
  end if
  open(newunit=fid1,file=TRIM(fname1),status='old',position='rewind')
  open(newunit=fid2,file=TRIM(fname2),status='replace')

  do while(.true.)
   read(fid1,'(A)',iostat=k) buffer
   if(k /= 0) exit
   if(buffer(1:1)=='C' .or. buffer(1:1)=='c' .or. buffer(1:1)=='!') then
     write(fid2,'(A)') TRIM(buffer)
     cycle
   end if
   j = LEN_TRIM(buffer)
   k = SCAN(buffer(1:j),'!')
   if(k /= 0) j = min(j,k-1)
   if(j > len_max) then
     if(buffer(j:j) == ',') then
       buffer(j:j) = ' '
       k = SCAN(buffer(1:len_max),',',.true.)
       write(fid2,'(A)') buffer(1:k)
       write(fid2,'(A)') '     *                '//TRIM(buffer(k+1:))//','
     else
       k = SCAN(buffer(1:len_max),',',.true.)
       write(fid2,'(A)') buffer(1:k)
       write(fid2,'(A)') '     *                '//TRIM(buffer(k+1:))
     end if
   else
    write(fid2,'(A)') TRIM(buffer)
   end if
  end do ! for while

  close(fid1,status='delete')
  close(fid2)
  k = RENAME(fname2,fname1)
 end do ! for i

 stop
end program main

