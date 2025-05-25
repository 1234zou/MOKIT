! Written by jxzou at 20171106.
! Modified by jxzou at 20180610.
! Updated by jxzou at 20220910: support modifying GAMESS-2022-R1

! This program is used to move the last few strings in a long line(> 72 charac-
!  ters) to the next line, such that GAMESS actvte.x will not recognize this
!  line as a long(>72) one.
! Note: since 2020-R2, actvte.x is removed. But this file is still needed

program modify_GMS2
 implicit none
 integer :: i, j, k, fid1, fid2, RENAME
 integer, parameter :: nfile = 23
 integer, parameter :: len_max = 72
 character(len=6) :: str6
 character(len=30) :: fname1, fname2
 character(len=200) :: buf
 character(len=10) :: filelist(nfile)
 logical :: alive
 data filelist /'comp','cphf','cprohf','fmoh2c','fmohss','grd1','grd2a','gvb',&
  'guess','hess','hss1c','hss2a','hss2b','hss2c','locpol','mexing','parley',&
  'prppop','qmfm','scflib','statpt','vector','vvos'/

 fname1 = ' '; fname2 = ' '; buf = ' '; str6 = ' '

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
   read(fid1,'(A)',iostat=k) buf
   if(k /= 0) exit
   if(buf(1:1)=='C' .or. buf(1:1)=='c' .or. buf(1:1)=='!') then
    write(fid2,'(A)') TRIM(buf)
    cycle
   end if

   j = LEN_TRIM(buf)
   k = SCAN(buf(1:j),'!')
   if(k /= 0) then
    j = min(j,k-1)
    ! in case there exists spaces, trim this string
    k = LEN_TRIM(buf(1:j))
    j = min(j,k)
   end if

   str6 = buf(7:12)
   call upper(str6)

   if(j>len_max .and. str6=='COMMON') then
    if(j==73 .and. buf(j:j)=="&") then
      write(fid2,'(A)') TRIM(buf)
    else if(buf(j:j) == ',') then
      buf(j:j) = ' '
      k = SCAN(buf(1:len_max),',',.true.)
      write(fid2,'(A)') buf(1:k)
      write(fid2,'(A)') '     *                '//TRIM(buf(k+1:))//','
    else
      k = SCAN(buf(1:len_max),',',.true.)
      write(fid2,'(A)') buf(1:k)
      write(fid2,'(A)') '     *                '//TRIM(buf(k+1:))
    end if
   else
    write(fid2,'(A)') TRIM(buf)
   end if
  end do ! for while

  close(fid1,status='delete')
  close(fid2)
  k = RENAME(fname2, fname1)
 end do ! for i
end program modify_GMS2

! transform a string into upper case
subroutine upper(buf)
 implicit none
 integer :: i, k
 character(len=*), intent(inout) :: buf

 do i = 1, LEN(buf), 1
  k = IACHAR(buf(i:i))
  if(k>=97 .and. k<=122) buf(i:i) = ACHAR(k-32)
 end do ! for i
end subroutine upper

