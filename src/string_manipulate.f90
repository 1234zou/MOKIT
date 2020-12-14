! written by jxzou at 20201208: move string manipulation subroutines into this file

! transform a string into upper case
subroutine upper(buf)
 implicit none
 integer :: i, k
 character(len=*), intent(inout) :: buf

 do i = 1, LEN(buf), 1
  k = ICHAR(buf(i:i))
  if(k>=97 .and. k<=122) buf(i:i) = CHAR(k-32)
 end do
 return
end subroutine upper

! transform a string into lower case
subroutine lower(buf)
 implicit none
 integer :: i, k
 character(len=*), intent(inout) :: buf

 do i = 1, LEN(buf), 1
  k = ICHAR(buf(i:i))
  if (k>=65 .and. k<=90) buf(i:i) = CHAR(k+32)
 end do
 return
end subroutine lower

! convert a (character) stype to (integer) itype
subroutine stype2itype(stype, itype)
 implicit none
 integer, parameter :: iout = 6
 integer, intent(out) :: itype
 character(len=1), intent(in) :: stype

 ! 'S', 'P', 'D', 'F', 'G', 'H', 'I'
 !  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7
 select case(stype)
 case('S')
  itype = 1
 case('P')
  itype = 2
 case('D')
  itype = 3
 case('F')
  itype = 4
 case('G')
  itype = 5
 case('H')
  itype = 6
 case('I')
  itype = 7
 case('L') ! 'L' is 'SP'
  itype = 0
 case default
  write(iout,'(A)') 'ERROR in subroutine stype2itype: stype out of range.'
  write(iout,'(A)') 'stype= '//TRIM(stype)
  stop
 end select

 return
end subroutine stype2itype

! check whether there exists DKH keywords in a given .inp file
subroutine check_DKH_in_inp(inpname, order)
 implicit none
 integer :: i, k, fid
 integer, intent(out) :: order
! -2: no DKH
! -1: RESC
!  0: DKH 0th-order
!  2: DKH2
!  4: DKH4 with SO
 integer, parameter :: iout = 6
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname
 character(len=1200) :: longbuf
 logical :: alive(6)

 longbuf = ' '
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  longbuf = TRIM(longbuf)//TRIM(buf)
  call upper(buf)
  if(index(buf,'$END') /= 0) exit
 end do ! for while
 close(fid)

 call upper(longbuf)

 order = -2
 if(index(longbuf,'RELWFN') == 0) then
  return
 else
  if(index(longbuf,'RELWFN=DK') == 0) then
   write(iout,'(A)') 'Warning in subroutine check_DKH_in_inp: unsupported&
                    & relativistic method detected.'
   write(iout,'(A)') 'Molcas/OpenMolcas does not support RELWFN=LUT-IOTC, IOTC,&
                    & RESC, or NESC. Only RELWFN=DK is supported.'
   write(iout,'(A)') 'The MO transferring will still be proceeded. But the result&
                    & may be non-sense.'
  end if
 end if

 order = 2
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  call upper(buf)
  if(index(buf,'$RELWFN') /= 0) exit
 end do ! for while
 close(fid)

 if(i == 0) then
  k = index(buf,'NORDER=')
  if(k /= 0) then
   read(buf(k+7:),*) order
  end if
 end if

 return
end subroutine check_DKH_in_inp

! convert a filename into which molpro requires, i.e.
! length <32 and in lowercase
subroutine convert2molpro_fname(fname, suffix)
 implicit none
 integer :: i, len1, len2
 character(len=240), intent(inout) :: fname
 character(len=*), intent(in) :: suffix

 if(fname(1:1) == ' ') fname = ADJUSTL(fname)

 len1 = INDEX(fname, '.', back=.true.) - 1
 if(len1 == -1) len1 = LEN_TRIM(fname)
 len2 = LEN(suffix)

 if(len1+len2 > 31) then
  fname = fname(1:32-len2)//suffix
 else
  fname = fname(1:len1)//suffix
 end if

 call lower(fname(1:32))
 return
end subroutine convert2molpro_fname

