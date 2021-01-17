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

! check whether there exists DKH keywords in a given GAMESS .inp file
subroutine check_DKH_in_gms_inp(inpname, order)
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
   write(iout,'(A)') 'Warning in subroutine check_DKH_in_gms_inp: unsupported&
                    & relativistic method detected.'
   write(iout,'(A)') '(Open)Molcas does not support RELWFN=LUT-IOTC, IOTC,&
                    & RESC, or NESC in GAMESS. Only RELWFN=DK is supported.'
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
end subroutine check_DKH_in_gms_inp

! check whether X2C appears in a given GAMESS .inp file
! Note: GAMESS does not support X2C, this is just for the utility bas_gms2molcas
!  to recognize X2C and pass it into (Open)Molcas .input file
subroutine check_X2C_in_gms_inp(inpname, X2C)
 implicit none
 integer :: i, fid
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname
 logical, intent(out) :: X2C

 X2C = .false. ! default

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 do i = 1, 5
  read(fid,'(A)') buf
  if(index(buf,'X2C') /= 0) X2C = .true.
  if(index(buf,'$END') /= 0) exit
 end do ! for i

 close(fid)
 return
end subroutine check_X2C_in_gms_inp

! convert a filename into which molpro requires, i.e.
! length <32 and in lowercase
subroutine convert2molpro_fname(fname, suffix)
 implicit none
 integer :: i, len1, len2
 integer, parameter :: iout = 6
 character(len=240), intent(inout) :: fname
 character(len=*), intent(in) :: suffix

 if(LEN_TRIM(fname) == 0) then
  write(iout,'(A)') 'ERROR in subroutine convert2molpro_fname: input fname is NULL.'
  stop
 end if

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

! add DKH2 related keywords into a given GAMESS .inp file,
! and switch the default SOSCF into DIIS
subroutine add_DKH2_into_gms_inp(inpname)
 implicit none
 integer :: i, k, fid1, fid2, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 inpname1 = TRIM(inpname)//'.tmp'
 open(newunit=fid1,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname1),status='replace')

 do i = 1, 3
  read(fid1,'(A)') buf

  if(index(buf, 'RELWFN=DK') /= 0) then
   close(fid1)
   close(fid2,status='delete')
   return
  end if

  k = index(buf,'$END')
  if(k /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for i

 write(fid2,'(A)') buf(1:k-1)//' RELWFN=DK $END'

 do while(.true.)
  read(fid1,'(A)') buf
  k = index(buf,'$END')
  if(k /= 0) exit
  if(index(buf,'$DATA') /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 if(k == 0) then
  write(fid2,'(A)') '$SCF DIRSCF=.TRUE. DIIS=.T. SOSCF=.F. $END'
 else
  write(fid2,'(A)') buf(1:k-1)//' DIIS=.T. SOSCF=.F. $END'
 end if

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
 return
end subroutine add_DKH2_into_gms_inp

! add DKH2 keyword into a given Gaussian .fch(k) file
subroutine add_DKH2_into_fch(fchname)
 implicit none
 integer :: i, j, k, nline, nterm, fid, fid1, RENAME
 integer, parameter :: iout = 6
 character(len=240) :: buf, fchname1
 character(len=240), intent(in) :: fchname
 character(len=1200) :: longbuf, longbuf1
 logical :: no_route, alive(3)

 buf = ' '
 fchname1 = ' '
 longbuf = ' '
 nterm = 0
 i = index(fchname, '.fch', back=.true.)
 fchname1 = fchname(1:i-1)//'.tmp'

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(fchname1),status='replace')

 no_route = .false.
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:6) == 'Charge') then
   no_route = .true.
   exit
  end if

  if(buf(1:5) == 'Route') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(no_route) then
  write(fid1,'(A5,38X,A)') 'Route','C   N=           3'
  write(fid1,'(A)') '#p int(nobasistransform,DKH2) nosymm'
 else
  if(i /= 0) then
   write(iout,'(A)') 'ERROR in subroutine add_DKH2_into_fch: incomplete .fch(k)  file.'
   write(iout,'(A)') "Neither 'Route' nor 'Charge' is detected in file "//TRIM(fchname)
   close(fid)
   close(fid1,status='delete')
   stop
  else
   k = index(buf, '='); read(buf(k+1:),*) nterm
   do while(.true.)
    read(fid,'(A)',iostat=i) buf
    if(i /= 0) exit
    if(buf(1:6) == 'Charge') exit
    longbuf = TRIM(longbuf)//TRIM(buf)
   end do ! for while

   if(i /= 0) then
    write(iout,'(A)') 'ERROR in subroutine add_DKH2_into_fch: incomplete .fch(k) file.'
    write(iout,'(A)') "No 'Charge' is detected in file "//TRIM(fchname)
    close(fid)
    close(fid1,status='delete')
    stop
   else
    longbuf1 = longbuf
    call upper(longbuf1)
    alive = [(index(longbuf1,'DKH2')/=0),(index(longbuf1,'DOUGLASKROLLHESS')/=0),&
             (index(longbuf1,'DKH')/=0 .and. index(longbuf1,'DKH4')==0 .and. &
              index(longbuf1,'NODKH')==0 .and. index(longbuf1,'DKHSO')==0)]
    if(ALL(alive .eqv. .false.)) then
     nterm = nterm + 1
     longbuf = TRIM(longbuf)//' int=DKH2'
    end if
    write(fid1,'(A5,38X,A,I2)') 'Route','C   N=          ', nterm
    k = LEN_TRIM(longbuf)
    nline = k/60
    if(k-60*nline > 0) nline = nline + 1
    do i = 1, nline, 1
     j = min(60*i,k)
     write(fid1,'(A)') longbuf(60*i-59:j)
    end do ! for i
   end if
  end if
 end if

 ! copy remaining content
 BACKSPACE(fid)
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(fchname1), TRIM(fchname))
 return
end subroutine add_DKH2_into_fch

! add X2C keyword into a given Gaussian .fch(k) file
! Note: Obviously, Gaussian cannot use X2C. This is just for other utilities to
!  recognize the X2C keyword, so that other utilities can add related keywords
!  when generating input files of other programs
subroutine add_X2C_into_fch(fchname)
 implicit none
 integer :: i, j, k, nline, nterm, fid, fid1, RENAME
 integer, parameter :: iout = 6
 character(len=240) :: buf, fchname1
 character(len=240), intent(in) :: fchname
 character(len=1200) :: longbuf, longbuf1
 logical :: no_route

 buf = ' '
 fchname1 = ' '
 longbuf = ' '
 nterm = 0
 i = index(fchname, '.fch', back=.true.)
 fchname1 = fchname(1:i-1)//'.tmp'

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(fchname1),status='replace')

 no_route = .false.
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:6) == 'Charge') then
   no_route = .true.
   exit
  end if

  if(buf(1:5) == 'Route') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(no_route) then
  write(fid1,'(A5,38X,A)') 'Route','C   N=           3'
  write(fid1,'(A)') '#p int(nobasistransform,X2C) nosymm'
 else
  if(i /= 0) then
   write(iout,'(A)') 'ERROR in subroutine add_X2C_into_fch: incomplete .fch(k)  file.'
   write(iout,'(A)') "Neither 'Route' nor 'Charge' is detected in file "//TRIM(fchname)
   close(fid)
   close(fid1,status='delete')
   stop
  else
   k = index(buf, '='); read(buf(k+1:),*) nterm
   do while(.true.)
    read(fid,'(A)',iostat=i) buf
    if(i /= 0) exit
    if(buf(1:6) == 'Charge') exit
    longbuf = TRIM(longbuf)//TRIM(buf)
   end do ! for while

   if(i /= 0) then
    write(iout,'(A)') 'ERROR in subroutine add_X2C_into_fch: incomplete .fch(k) file.'
    write(iout,'(A)') "No 'Charge' is detected in file "//TRIM(fchname)
    close(fid)
    close(fid1,status='delete')
    stop
   else
    longbuf1 = longbuf
    call upper(longbuf1)
    j = index(longbuf1, 'DKH'); k = index(longbuf1, 'NODKH')
    if(j/=0 .and. k==0) then
     longbuf(j:j+2) = 'X2C'
    else
     nterm = nterm + 1
     longbuf = TRIM(longbuf)//' int=X2C'
    end if
    write(fid1,'(A5,38X,A,I2)') 'Route','C   N=          ', nterm
    k = LEN_TRIM(longbuf)
    nline = k/60
    if(k-60*nline > 0) nline = nline + 1
    do i = 1, nline, 1
     j = min(60*i,k)
     write(fid1,'(A)') longbuf(60*i-59:j)
    end do ! for i
   end if
  end if
 end if

 ! copy remaining content
 BACKSPACE(fid)
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(fchname1), TRIM(fchname))
 return
end subroutine add_X2C_into_fch

! add '.x2c()' into a given PySCF .py file
subroutine add_X2C_into_py(pyname)
 implicit none
 integer :: i, fid, fid1, RENAME
 integer, parameter :: iout = 6
 character(len=240) :: buf, pyname1
 character(len=240), intent(in) :: pyname

 pyname1 = TRIM(pyname)//'.tmp'
 open(newunit=fid,file=TRIM(pyname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(pyname1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:4) == 'mf =') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine add_X2C_into_py: 'mf =' not found in&
                   & file "//TRIM(pyname)
  stop
 end if

 write(fid1,'(A)') TRIM(buf)//'.x2c()'

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(pyname1), TRIM(pyname))
 return
end subroutine add_X2C_into_py

! detect the number of columns of data in a string buf
function detect_ncol_in_buf(buf) result(ncol)
 implicit none
 integer :: i, ncol
 character(len=24), allocatable :: sbuf(:)
 character(len=240), intent(in) :: buf

 if(LEN_TRIM(buf) == 0) then
  ncol = 0
  return
 end if

 ncol = 1
 do while(.true.)
  allocate(sbuf(ncol))
  read(buf,*,iostat=i) sbuf(1:ncol)
  deallocate(sbuf)
  if(i /= 0) exit
  ncol = ncol + 1
 end do ! for while

 ncol = ncol - 1
 return
end function detect_ncol_in_buf

