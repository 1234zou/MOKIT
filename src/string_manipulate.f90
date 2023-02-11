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
end subroutine lower

! convert a (character) stype to (integer) itype
subroutine stype2itype(stype, itype)
 implicit none
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
  write(6,'(A)') 'ERROR in subroutine stype2itype: stype out of range.'
  write(6,'(A)') 'stype= '//TRIM(stype)
  stop
 end select
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
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname
 character(len=1200) :: longbuf

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
   write(6,'(A)') 'Warning in subroutine check_DKH_in_gms_inp: unsupported&
                    & relativistic method detected.'
   write(6,'(A)') '(Open)Molcas does not support RELWFN=LUT-IOTC, IOTC,&
                    & RESC, or NESC in GAMESS. Only RELWFN=DK is supported.'
   write(6,'(A)') 'The MO transferring will still be proceeded. But the result&
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
end subroutine check_X2C_in_gms_inp

subroutine check_sph_in_gjf(gjfname, sph)
 implicit none
 integer :: i, fid
 character(len=240) :: buf
 character(len=1200) :: longbuf
 character(len=240), intent(in) :: gjfname
 logical, intent(out) :: sph

 sph = .true.
 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:1) == '#') exit
 end do ! for while

 longbuf = buf
 do i = 1, 5
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
  longbuf = TRIM(longbuf)//TRIM(buf)
 end do ! for i

 close(fid)
 if(index(longbuf,'6D')>0 .or. index(longbuf,'6d')>0) sph = .false.
end subroutine check_sph_in_gjf

! convert a filename into which molpro requires, i.e.
! length <32 and in lowercase
subroutine convert2molpro_fname(fname, suffix)
 implicit none
 integer :: len1, len2
 character(len=240), intent(inout) :: fname
 character(len=*), intent(in) :: suffix

 if(LEN_TRIM(fname) == 0) then
  write(6,'(A)') 'ERROR in subroutine convert2molpro_fname: input fname is NULL.'
  stop
 end if

 if(fname(1:1) == ' ') fname = ADJUSTL(fname)

 len1 = INDEX(fname, '.', back=.true.) - 1
 if(len1 == -1) len1 = LEN_TRIM(fname)
 len2 = LEN(suffix)

 if(len1+len2 > 15) then
  fname = fname(1:16-len2)//suffix
 else
  fname = fname(1:len1)//suffix
 end if

 call lower(fname(1:16))
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
end subroutine add_DKH2_into_gms_inp

! add DKH2 keyword into a given Gaussian .fch(k) file
subroutine add_DKH2_into_fch(fchname)
 implicit none
 integer :: i, j, k, nline, nterm, fid, fid1, RENAME
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
   write(6,'(A)') 'ERROR in subroutine add_DKH2_into_fch: incomplete .fch(k)  file.'
   write(6,'(A)') "Neither 'Route' nor 'Charge' is detected in file "//TRIM(fchname)
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
    write(6,'(A)') 'ERROR in subroutine add_DKH2_into_fch: incomplete .fch(k) file.'
    write(6,'(A)') "No 'Charge' is detected in file "//TRIM(fchname)
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
end subroutine add_DKH2_into_fch

! add X2C keyword into a given Gaussian .fch(k) file
! Note: Obviously, Gaussian cannot use X2C. This is just for other utilities to
!  recognize the X2C keyword, so that other utilities can add related keywords
!  when generating input files of other programs
subroutine add_X2C_into_fch(fchname)
 implicit none
 integer :: i, j, k, nline, nterm, fid, fid1, RENAME
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
   write(6,'(A)') 'ERROR in subroutine add_X2C_into_fch: incomplete .fch(k)  file.'
   write(6,'(A)') "Neither 'Route' nor 'Charge' is detected in file "//TRIM(fchname)
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
    write(6,'(A)') 'ERROR in subroutine add_X2C_into_fch: incomplete .fch(k) file.'
    write(6,'(A)') "No 'Charge' is detected in file "//TRIM(fchname)
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
end subroutine add_X2C_into_fch

! add '.x2c1e()' into a given PySCF .py file
subroutine add_X2C_into_py(pyname)
 implicit none
 integer :: i, fid, fid1, RENAME
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
  write(6,'(A)') "ERROR in subroutine add_X2C_into_py: 'mf =' not found in&
                   & file "//TRIM(pyname)
  stop
 end if

 write(fid1,'(A)') TRIM(buf)//'.x2c1e()'

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(pyname1), TRIM(pyname))
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
end function detect_ncol_in_buf

! modify the memory in a given .inp file
subroutine modify_memory_in_gms_inp(inpname, mem, nproc)
 implicit none
 integer :: i, fid1, fid2, RENAME
 integer, intent(in) :: mem, nproc
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 inpname1 = TRIM(inpname)//'.tmp'
 open(newunit=fid1,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(index(buf,'MWORDS') /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine modify_memory_in_gms_inp: no 'MWORDS' fou&
                 &nd in file "//TRIM(inpname)
  close(fid1)
  close(fid2,status='delete')
  stop
 end if

 write(fid2,'(A,I0,A)') ' $SYSTEM MWORDS=',FLOOR(DBLE(mem)*1d3/(8d0*DBLE(nproc))),' $END'

 ! copy the remaining content
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine modify_memory_in_gms_inp

! modify memory in a given PSI4 input file
! Note: input mem is in unit GB
subroutine modify_memory_in_psi4_inp(inpname, mem)
 implicit none
 integer :: i, fid, fid1, RENAME
 integer, intent(in) :: mem
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 inpname1 = TRIM(inpname)//'.t'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 read(fid,'(A)') buf
 write(fid1,'(A)') TRIM(buf)
 read(fid,'(A)') buf
 write(fid1,'(A,I0,A)') 'memory ', mem, ' GB'

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine modify_memory_in_psi4_inp

! modify memory in a given Q-Chem input file
! Note: input mem is in unit GB
subroutine modify_memory_in_qchem_inp(mem, inpname)
 implicit none
 integer :: i, fid, fid1, RENAME
 integer, intent(in) :: mem
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 i = index(inpname, '.in', back=.true.)
 inpname1 = inpname(1:i-1)//'.t'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:7) == 'mem_tot') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 write(fid1,'(A,I0)') 'mem_total ',mem*1000 ! in MB

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine modify_memory_in_qchem_inp

! add given/specified RIJK basis set into a PSI4 input file
subroutine add_RIJK_bas_into_psi4_inp(inpname, RIJK_bas)
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=240) :: inpname1
 character(len=21), intent(in) :: RIJK_bas
 character(len=240), intent(in) :: inpname

 if(LEN_TRIM(RIJK_bas) == 0) then
  write(6,'(A)') 'ERROR in subroutine add_RIJK_bas_into_psi4_inp:'
  stop
 end if

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine add_RIJK_bas_into_psi4_inp

! add given/specified RIJK basis set into an ORCA input file
subroutine add_RIJK_bas_into_orca_inp(inpname, RIJK_bas)
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=240) :: inpname1
 character(len=21), intent(in) :: RIJK_bas
 character(len=240), intent(in) :: inpname

 if(LEN_TRIM(RIJK_bas) == 0) then
  write(6,'(A)') 'ERROR in subroutine add_RIJK_bas_into_orca_inp: input RI&
                   & basis set is null string.'
  stop
 end if
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine add_RIJK_bas_into_orca_inp

! detect whether there exists the charge keyword in a given .gjf file
function detect_charge_key_in_gjf(gjfname) result(has_charge)
 implicit none
 integer :: i, j, nblank, fid
 character(len=240) :: buf
 character(len=240), intent(in) :: gjfname
 logical :: has_charge

 has_charge = .true.; nblank = 0
 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:1) == '#') then
   call lower(buf)
   if(index(buf,'charge') > 0) then
    close(fid)
    return
   end if
  end if
  if(LEN_TRIM(buf) == 0) nblank = nblank + 1
  if(nblank == 1) exit
 end do ! for while

 read(fid,'(A)') buf
 close(fid)
 i = index(buf,'{')
 j = index(buf,'}')
 if(i>0 .and. j>0) then
  call lower(buf(i+1:j-1))
  if(index(buf(i+1:j-1),'charge') > 0) return
 end if
 has_charge = .false.
end function detect_charge_key_in_gjf

! copy mixed/user-defined basis set in a given .gjf file to a .bas file
subroutine record_gen_basis_in_gjf(gjfname, basname, add_path)
 implicit none
 integer :: i, nblank0, nblank, fid1, fid2
 character(len=240) :: buf
 character(len=240), intent(in) :: gjfname
 character(len=240), intent(out) :: basname
 logical, intent(in) :: add_path
 logical, external :: detect_charge_key_in_gjf
 logical :: nobasis

 if(detect_charge_key_in_gjf(gjfname)) then
  nblank0 = 4
 else
  nblank0 = 3
 end if
 i = index(gjfname, '.gjf', back=.true.)
 basname = gjfname(1:i-1)//'.bas'

 open(newunit=fid1,file=TRIM(gjfname),status='old',position='rewind')
 nblank = 0
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) nblank = nblank + 1
  if(nblank == nblank0) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine record_gen_basis_in_gjf: incomplete&
                 & file: '//TRIM(gjfname)
  close(fid1)
  stop
 end if

 read(fid1,'(A)',iostat=i) buf
 nobasis = .false.
 if(i /= 0) nobasis = .true.
 if((.not.nobasis) .and. LEN_TRIM(buf)==0) nobasis = .true.

 if(nobasis) then
  write(6,'(A)') 'ERROR in subroutine record_gen_basis_in_gjf: no mixed/user&
                 &-defined basis'
  write(6,'(A)') 'set detected in file '//TRIM(gjfname)
  close(fid1)
  stop
 end if

 buf = ADJUSTL(buf)
 if(buf(1:1) == '-') then
  i = IACHAR(buf(2:2))
 else
  i = IACHAR(buf(1:1))
 end if

 if(.not. ((i>96 .and. i<123) .or. (i>64 .and. i<91))) then
  write(6,'(A)') 'ERROR in subroutine record_gen_basis_in_gjf: the first charac&
                 &ter in mixed/user-defined'
  write(6,'(A)') 'basis set is neither a-z, nor A-Z. This is not an element sym&
                 &bol. This format of basis'
  write(6,'(A)') 'set cannot be recognized by automr. Problematic file: '//&
                  TRIM(gjfname)
  close(fid1)
  stop
 end if

 open(newunit=fid2,file=TRIM(basname),status='replace')
 write(fid2,'(A)') TRIM(buf)

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 write(fid2,'(/)')
 close(fid1)
 close(fid2)

 call add_hyphen_for_elem_in_basfile(basname)
 if(add_path) call add_mokit_path_to_genbas(basname)
end subroutine record_gen_basis_in_gjf

! add '-' symbol before elements, in a .bas file
subroutine add_hyphen_for_elem_in_basfile(basname)
 implicit none
 integer :: i, j, nbat1, nbat2, fid, fid1, RENAME
 character(len=7) :: str
 character(len=240) :: buf0, buf, basname1
 character(len=240), intent(in) :: basname

 basname1 = TRIM(basname)//'.t'
 open(newunit=fid,file=TRIM(basname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(basname1),status='replace')
 buf0 = '****'

 ! deal with the basis set data
 do while(.true.)
  read(fid,'(A)') buf

  if(buf0(1:4) == '****') then
   if(LEN_TRIM(buf) == 0) exit
   call add_hyphen_for_elem_in_buf(buf)
   write(fid1,'(A)') TRIM(buf)
  else
   write(fid1,'(A)') TRIM(buf)
  end if

  buf0 = buf
 end do ! for i

 write(fid1,'(/)',advance='no')

 ! deal with the pseudo potential data
 do while(.true.)
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit

  call add_hyphen_for_elem_in_buf(buf)
  write(fid1,'(A)') TRIM(buf)

  read(fid,'(A)') buf
  write(fid1,'(A)') TRIM(buf)
  read(buf,*,iostat=i) str, nbat1

  if(i == 0) then ! ECP/PP data, not name
   nbat1 = nbat1 + 1
   do i = 1, nbat1, 1
    read(fid,'(A)') buf
    write(fid1,'(A)') TRIM(buf)
    read(fid,'(A)') buf
    write(fid1,'(A)') TRIM(buf)
    read(buf,*) nbat2
    do j = 1, nbat2, 1
     read(fid,'(A)') buf
     write(fid1,'(A)') TRIM(buf)
    end do ! for j
   end do ! for i
  end if
 end do ! for while

 write(fid1,'(/)')
 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(basname1), TRIM(basname))
end subroutine add_hyphen_for_elem_in_basfile

function get_mokit_root() result(mokit_root)
 implicit none
 integer :: i, fid
 character(len=240) :: home, mokit_root !, buf
 mokit_root = ' '
 call getenv('MOKIT_ROOT', mokit_root)
 if (len_trim(mokit_root) < 1) then
  call getenv('HOME', home)
  open(newunit=fid,file=TRIM(home)//'/.mokitrc',status='old',position='rewind')
  read(fid,'(A)',iostat=i) mokit_root
  if (len_trim(mokit_root) < 1) then
    write(6,'(/,A)') 'ERROR in subroutine get_mokit_root: invalid MOKIT_ROOT'
    stop
  end if
  close(fid)
 end if
end function get_mokit_root

! add MOKIT_ROOT path into basis sets like ANO-RCC-VDZP, DKH-def2-SVP in file
! basname because MOKIT has these basis sets in $MOKIT_ROOT/basis/
subroutine add_mokit_path_to_genbas(basname)
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=12) :: sbuf
 character(len=240) :: buf, mokit_root, get_mokit_root, basname1
 character(len=240), intent(in) :: basname

 !mokit_root = ' '
 !call getenv('MOKIT_ROOT', mokit_root)
 mokit_root = get_mokit_root()
 basname1 = TRIM(basname)//'.t'

 open(newunit=fid,file=TRIM(basname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(basname1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) exit
  sbuf = TRIM(buf(1:12))
  call upper(sbuf)
  if(sbuf(1:6)=='PCSSEG' .or. sbuf(1:7)=='ANO-RCC' .or. sbuf(1:8)=='DKH-DEF2' &
     .or. sbuf(1:9)=='ZORA-DEF2' .or. sbuf(1:11)=='MA-DKH-DEF2' .or. &
     sbuf(1:12)=='MA-ZORA-DEF2') then
   buf = '@'//TRIM(mokit_root)//'/basis/'//TRIM(buf)
  end if
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 write(fid1,'(/)')
 close(fid1)
 i = RENAME(TRIM(basname1), TRIM(basname))
end subroutine add_mokit_path_to_genbas

subroutine create_basfile(basfile, basis)
 implicit none
 integer :: fid
 character(len=240), intent(in) :: basfile
 character(len=*), intent(in) :: basis

 open(newunit=fid,file=TRIM(basfile),status='replace')
 write(fid,'(A)') TRIM(basis)
 write(fid,'(/)')
 close(fid)
 call add_mokit_path_to_genbas(basfile)
end subroutine create_basfile

! add '-' symbol for each element, in a buf
subroutine add_hyphen_for_elem_in_buf(buf)
 implicit none
 integer :: i, j
 integer, parameter :: max_nelem = 21
 character(len=3) :: str(max_nelem)
 character(len=240), intent(inout) :: buf

 do i = 1, max_nelem, 1
  read(buf,*,iostat=j) str(1:i)
  if(j /= 0) exit
 end do ! for i

 i = i - 1
 if(TRIM(str(i)) == '0') i = i - 1
 forall(j=1:i, str(j)(1:1)/='-') str(j) = '-'//TRIM(str(j))

 buf = TRIM(str(1))
 do j = 2, i, 1
  buf = TRIM(buf)//' '//TRIM(str(j))
 end do ! for i

 buf = TRIM(buf)//' 0'
end subroutine add_hyphen_for_elem_in_buf

! copy mixed/user-defined basis set content from file basname to gjfname
subroutine copy_gen_basis_bas2gjf(basname, gjfname)
 implicit none
 integer :: i, nblank, fid1, fid2
 character(len=240) :: buf
 character(len=240), intent(in) :: basname, gjfname

 if(LEN_TRIM(basname) == 0) return

 open(newunit=fid1,file=TRIM(gjfname),status='old',position='rewind')
 nblank = 0
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) nblank = nblank + 1
  if(nblank == 3) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine copy_gen_basis_bas2gjf: incomplete&
                 & file: '//TRIM(gjfname)
  close(fid1)
  stop
 end if

 open(newunit=fid2,file=TRIM(basname),status='old',position='rewind')
 do while(.true.)
  read(fid2,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid1)
 close(fid2)
end subroutine copy_gen_basis_bas2gjf

! read the version of dispersion correction from a .gjf file
subroutine read_disp_ver_from_gjf(gjfname, itype)
 implicit none
 integer :: i, fid
 integer, intent(out) :: itype
 character(len=240) :: buf
 character(len=240), intent(in) :: gjfname

 itype = 0
 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:1) == '#') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_disp_ver_from_gjf: no '#' symbol&
                   & found in file "//TRIM(gjfname)
  stop
 end if

 call lower(buf)
 if(index(buf,'em=gd3bj')>0 .or. index(buf,'empiricaldispersion=gd3bj')>0) then
  itype = 2
 else if(index(buf,'em=gd3')>0 .or. index(buf,'empiricaldispersion=gd3')>0) then
  itype = 1
 end if
end subroutine read_disp_ver_from_gjf

! print Fock operator coupling coefficients for ROGVB when nopen>=3
subroutine prt_gvb_couple_coeff(fid, nopen)
 implicit none
 integer :: i, ia
 integer, intent(in) :: fid, nopen
 character(len=3), allocatable :: f(:), alpha(:)
 character(len=4), allocatable :: beta(:)

 allocate(f(nopen))
 f = '0.5'
 ia = nopen*(nopen+3)/2
 allocate(alpha(ia))
 alpha = '0.5'
 forall(i = 1:nopen) alpha(i*(i+1)/2) = '1.0'
 allocate(beta(ia))
 beta = '-0.5'
 write(fid,'(A,10(A1,A3))') '  F(1)=1.0', (',',f(i),i=1,nopen)
 write(fid,'(A,15(A1,A3))') '  ALPHA(1)=2.0', (',',alpha(i),i=1,ia)
 write(fid,'(A,12(A1,A4))') '  BETA(1)=-1.0', (',',beta(i),i=1,ia)

 deallocate(f, alpha, beta)
end subroutine prt_gvb_couple_coeff

! copy GVB CI coefficients (or called pair coefficients) from a .dat file into
! another one 
subroutine copy_and_add_pair_coeff(addH_dat, datname, nopen)
 implicit none
 integer :: i, j, npair, fid1, fid2, fid3, RENAME
 integer, intent(in) :: nopen
 character(len=240) :: buf, new_dat
 character(len=240), intent(in) :: addH_dat, datname

 i = index(addH_dat, '.dat', back=.true.)
 new_dat = addH_dat(1:i-1)//'.t'
 open(newunit=fid1,file=TRIM(addH_dat),status='old',position='rewind')
 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(2:6) == '$DATA') exit
 end do ! for while

 open(newunit=fid2,file=TRIM(new_dat),status='replace')
 write(fid2,'(A)') ' $DATA'

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(2:5) == '$VEC') exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 open(newunit=fid3,file=TRIM(datname),status='old',position='rewind')
 do while(.true.)
  read(fid3,'(A)') buf
  if(buf(2:5) == '$SCF') exit
 end do ! for while

 i = index(buf, 'CICOEF')
 if(i > 0) then
  write(fid2,'(A)') TRIM(buf)
  npair = 1
 else
  npair = 0
 end if

 do while(.true.)
  read(fid3,'(A)') buf
  i = index(buf, 'CICOEF')
  if(i > 0) then
   j = index(buf, '$END')
   if(j > 0) then
    buf(j:j+3) = '    '
    write(fid2,'(A)') TRIM(buf)
    exit
   else
    write(fid2,'(A)') TRIM(buf)
   end if
   npair = npair + 1
  else
   exit
  end if
 end do ! for while

 close(fid3)
 do i = 1, nopen, 1
  j = 2*(npair+i) - 1
  write(fid2,'(3X,A,I3,A)') 'CICOEF(',j,')= 0.7071067811865476,-0.7071067811865476'
 end do ! for i
 write(fid2,'(A,/,A)') ' $END', ' $VEC'

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(TRIM(new_dat), TRIM(addH_dat))
end subroutine copy_and_add_pair_coeff

subroutine add_force_key2py_script(mem, pyname)
 implicit none
 integer :: fid
 integer, intent(in) :: mem ! GB
 character(len=240), intent(in) :: pyname

 open(newunit=fid,file=TRIM(pyname),status='old',position='append')
 write(fid,'(A)') 'from pyscf import grad'
 write(fid,'(A)') 'mcg = mc.Gradients()'
 write(fid,'(A,I0,A)') 'mcg.max_memory = ',mem*1000,' # MB'
 write(fid,'(A)') 'mcg.kernel()'
 close(fid)
end subroutine add_force_key2py_script

! standardize a set of elements, e.g. he -> He
subroutine standardize_elem(natom, elem)
 implicit none
 integer :: i, j
 integer, intent(in) :: natom
 character(len=2), intent(inout) :: elem(natom)

 do i = 1, natom, 1
  elem(i) = ADJUSTL(elem(i))

  j = IACHAR(elem(i)(1:1))
  if(j>96 .and. j<123) elem(i)(1:1) = ACHAR(j-32)

  j = IACHAR(elem(i)(2:2))
  if(j>64 .and. j<91) elem(i)(2:2) = ACHAR(j+32)
 end do ! for i
end subroutine standardize_elem

