! written by jxzou at 20200809: wrappers of utilities

module util_wrapper
 implicit none
 integer, parameter :: iout = 6
contains

! wrapper of the Gaussian utility formchk
subroutine formchk(chkname, fchname)
 implicit none
 integer :: i, system
 character(len=240) :: fchname0
 character(len=500) :: buf
 character(len=240), intent(in) :: chkname
 character(len=240), optional :: fchname
 logical :: alive

 inquire(file=TRIM(chkname), exist=alive)
 if(.not. alive) then
  write(iout,'(A)') 'ERROR in subroutine formchk: file does not exist!'
  write(iout,'(A)') 'chkname='//TRIM(chkname)
  stop
 end if

 if(present(fchname)) then
  fchname0 = fchname
 else
  i = index(chkname, '.chk', back=.true.)
  if(i == 0) then
   write(iout,'(A)') 'ERROR in subroutine formchk: .chk suffix not found!'
   write(iout,'(A)') 'chkname='//TRIM(chkname)
   stop
  end if
  fchname0 = chkname(1:i-1)//'.fch'
 end if
 buf = 'formchk '//TRIM(chkname)//' '//TRIM(fchname0)

#ifdef _WIN32
 i = system(TRIM(buf)//' > NUL')
#else
 i = system(TRIM(buf)//' > /dev/null')
#endif

 if(i /= 0) then
  write(iout,'(/,A)') 'ERROR in subroutine formchk: failed to call Gaussian utility formchk.'
  write(iout,'(A)') 'The file '//TRIM(chkname)//' may be incomplete, or Gaussian&
                   & utility formchk does not exist.'
  stop
 end if

 return
end subroutine formchk

! wrapper of the Gaussian utility unfchk
subroutine unfchk(fchname, chkname)
 implicit none
 integer :: i, system
 character(len=240) :: chkname0
 character(len=500) :: buf
 character(len=240), intent(in) :: fchname
 character(len=240), optional :: chkname

 if(present(chkname)) then
  chkname0 = chkname
 else
  i = index(fchname, '.fch', back=.true.)
  chkname0 = fchname(1:i-1)//'.chk'
 end if
 buf = 'unfchk '//TRIM(fchname)//' '//TRIM(chkname0)

#ifdef _WIN32
 i = system(TRIM(buf)//' > NUL')
#else
 i = system(TRIM(buf)//' > /dev/null')
#endif

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine formchk: failed to call Gaussian utility unfchk.'
  write(iout,'(A)') 'The file '//TRIM(fchname)//' may be incomplete, or Gaussian&
                   & utility unfchk does not exist.'
  stop
 end if

 return
end subroutine unfchk

! wrapper of the ORCA utility orca_2mkl, only .gbw -> .mkl
subroutine gbw2mkl(gbwname, mklname)
 implicit none
 integer :: i, k, system, RENAME
 character(len=240), intent(in) :: gbwname
 character(len=240), optional :: mklname
 logical :: alive

 inquire(file=TRIM(gbwname), exist=alive)
 if(.not. alive) then
  write(iout,'(A)') 'ERROR in subroutine gbw2mkl: file does not exist!'
  write(iout,'(A)') 'gbwname='//TRIM(gbwname)
  stop
 end if

 k = index(gbwname, '.gbw', back=.true.)
#ifdef _WIN32
 i = system('orca_2mkl '//gbwname(1:k-1)//' -mkl > NUL')
#else
 i = system('orca_2mkl '//gbwname(1:k-1)//' -mkl > /dev/null')
#endif

 if(i /= 0) call prt_orca_2mkl_error(gbwname)

 if(present(mklname)) then
  if(TRIM(mklname) /= gbwname(1:k-1)//'.mkl') then
   i = RENAME(gbwname(1:k-1)//'.mkl', TRIM(mklname))
  end if
 end if

 return
end subroutine gbw2mkl

! wrapper of the ORCA utility orca_2mkl, only .mkl -> .gbw
subroutine mkl2gbw(mklname, gbwname)
 implicit none
 integer :: i, k, system, RENAME
 character(len=240), intent(in) :: mklname
 character(len=240), optional :: gbwname

 k = index(mklname, '.mkl', back=.true.)
#ifdef _WIN32
 i = system('orca_2mkl '//mklname(1:k-1)//' -gbw > NUL')
#else
 i = system('orca_2mkl '//mklname(1:k-1)//' -gbw > /dev/null')
#endif

 if(i /= 0) call prt_orca_2mkl_error(mklname)

 if(present(gbwname)) then
  if(TRIM(gbwname) /= mklname(1:k-1)//'.gbw') then
   i = RENAME(mklname(1:k-1)//'.gbw', TRIM(gbwname))
  end if
 end if
end subroutine mkl2gbw

subroutine prt_orca_2mkl_error(fname)
 implicit none
 character(len=240), intent(in) :: fname

 write(6,'(/,A)') 'ERROR: failed to call ORCA utility orca_2mkl. Three&
                 & possible reasons:'
 write(6,'(A)') '(1) Your ORCA environment variables are incorrect.'
 write(6,'(A)') '(2) ORCA utility orca_2mkl does not exist.'
 write(6,'(A)') '(3) The file '//TRIM(fname)//' may be incomplete.'
 stop
end subroutine prt_orca_2mkl_error

! wrapper of the utility fch2psi
subroutine fch2psi_wrap(fchname)
 implicit none
 integer :: i, system
 character(len=240), intent(in) :: fchname

#ifdef _WIN32
 i = system('fch2psi '//TRIM(fchname)//' > NUL')
#else
 i = system('fch2psi '//TRIM(fchname)//' > /dev/null')
#endif

 if(i /= 0) then
  write(iout,'(/,A)') 'ERROR in subroutine fch2psi_wrap: failed to call utility&
                     & fch2psi. Two possible reasons:'
  write(iout,'(A)') '(1) You might forget to compile the utility fch2psi.'
  write(iout,'(A)') '(2) The file '//TRIM(fchname)//' may be incomplete.'
  write(iout,'(A)') '(3) There might exist a bug in the utility fch2psi.'
  stop
 end if

 return
end subroutine fch2psi_wrap

! wrapper of the utility fch2inp
subroutine fch2inp_wrap(fchname, gvb, npair, nopen)
 implicit none
 integer :: i, system
 integer, intent(in) :: npair, nopen
 character(len=240) :: buf = ' ', buf2 = ' '
 character(len=240), intent(in) :: fchname
 logical, intent(in) :: gvb

 if(gvb) then
  if(npair<0 .or. nopen<0) then
   write(iout,'(A)') 'ERROR in subroutine fch2inp_wrap: gvb=.True. but npair<0&
                    & and/or nopen<0.'
   write(iout,'(2(A,I0))') 'npair = ', npair, ', nopen = ', nopen
   stop
  end if

  write(buf,'(A,I0)') 'fch2inp '//TRIM(fchname)//' -gvb ', npair
  buf = ADJUSTL(buf)
  if(nopen == 0) then
   buf2 = buf
  else ! nopen > 0
   write(buf2,'(A,I0)') ' -open ',nopen
   buf2 = ADJUSTL(buf2)
   buf2 = TRIM(buf)//TRIM(buf2)
  end if

#ifdef _WIN32
  i = system(TRIM(buf2)//' > NUL')
#else
  i = system(TRIM(buf2)//' > /dev/null')
#endif

 else ! R(O)HF, UHF, CAS
#ifdef _WIN32
  i = system('fch2inp '//TRIM(fchname)//' > NUL')
#else
  i = system('fch2inp '//TRIM(fchname)//' > /dev/null')
#endif
 end if

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine fch2inp_wrap: failed to call utility fch2inp.'
  write(iout,'(A)') 'Filename = '//TRIM(fchname)
  write(iout,'(2(A,I0),A,L1)') 'npair= ',npair,', nopen= ',nopen,', gvb= ',gvb
  stop
 end if

 return
end subroutine fch2inp_wrap

subroutine mkl2fch_wrap(mklname, fchname, prt_no)
 implicit none
 integer :: i, system
 character(len=240), intent(in) :: mklname, fchname
 character(len=500) :: buf
 logical, intent(in) :: prt_no

 buf = 'mkl2fch '//TRIM(mklname)//' '//TRIM(fchname)
 if(prt_no) buf = TRIM(buf)//' -no'

#ifdef _WIN32
 i = system(TRIM(buf)//' > NUL')
#else
 i = system(TRIM(buf)//' > /dev/null')
#endif

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine mkl2fch_wrap: failed to call utility mkl2fch.'
  write(iout,'(A)') 'mklname='//TRIM(mklname)//', fchname='//TRIM(fchname)
  write(iout,'(A,L1)') 'prt_no=', prt_no
  stop
 end if
 return
end subroutine mkl2fch_wrap

subroutine fch2mkl_wrap(fchname, mklname)
 implicit none
 integer :: i, system, RENAME
 character(len=240) :: mklname1
 character(len=240), intent(in) :: fchname
 character(len=240), optional :: mklname

#ifdef _WIN32
 i = system('fch2mkl '//TRIM(fchname)//' > NUL')
#else
 i = system('fch2mkl '//TRIM(fchname)//' > /dev/null')
#endif

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine fch2mkl_wrap: failed to call utility fch2mkl.'
  write(iout,'(A)') 'fchname='//TRIM(fchname)
  stop
 end if

 if(present(mklname)) then
  i = index(fchname, '.fch')
  mklname1 = fchname(1:i-1)//'_o.mkl'
  if(TRIM(mklname) /= TRIM(mklname1)) then
   i = RENAME(TRIM(mklname1), TRIM(mklname))
  end if
 end if
end subroutine fch2mkl_wrap

subroutine fch2gbw(fchname, gbwname)
 implicit none
 integer :: i
 character(len=240) :: mklname, inpname
 character(len=240), intent(in) :: fchname
 character(len=240), optional :: gbwname

 i = index(fchname, '.fch', back=.true.)
 inpname = fchname(1:i-1)//'_o.inp'

 if(present(gbwname)) then
  i = index(gbwname, '.gbw', back=.true.)
  mklname = gbwname(1:i-1)//'.mkl'
 else
  mklname = fchname(1:i-1)//'.mkl'
 end if

 call fch2mkl_wrap(fchname, mklname)
 open(newunit=i,file=TRIM(inpname))
 close(unit=i,status='delete')

 call mkl2gbw(mklname)
 open(newunit=i,file=TRIM(mklname))
 close(unit=i,status='delete')
end subroutine fch2gbw

subroutine chk2gbw(chkname)
 implicit none
 integer :: i
 character(len=240) :: fchname, inpname, mklname, gbwname
 character(len=240), intent(in) :: chkname

 i = index(chkname, '.chk')
 fchname = chkname(1:i-1)//'.fch'
 inpname = chkname(1:i-1)//'_o.inp'
 mklname = chkname(1:i-1)//'_o.mkl'
 gbwname = chkname(1:i-1)//'.gbw'
 call formchk(chkname)

 call fch2mkl_wrap(fchname)
 open(newunit=i,file=TRIM(fchname),status='old')
 close(unit=i,status='delete')
 open(newunit=i,file=TRIM(inpname),status='old')
 close(unit=i,status='delete')

 call mkl2gbw(mklname, gbwname)
 open(newunit=i,file=TRIM(mklname),status='old')
 close(unit=i,status='delete')
end subroutine chk2gbw

subroutine fch_u2r_wrap(fchname)
 implicit none
 integer :: i, system
 character(len=240), intent(in) :: fchname

 i = system('fch_u2r '//TRIM(fchname))
 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine fch_u2r_wrap: failed to call utility&
                   & fch_u2r.'
  write(iout,'(A)') 'fchname='//TRIM(fchname)
  stop
 end if
end subroutine fch_u2r_wrap

subroutine fch2dal_wrap(fchname, dalname)
 implicit none
 integer :: i, system, RENAME
 character(len=240) :: molname, dalname1, molname1
 character(len=240), intent(in) :: fchname
 character(len=240), optional :: dalname

 i = system('fch2dal '//TRIM(fchname))
 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine fch2dal_wrap: failed to call utility&
                   & fch2dal.'
  write(iout,'(A)') 'fchname='//TRIM(fchname)
  stop
 end if

 if(present(dalname)) then
  i = index(dalname, '.dal')
  molname = dalname(1:i-1)//'.mol'
  i = index(fchname, '.fch')
  dalname1 = fchname(1:i-1)//'.dal'
  molname1 = fchname(1:i-1)//'.mol'
  i = RENAME(TRIM(dalname1), TRIM(dalname))
  i = RENAME(TRIM(molname1), TRIM(molname))
 end if

end subroutine fch2dal_wrap

subroutine bas_fch2py_wrap(fchname, pyname)
 implicit none
 integer :: i, system, RENAME
 character(len=240) :: pyname0
 character(len=240), intent(in) :: fchname
 character(len=240), optional :: pyname

#ifdef _WIN32
 i = system('bas_fch2py '//TRIM(fchname)//' > NUL')
#else
 i = system('bas_fch2py '//TRIM(fchname)//' > /dev/null')
#endif

 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine bas_fch2py_wrap: failed to call utility&
                & bas_fch2py.'
  write(6,'(A)') 'fchname='//TRIM(fchname)
  stop
 end if

 if(present(pyname)) then
  i = index(fchname, '.fch', back=.true.)
  pyname0 = fchname(1:i-1)//'.py'
  i = RENAME(TRIM(pyname0), TRIM(pyname))
 end if
end subroutine bas_fch2py_wrap

end module util_wrapper

