! written by jxzou at 20200809: wrappers of utilities

module util_wrapper
 implicit none

contains

! wrapper of the Gaussian utility formchk
subroutine formchk(chkname, fchname)
 implicit none
 integer :: i, system
 integer, parameter :: iout = 6
 character(len=240), intent(in) :: chkname
 character(len=240), optional :: fchname

 if(.not. present(fchname)) then
  i = index(chkname, '.chk', back=.true.)
  fchname = chkname(1:i-1)//'.fch'
 end if

#ifdef _WIN32
 i = system('formchk '//TRIM(chkname)//' '//TRIM(fchname)//' > NUL')
#else
 i = system('formchk '//TRIM(chkname)//' '//TRIM(fchname)//' > /dev/null')
#endif

 if(i /= 0) then
  write(iout,'(/,A)') 'ERROR in subroutine formchk: call Gaussian utility formchk failed.'
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
 integer, parameter :: iout = 6
 character(len=240), intent(in) :: fchname
 character(len=240), optional :: chkname

 if(.not. present(chkname)) then
  i = index(fchname, '.fch', back=.true.)
  chkname = fchname(1:i-1)//'.chk'
 end if

#ifdef _WIN32
 i = system('unfchk '//TRIM(fchname)//' '//TRIM(chkname)//' > NUL')
#else
 i = system('unfchk '//TRIM(fchname)//' '//TRIM(chkname)//' > /dev/null')
#endif

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine formchk: call Gaussian utility unfchk failed.'
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
 integer, parameter :: iout = 6
 character(len=240), intent(in) :: gbwname
 character(len=240), optional :: mklname

 k = index(gbwname, '.gbw', back=.true.)
#ifdef _WIN32
 i = system('orca_2mkl '//gbwname(1:k-1)//' -mkl > NUL')
#else
 i = system('orca_2mkl '//gbwname(1:k-1)//' -mkl > /dev/null')
#endif

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine gbw2mkl: call ORCA utility orca_2mkl failed.'
  write(iout,'(A)') 'The file '//TRIM(gbwname)//' may be incomplete, or ORCA utility&
                   & orca_2mkl does not exist.'
  stop
 end if

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
 integer, parameter :: iout = 6
 character(len=240), intent(in) :: mklname
 character(len=240), optional :: gbwname

 k = index(mklname, '.mkl', back=.true.)
#ifdef _WIN32
 i = system('orca_2mkl '//mklname(1:k-1)//' -gbw > NUL')
#else
 i = system('orca_2mkl '//mklname(1:k-1)//' -gbw > /dev/null')
#endif

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine mkl2gbw: call ORCA utility orca_2mkl failed.'
  write(iout,'(A)') 'The file '//TRIM(mklname)//' may be incomplete, or ORCA utility&
                   & orca_2mkl does not exist.'
  stop
 end if

 if(present(gbwname)) then
  if(TRIM(gbwname) /= mklname(1:k-1)//'.gbw') then
   i = RENAME(mklname(1:k-1)//'.gbw', TRIM(gbwname))
  end if
 end if

 return
end subroutine mkl2gbw

! wrapper of the utility fch2inp
subroutine fch2inp_wrap(fchname, gvb, npair, nopen)
 implicit none
 integer :: i, system
 integer, parameter :: iout = 6
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

end module util_wrapper

