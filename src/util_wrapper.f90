! written by jxzou at 20200809: wrappers of utilities

module util_wrapper
 implicit none
contains

! wrapper of the Gaussian utility formchk
subroutine formchk(chkname, fchname)
 implicit none
 integer :: i, SYSTEM
 character(len=240) :: fchname0
 character(len=500) :: buf
 character(len=240), intent(in) :: chkname
 character(len=240), intent(in), optional :: fchname
 logical :: alive

 inquire(file=TRIM(chkname), exist=alive)
 if(.not. alive) then
  write(6,'(/,A)') 'ERROR in subroutine formchk: file does not exist!'
  write(6,'(A)') 'chkname='//TRIM(chkname)
  stop
 end if

 if(PRESENT(fchname)) then
  fchname0 = fchname
 else
  i = INDEX(chkname, '.chk', back=.true.)
  fchname0 = chkname(1:i-1)//'.fch'
 end if
 buf = 'formchk '//TRIM(chkname)//' '//TRIM(fchname0)

#ifdef _WIN32
 i = SYSTEM(TRIM(buf)//' > NUL')
#else
 i = SYSTEM(TRIM(buf)//' > /dev/null')
#endif

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine formchk: failed to call Gaussian utilit&
                   &y formchk.'
  write(6,'(A)') 'The file '//TRIM(chkname)//' may be problematic, or the utili&
                 &ty formchk does not exist.'
  stop
 end if
end subroutine formchk

! wrapper of the Gaussian utility unfchk
subroutine unfchk(fchname, chkname)
 implicit none
 integer :: i, SYSTEM
 character(len=240) :: chkname0
 character(len=500) :: buf
 character(len=240), intent(in) :: fchname
 character(len=240), intent(in), optional :: chkname

 if(PRESENT(chkname)) then
  chkname0 = chkname
 else
  i = INDEX(fchname, '.fch', back=.true.)
  chkname0 = fchname(1:i-1)//'.chk'
 end if
 buf = 'unfchk '//TRIM(fchname)//' '//TRIM(chkname0)

#ifdef _WIN32
 i = SYSTEM(TRIM(buf)//' > NUL')
#else
 i = SYSTEM(TRIM(buf)//' > /dev/null')
#endif

 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine formchk: failed to call Gaussian utility unfchk.'
  write(6,'(A)') 'The file '//TRIM(fchname)//' may be incomplete, or Gaussian&
                   & utility unfchk does not exist.'
  stop
 end if
end subroutine unfchk

! wrapper of the ORCA utility orca_2mkl, only .gbw -> .mkl
subroutine gbw2mkl(gbwname, mklname)
 implicit none
 integer :: i, k, SYSTEM, RENAME
 character(len=240), intent(in) :: gbwname
 character(len=240), intent(in), optional :: mklname
 logical :: alive

 inquire(file=TRIM(gbwname), exist=alive)
 if(.not. alive) then
  write(6,'(/,A)') 'ERROR in subroutine gbw2mkl: file does not exist!'
  write(6,'(A)') 'Input gbwname='//TRIM(gbwname)
  stop
 end if

 k = INDEX(gbwname, '.gbw', back=.true.)
#ifdef _WIN32
 i = SYSTEM('orca_2mkl '//gbwname(1:k-1)//' -mkl > NUL')
#else
 i = SYSTEM('orca_2mkl '//gbwname(1:k-1)//' -mkl > /dev/null')
#endif

 if(i /= 0) call prt_orca_2mkl_error(gbwname)

 if(PRESENT(mklname)) then
  if(TRIM(mklname) /= gbwname(1:k-1)//'.mkl') then
   i = RENAME(gbwname(1:k-1)//'.mkl', TRIM(mklname))
  end if
 end if
end subroutine gbw2mkl

! wrapper of the ORCA utility orca_2mkl, only .mkl -> .gbw
subroutine mkl2gbw(mklname, gbwname)
 implicit none
 integer :: i, k, SYSTEM, RENAME
 character(len=240), intent(in) :: mklname
 character(len=240), intent(in), optional :: gbwname

 k = INDEX(mklname, '.mkl', back=.true.)
#ifdef _WIN32
 i = SYSTEM('orca_2mkl '//mklname(1:k-1)//' -gbw > NUL')
#else
 i = SYSTEM('orca_2mkl '//mklname(1:k-1)//' -gbw > /dev/null')
#endif

 if(i /= 0) call prt_orca_2mkl_error(mklname)

 if(PRESENT(gbwname)) then
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
subroutine fch2psi_wrap(fchname, inpname)
 implicit none
 integer :: i, SYSTEM, RENAME
 character(len=240) :: inpname1
 character(len=240), intent(in) :: fchname
 character(len=240), intent(in), optional :: inpname

 i = SYSTEM('fch2psi '//TRIM(fchname))
 if(i /= 0) call prt_call_util_error('fch2psi', fchname)
 if(PRESENT(inpname)) then
  i = INDEX(fchname, '.fch')
  inpname1 = fchname(1:i-1)//'_psi.inp'
  if(TRIM(inpname) /= TRIM(inpname1)) then
   i = RENAME(TRIM(inpname1), TRIM(inpname))
  end if
 end if
end subroutine fch2psi_wrap

! wrapper of the utility fch2inp
subroutine fch2inp_wrap(fchname, gvb, npair, nopen, prt)
 implicit none
 integer :: i, SYSTEM
 integer, intent(in) :: npair, nopen
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
 logical :: alive
 logical, intent(in) :: gvb
 logical, intent(in), optional :: prt

 buf = ' '

 if(gvb) then
  if(npair<0 .or. nopen<0) then
   write(6,'(/,A)') 'ERROR in subroutine fch2inp_wrap: gvb=.True. but npair<0 a&
                    &nd/or nopen<0.'
   write(6,'(2(A,I0))') 'npair = ', npair, ', nopen = ', nopen
   stop
  end if

  if(nopen == 0) then
   write(buf,'(A,I0)') 'fch2inp '//TRIM(fchname)//' -gvb ', npair
  else
   write(buf,'(2(A,I0))') 'fch2inp '//TRIM(fchname)//' -gvb ',npair,' -open ',nopen
  end if
 else ! for R(O)HF, UHF, CAS
  write(buf,'(A)') 'fch2inp '//TRIM(fchname)
 end if

 alive = .true.
 if(present(prt)) then
  if(.not. prt) alive = .false.
 end if

 if(alive) then
  write(6,'(A)') '$'//TRIM(buf)
  i = SYSTEM(TRIM(buf))
 else
#ifdef _WIN32
  i = SYSTEM(TRIM(buf)//' > NUL')
#else
  i = SYSTEM(TRIM(buf)//' > /dev/null')
#endif
 end if

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine fch2inp_wrap: failed to call utility fc&
                   &h2inp.'
  write(6,'(A)') 'Filename = '//TRIM(fchname)
  write(6,'(2(A,I0),A,L1)') 'npair= ',npair,', nopen= ',nopen,', gvb= ',gvb
  stop
 end if
end subroutine fch2inp_wrap

subroutine mkl2fch_wrap(mklname, fchname, ino, irel)
 implicit none
 integer :: i, SYSTEM
 integer, intent(in), optional :: ino, irel
 character(len=240), intent(in) :: mklname, fchname
 character(len=500) :: buf

 buf = 'mkl2fch '//TRIM(mklname)//' '//TRIM(fchname)

 if(PRESENT(ino)) then
  select case(ino)
  case(0) ! do nothing
  case(1)
   buf = TRIM(buf)//' -no'
  case(2)
   buf = TRIM(buf)//' -nso'
  case default
   write(6,'(/,A,I0)') 'ERROR in subroutine mkl2fch_wrap: invalid ino=', ino
   write(6,'(A)') 'Only {0,1,2} are allowed.'
   stop
  end select
 end if

 if(PRESENT(irel)) then
  select case(irel)
  case(-3) ! sfX2C
   buf = TRIM(buf)//' -sfx2c'
  case(-1) ! do nothing
  case(2)  ! DKH2
   buf = TRIM(buf)//' -dkh2'
  case default
   write(6,'(/,A,I0)') 'ERROR in subroutine mkl2fch_wrap: invalid irel=',irel
   write(6,'(A)') 'Only {-3,-1,2} are allowed.'
   stop
  end select
 end if

#ifdef _WIN32
 i = SYSTEM(TRIM(buf)//' > NUL')
#else
 i = SYSTEM(TRIM(buf)//' > /dev/null')
#endif

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine mkl2fch_wrap: failed to call utility mk&
                   &l2fch.'
  write(6,'(A)') 'mklname = '//TRIM(mklname)
  write(6,'(A)') 'fchname = '//TRIM(fchname)
  if(PRESENT(ino)) write(6,'(A,I0)') 'ino = ', ino
  if(PRESENT(irel)) write(6,'(A,I0)') 'irel = ', irel
  stop
 end if
end subroutine mkl2fch_wrap

subroutine fch2mkl_wrap(fchname, mklname)
 implicit none
 integer :: i, SYSTEM, RENAME
 character(len=240) :: mklname1
 character(len=240), intent(in) :: fchname
 character(len=240), intent(in), optional :: mklname

 i = SYSTEM('fch2mkl '//TRIM(fchname))
 if(i /= 0) call prt_call_util_error('fch2mkl', fchname)

 if(PRESENT(mklname)) then
  i = INDEX(fchname, '.fch')
  mklname1 = fchname(1:i-1)//'_o.mkl'
  if(TRIM(mklname) /= TRIM(mklname1)) then
   i = RENAME(TRIM(mklname1), TRIM(mklname))
  end if
 end if
end subroutine fch2mkl_wrap

subroutine chk2gbw(chkname)
 implicit none
 integer :: i
 character(len=240) :: fchname, inpname, mklname, gbwname
 character(len=240), intent(in) :: chkname

 i = INDEX(chkname, '.chk')
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

subroutine fch_u2r_wrap(fchname, new_fch)
 implicit none
 integer :: i, SYSTEM, RENAME
 character(len=240) :: rfch
 character(len=240), intent(in) :: fchname
 character(len=240), intent(in), optional :: new_fch

 i = SYSTEM('fch_u2r '//TRIM(fchname))
 if(i /= 0) call prt_call_util_error('fch_u2r', fchname)

 if(PRESENT(new_fch)) then
  i = INDEX(fchname, '.fch', back=.true.)
  rfch = fchname(1:i-1)//'_r.fch'
  i = RENAME(TRIM(rfch), TRIM(new_fch))
 end if
end subroutine fch_u2r_wrap

subroutine fch2dal_wrap(fchname, dalname)
 implicit none
 integer :: i, SYSTEM, RENAME
 character(len=240) :: molname, dalname1, molname1
 character(len=240), intent(in) :: fchname
 character(len=240), intent(in), optional :: dalname

 i = SYSTEM('fch2dal '//TRIM(fchname))
 if(i /= 0) call prt_call_util_error('fch2dal', fchname)

 if(PRESENT(dalname)) then
  i = INDEX(dalname, '.dal')
  molname = dalname(1:i-1)//'.mol'
  i = INDEX(fchname, '.fch')
  dalname1 = fchname(1:i-1)//'.dal'
  molname1 = fchname(1:i-1)//'.mol'
  i = RENAME(TRIM(dalname1), TRIM(dalname))
  i = RENAME(TRIM(molname1), TRIM(molname))
 end if

end subroutine fch2dal_wrap

! wrapper of utility fch2qchem
subroutine fch2qchem_wrap(fchname, npair, inpname)
 implicit none
 integer :: i, SYSTEM, RENAME
 integer, intent(in) :: npair
 character(len=240) :: inpname0, dirname
 character(len=240), intent(in) :: fchname
 character(len=240), intent(in), optional :: inpname
 character(len=260) :: buf
 character(len=480) :: scr_dir, scr_dir0

 buf = 'fch2qchem '//TRIM(fchname)
 if(npair > 0) write(buf,'(A,I0)') TRIM(buf)//' -gvb ', npair

#ifdef _WIN32
 i = SYSTEM(TRIM(buf)//' > NUL')
#else
 i = SYSTEM(TRIM(buf)//' > /dev/null')
#endif

 if(i /= 0) call prt_call_util_error('fch2qchem', fchname)

 if(PRESENT(inpname)) then
  dirname = ' '
  call getenv('QCSCRATCH', dirname)

  i = INDEX(fchname, '.fch', back=.true.)
  scr_dir0 = TRIM(dirname)//'/'//fchname(1:i-1)
  inpname0 = fchname(1:i-1)//'.in'
  i = RENAME(TRIM(inpname0), TRIM(inpname))

  i = INDEX(inpname, '.in', back=.true.)
  scr_dir = TRIM(dirname)//'/'//inpname(1:i-1)
  call remove_dir(TRIM(scr_dir))
  i = RENAME(TRIM(scr_dir0), TRIM(scr_dir))
 end if
end subroutine fch2qchem_wrap

subroutine bas_fch2py_wrap(fchname, dft, pyname)
 implicit none
 integer :: i, SYSTEM, RENAME
 character(len=240) :: pyname0
 character(len=240), intent(in) :: fchname
 character(len=240), intent(in), optional :: pyname
 character(len=256) :: buf
 logical, intent(in) :: dft

 buf = 'bas_fch2py '//TRIM(fchname)
 if(dft) buf = TRIM(buf)//' -dft'
 i = SYSTEM(TRIM(buf))
 if(i /= 0) call prt_call_util_error('bas_fch2py', fchname)

 if(PRESENT(pyname)) then
  i = INDEX(fchname, '.fch', back=.true.)
  pyname0 = fchname(1:i-1)//'.py'
  i = RENAME(TRIM(pyname0), TRIM(pyname))
 end if
end subroutine bas_fch2py_wrap

subroutine fch2com_wrap(fchname, inpname)
 implicit none
 integer :: i, SYSTEM, RENAME
 character(len=240) :: inpname1
 character(len=240), intent(in) :: fchname
 character(len=240), intent(in), optional :: inpname

 i = SYSTEM('fch2com '//TRIM(fchname))
 if(i /= 0) call prt_call_util_error('fch2com', fchname)

 if(PRESENT(inpname)) then
  i = INDEX(fchname, '.fch')
  inpname1 = fchname(1:i-1)//'.com'
  if(TRIM(inpname) /= TRIM(inpname1)) then
   i = RENAME(TRIM(inpname1), TRIM(inpname))
  end if
 end if
end subroutine fch2com_wrap

subroutine fch2inporb_wrap(fchname, prt_no, inpname)
 implicit none
 integer :: i, k, SYSTEM, RENAME
 character(len=240) :: old_inp, old_orb, orbname
 character(len=240), intent(in) :: fchname
 character(len=240), intent(in), optional :: inpname
 character(len=255) :: buf
 logical, intent(in) :: prt_no

 buf = 'fch2inporb '//TRIM(fchname)
 if(prt_no) buf = TRIM(buf)//' -no'
 i = SYSTEM(TRIM(buf))
 if(i /= 0) call prt_call_util_error('fch2inporb', fchname)

 if(PRESENT(inpname)) then
  k = INDEX(fchname, '.fch')
  old_inp = fchname(1:k-1)//'.input'
  if(TRIM(inpname) /= TRIM(old_inp)) then
   old_orb = fchname(1:k-1)//'.INPORB'
   call find_specified_suffix(inpname, '.input', i)
   orbname = inpname(1:i-1)//'.INPORB'
   i = RENAME(TRIM(old_inp), TRIM(inpname))
   i = RENAME(TRIM(old_orb), TRIM(orbname))
   call modify_orbname_in_molcas_inp(inpname, orbname)
  end if
 end if
end subroutine fch2inporb_wrap

subroutine orb2fch_wrap(orbname, fchname, prt_no)
 implicit none
 integer :: i, SYSTEM
 character(len=500) :: buf
 character(len=240), intent(in) :: orbname, fchname
 logical, intent(in) :: prt_no

 buf = 'orb2fch '//TRIM(orbname)//' '//TRIM(fchname)
 if(prt_no) buf = TRIM(buf)//' -no'
 i = SYSTEM(TRIM(buf))
 if(i /= 0) call prt_call_util_error('orb2fch', orbname)
end subroutine orb2fch_wrap

! modify the orbital filename in a (Open)Molcas .input file
subroutine modify_orbname_in_molcas_inp(inpname, orbfile)
 implicit none
 integer :: i, j, fid, fid1, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname, orbfile

 call find_specified_suffix(inpname, '.inp', i)
 inpname1 = inpname(1:i-1)//'.t'

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:7) == 'FILEORB') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(i /= 0) then
  close(fid)
  close(fid1,status='delete')
  write(6,'(/,A)') "ERROR in subroutine modify_orbname_in_molcas_inp: keyword '&
                   &FILEORB' not found"
  write(6,'(A)') 'in file '//TRIM(inpname)
  stop
 end if

 j = INDEX(buf, '=')
 buf = buf(1:j)//' '//TRIM(orbfile)
 write(fid1,'(A)') TRIM(buf)

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine modify_orbname_in_molcas_inp

subroutine fch2cfour_wrap(fchname)
 implicit none
 integer :: i, SYSTEM
 character(len=240), intent(in) :: fchname

 i = SYSTEM('fch2cfour '//TRIM(fchname))
 if(i /= 0) call prt_call_util_error('fch2cfour', fchname)
end subroutine fch2cfour_wrap

subroutine dat2fch_wrap(datname, fchname)
 implicit none
 integer :: i, SYSTEM
 character(len=240) :: fchname0
 character(len=240), intent(in) :: datname
 character(len=240), intent(in), optional :: fchname
 character(len=500) :: buf

 if(present(fchname)) then
  fchname0 = fchname
 else
  call find_specified_suffix(datname, '.dat', i)
  fchname0 = datname(1:i-1)//'.fch'
 end if
 buf = 'dat2fch '//TRIM(datname)//' '//TRIM(fchname0)
 write(6,'(A)') '$'//TRIM(buf)

#ifdef _WIN32
 i = SYSTEM(TRIM(buf)//' > NUL')
#else
 i = SYSTEM(TRIM(buf)//' > /dev/null')
#endif

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine dat2fch_wrap: failed to call dat2fch.'
  write(6,'(A)') 'datname='//TRIM(datname)
  write(6,'(A)') 'fchname='//TRIM(fchname0)
  stop
 end if
end subroutine dat2fch_wrap

! call `orca_2mkl` to convert .gbw -> .molden
subroutine gbw2molden(gbwname, molden)
 implicit none
 integer :: i, SYSTEM
 character(len=240) :: molden0
 character(len=240), intent(in) :: gbwname
 character(len=240), intent(in), optional :: molden
 character(len=511) :: buf

 call find_specified_suffix(gbwname, '.gbw', i)
 if(present(molden)) then
  molden0 = molden
 else
  molden0 = gbwname(1:i-1)//'.molden'
 end if
 buf = 'orca_2mkl '//TRIM(gbwname)//' '//TRIM(molden0)//' -molden -anyorbs >'

#ifdef _WIN32
 i = SYSTEM(TRIM(buf)//' NUL')
#else
 i = SYSTEM(TRIM(buf)//' /dev/null')
#endif

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine gbw2molden: failed to call orca_2mkl.'
  write(6,'(A)') 'Please check your ORCA environment variables.'
  stop
 end if
end subroutine gbw2molden

! wrapper for molden2fch
subroutine molden2fch_wrap(molden, fchname, prog, natorb)
 implicit none
 integer :: i, SYSTEM, RENAME
 character(len=7), intent(in) :: prog ! lower case
 character(len=240) :: fchname0
 character(len=240), intent(in) :: molden, fchname
 character(len=500) :: buf
 logical, intent(in) :: natorb

 call find_specified_suffix(molden, '.molden', i)
 fchname0 = molden(1:i-1)//'.fch'

 buf = 'molden2fch '//TRIM(molden)//' -'//TRIM(prog)
 if(natorb) buf = TRIM(buf)//' -no'

#ifdef _WIN32
 i = SYSTEM(TRIM(buf)//' > NUL')
#else
 i = SYSTEM(TRIM(buf)//' > /dev/null')
#endif
 if(i /= 0) call prt_call_util_error('molden2fch', molden)

 if(TRIM(fchname0) /= TRIM(fchname)) then
  i = RENAME(TRIM(fchname0), TRIM(fchname))
 end if
end subroutine molden2fch_wrap

! wrapper for subroutine gvb_exclude_XH_A
subroutine gvb_exclude_XH_A_wrap(datname, gmsname, reverted, new_inp)
 implicit none
 integer :: i, fid, SYSTEM
 character(len=500) :: buf
 character(len=18), parameter :: txt = 'gvb_exclude_XH.txt'
 character(len=240), intent(in) :: datname, gmsname
 character(len=240), intent(out) :: new_inp
 logical, intent(in) :: reverted

 buf = 'gvb_exclude_XH '//TRIM(datname)//' '//TRIM(gmsname)
 if(reverted) buf = TRIM(buf)//' r' ! reverted use of this utility
 write(6,'(/,A)') '$'//TRIM(buf)
 i = SYSTEM(TRIM(buf)//' >'//txt//" 2>&1")

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine gvb_exclude_XH_A_wrap: failed to call&
                  & utility gvb_exclude_XH.'
  write(6,'(A)') 'Did you delete it or forget to compiled it?'
  stop
 end if

 open(newunit=fid,file=txt,status='old',position='rewind')
 do i = 1, 3
  read(fid,'(A)') buf
  write(6,'(A)') TRIM(buf)
 end do ! for i
 close(fid,status='delete')

 i = INDEX(buf, ':')
 read(buf(i+1:),*) new_inp
end subroutine gvb_exclude_XH_A_wrap

subroutine prt_call_util_error(utilname, fname)
 implicit none
 character(len=*), intent(in) :: utilname
 character(len=240), intent(in) :: fname

 write(6,'(/,A)') 'ERROR in subroutine '//utilname//'_wrap: failed to call util&
                  &ity '//utilname//'.'
 write(6,'(A)') 'fname='//TRIM(fname)
 stop
end subroutine prt_call_util_error

end module util_wrapper

