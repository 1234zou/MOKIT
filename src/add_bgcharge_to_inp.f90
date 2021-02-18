! written by jxzou at 20200729: add background point charges into input files of
!  various programs (.py, .inp, .input, etc)
! updated by jxzou at 20210216: add support of PSI4

! In the following comments, we use 'n' for nuclear, 'e' for background point charges,
!  and n-e for nuclear-charge interactions, e-e for self energy of point charges

! Note: the e-e is
!       (1) taken into account in energy of Gaussian
!       (2) taken into account in energy of ORCA when using Q in .inp file
!       (3) not taken into account in energy of PySCF, OpenMolcas

! When using $QUANPO in GAMESS,
! For GAMESS CASCI : e-e and n-e not taken into account in CASCI energy
! For GAMESS GVB   : e-e and n-e not taken into account in 1st cycle energy of GVB
! For GAMESS GVB   : e-e and n-e taken into account in >1st cycle energy of GVB
! For GAMESS CASSCF: e-e and n-e taken into account in each cycle energy

! These findings were tested on GAMESS 2017.

! The e-e and n-e energy are considered in energies in output of AutoMR.
program main
 implicit none
 integer :: i
 integer, parameter :: iout = 6
 character(len=38), parameter :: error_warn='ERROR in program add_bgcharge_to_inp: '
 character(len=240) :: chgname, inpname

 i = iargc()
 if(i /= 2) then
  write(iout,'(/,A)') error_warn//'wrong command line arguments!'
  write(iout,'(/,A)') 'Format: add_bgcharge_to_inp chgname inpname'
  write(iout,'(A)')   'Example 1 (PySCF) : add_bgcharge_to_inp a.chg a.py'
  write(iout,'(A)')   'Example 2 (GAMESS): add_bgcharge_to_inp a.chg a.inp'
  write(iout,'(A)')   'Example 3 (Molcas): add_bgcharge_to_inp a.chg a.input'
  write(iout,'(A)')   'Example 4 (ORCA)  : add_bgcharge_to_inp a.chg a.inp'
  write(iout,'(A)')   'Example 5 (Molpro): add_bgcharge_to_inp a.chg a.com'
  write(iout,'(A)')   'Example 6 (BDF)   : add_bgcharge_to_inp a.chg a.inp'
  write(iout,'(A,/)') 'Example 7 (PSI4)  : add_bgcharge_to_inp a.chg a.inp'
  stop
 end if

 call getarg(1, chgname)
 call require_file_exist(chgname)
 call getarg(2, inpname)
 call require_file_exist(inpname)

 call add_bgcharge_to_inp(chgname, inpname)
 stop
end program main

! add background charges into input files of various programs (.py, .inp, .input, etc)
subroutine add_bgcharge_to_inp(chgname, inpname)
 implicit none
 integer :: i, j, n, fid
 integer, parameter :: iout = 6
 real(kind=8), allocatable :: charge(:,:)
 character(len=41),parameter::error_warn='ERROR in subroutine add_bgcharge_to_inp: '
 character(len=240) :: buf
 character(len=240), intent(in) :: chgname, inpname

 open(newunit=fid,file=TRIM(chgname),status='old',position='rewind')
 read(fid,*,iostat=i) n
 if(i /= 0) then
  write(iout,'(A)') error_warn//'failed to read number of charges!'
  stop
 end if

 if(n < 1) then
  write(iout,'(A)') error_warn//'number of charges less than 1!'
  stop
 end if

 allocate(charge(4,n), source=0d0)
 do i = 1, n, 1
  read(fid,*,iostat=j) charge(1:4,i)
  if(j /= 0) exit
 end do ! for i
 close(fid)

 if(j /= 0) then
  write(iout,'(A)') error_warn//' missing charges.'
  stop
 end if

 i = index(inpname, '.', back=.true.)

 select case(TRIM(inpname(i+1:)))
 case('gjf')
  call add_bgcharge_to_gjf(inpname, n, charge)
 case('py')    ! PySCF .py file
  call add_bgcharge_to_py(inpname, n, charge)
 case('input') ! OpenMolcas .input file
  call add_bgcharge_to_molcas_input(inpname, n, charge)
 case('inp')
  open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
  read(fid,'(A)') buf
  close(fid)

  call lower(buf)
  if(buf(2:8) == '$contrl') then  ! GAMESS .inp file
   call add_bgcharge_to_gms_inp(inpname, n, charge)
  else if(buf(1:3) == '***') then ! Molpro .inp file
   call add_bgcharge_to_molpro_inp(inpname, n, charge)
  else if(buf(1:8) == '$compass') then ! BDF .inp file
   call add_bgcharge_to_bdf_inp(inpname, n, charge)
  else if(buf(1:4)=='%pal' .or. buf(1:8)=='%maxcore') then ! ORCA .inp file
   call add_bgcharge_to_orca_inp(inpname, n, charge)
  else ! PSI4 .inp file
   call add_bgcharge_to_psi4_inp(inpname, n, charge)
  end if

 case('com')
  open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
  read(fid,'(A)') buf
  close(fid)

  if(buf(1:3) == '***') then   ! Molpro .inp file
   call add_bgcharge_to_molpro_inp(inpname, n, charge) 
  else                         ! Gaussian .com file
   call add_bgcharge_to_gjf(inpname, n, charge)
  end if
 case default
  write(iout,'(A)') error_warn//'filetype not supported!'
  stop
 end select

 return
end subroutine add_bgcharge_to_inp

! add background charges into a Gaussian .gjf file
subroutine add_bgcharge_to_gjf(gjfname, n, charge)
 implicit none
 integer :: i, nblank, fid1, fid2, RENAME
 integer, intent(in) :: n
 real(kind=8), intent(in) :: charge(4,n)
 character(len=240) :: gjfname1, buf
 character(len=240), intent(in) :: gjfname

 gjfname1 = TRIM(gjfname)//'.tmp'
 open(newunit=fid1,file=TRIM(gjfname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(gjfname1),status='replace')

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:1) == '#') exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 write(fid2,'(A,/)') TRIM(buf)//' charge'

 if(index(buf,'geom=allcheck') /= 0) then
  ! do nothing
 else if(index(buf,'geom=check') /= 0) then
  do i = 1, 5, 1
   read(fid2,'(A)') buf
   write(fid2,'(A)') TRIM(buf)
  end do ! for i
 else   ! skip title, charge and Cartesian coordinates
  read(fid1,'(A)') buf
  nblank = 0
  do while(.true.)
   read(fid1,'(A)') buf
   write(fid2,'(A)') TRIM(buf)
   if(LEN_TRIM(buf) == 0) nblank = nblank + 1
   if(nblank == 2) exit
  end do ! for while
 end if

 do i = 1, n, 1
  write(fid2,'(4F17.8)') charge(:,i)
 end do ! for i

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:1) == '#') then
   write(fid2,'(A)') TRIM(buf)//' charge=check'
  else
   write(fid2,'(A)') TRIM(buf)
  end if
 end do ! for while

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(TRIM(gjfname1), TRIM(gjfname))
 return
end subroutine add_bgcharge_to_gjf

! add background charges into a PySCF .py input file
subroutine add_bgcharge_to_py(pyname, n, charge)
 implicit none
 integer :: i, fid1, fid2, RENAME
 integer, intent(in) :: n
 real(kind=8), intent(in) :: charge(4,n)
 character(len=240) :: buf, pyname1
 character(len=240), intent(in) :: pyname

 pyname1 = TRIM(pyname)//'.tmp'
 open(newunit=fid1,file=TRIM(pyname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(pyname1),status='replace')
 write(fid2,'(A)') 'from pyscf import qmmm'

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(1:9) == 'mf.kernel') exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 write(fid2,'(A)') 'bgcoord = ['
 do i = 1, n, 1
  write(fid2,'(3(A1,F17.8),A2)') '[',charge(1,i),',',charge(2,i),',',charge(3,i),'],'
 end do ! for i
 write(fid2,'(A)') ']'

 write(fid2,'(A)') 'bgcharge = ['
 write(fid2,'(8(E15.8,A1))') (charge(4,i),',',i=1,n)
 write(fid2,'(A)') ']'

 write(fid2,'(A)') 'mf = qmmm.mm_charge(mf, bgcoord, bgcharge)'
 write(fid2,'(A)') 'mf.kernel()'

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(TRIM(pyname1), TRIM(pyname))
 return
end subroutine add_bgcharge_to_py

! add background charges into a GAMESS .inp input file
subroutine add_bgcharge_to_gms_inp(inpname, n, charge)
 implicit none
 integer :: i, fid1, fid2, natom, RENAME
 integer, intent(in) :: n
 integer, allocatable :: nuc(:)
 real(kind=8), intent(in) :: charge(4,n)
 real(kind=8), allocatable :: coor(:,:)
 character(len=2), allocatable :: elem(:)
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 call read_natom_from_gms_inp(inpname, natom)
 allocate(coor(3,natom), nuc(natom), elem(natom))
 call read_elem_nuc_coor_from_gms_inp(inpname, natom, elem, nuc, coor)

 inpname1 = TRIM(inpname)//'.tmp'
 open(newunit=fid1,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(2:6) == '$DATA') exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 write(fid2,'(A)') ' $QUANPO NFFTYP=0 $END'
 write(fid2,'(A)') ' $FFDATA'
 write(fid2,'(A)') '  COORDINATES   NUC   X   Y   Z'
 do i = 1, natom, 1
  write(fid2,'(3X,A2,1X,I0,A2,3(1X,F17.8))') elem(i), nuc(i),'.0', coor(1:3,i)
 end do ! for i
 do i = 1, n, 1
  write(fid2,'(3X,A,3(1X,F17.8))') 'Q   0.0',charge(1:3,i)
 end do ! for i
 write(fid2,'(A)') '  STOP'

 write(fid2,'(A)') '  PARAMETERS   MASS   Q   POL   RMIN/2  EPSILON  RMIN/2  EPSILON'
 do i = 1, natom, 1
  write(fid2,'(3X,A2,9X,A)') elem(i), '0.0     0.0      0.0 0.0 0.0 0.0 0.0'
 end do ! for i
 do i = 1, n, 1
  write(fid2,'(3X,A,2X,F11.6,A)') 'Q          0.0',charge(4,i),' 0.0 0.0 0.0 0.0 0.0'
 end do ! for i
 write(fid2,'(A)') '  STOP'
 write(fid2,'(A)') ' $END'

 deallocate(coor, elem, nuc)
 BACKSPACE(fid1)

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
 return
end subroutine add_bgcharge_to_gms_inp

! add background charges into a ORCA .inp input file
subroutine add_bgcharge_to_orca_inp(inpname, n, charge)
 implicit none
 integer :: i, iend, fid1, fid2, RENAME
 integer, intent(in) :: n
 integer, parameter :: iout = 6
 real(kind=8), intent(in) :: charge(4,n)
 character(len=240) :: buf, mklname, mklname1
 character(len=240), intent(in) :: inpname
 logical :: begin

 i = index(inpname, '.inp', back=.true.)
 mklname = inpname(1:i-1)//'.mkl'
 mklname1 = TRIM(mklname)//'.tmp'
 open(newunit=fid1,file=TRIM(inpname),status='old',position='append')

 iend = 0
 do while(.true.)
  BACKSPACE(fid1)
  BACKSPACE(fid1)
  read(fid1,'(A)') buf
  if(index(buf,'end') /= 0) iend = iend + 1
  if(iend == 3) exit
 end do ! for while

 do i = 1, n, 1
  write(fid1,'(A,I0,A1,1X,F11.6,3(1X,F17.8))') '  Q(',i,')',charge(4,i),charge(1:3,i)
 end do ! for i
 write(fid1,'(A)') ' end'
 write(fid1,'(A)') 'end'
 close(fid1)

 open(newunit=fid1,file=TRIM(mklname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(mklname1),status='replace')

 begin = .false.
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:6) == '$COORD') begin = .true.
  if(begin .and. buf(1:4)=='$END') exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine add_bgcharge_to_orca_inp: section&
                   & '$COORD' in file '"//TRIM(mklname)//"' is incomplete."
  close(fid1)
  close(fid2,status='delete')
  stop
 end if

 do i = 1, n, 1
  write(fid2,'(A,3(1X,F17.8))') ' 113', charge(1:3,i)
 end do ! for i
 write(fid2,'(A)') '$END'

 begin = .false.
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:8) == '$CHARGES') begin = .true.
  if(begin .and. buf(1:4)=='$END') exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine add_bgcharge_to_orca_inp: section&
                  & '$CHARGES' in file '"//TRIM(mklname)//"' is incomplete."
  close(fid1)
  close(fid2,status='delete')
  stop
 end if

 do i = 1, n, 1
  write(fid2,'(F11.6)') charge(4,i)
 end do ! for i
 write(fid2,'(A)') '$END'

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(TRIM(mklname1), TRIM(mklname))
 return
end subroutine add_bgcharge_to_orca_inp

! add background charges into a (Open)Molcas .input file
subroutine add_bgcharge_to_molcas_input(input, n, charge)
 implicit none
 integer :: i, fid1, fid2, RENAME
 integer, intent(in) :: n
 integer, parameter :: iout = 6
 real(kind=8), intent(in) :: charge(4,n)
 real(kind=8), parameter :: Bohr_const = 0.52917721092d0
 character(len=240) :: buf, input1
 character(len=240), intent(in) :: input

 input1 = TRIM(input)//'.tmp'
 open(newunit=fid1,file=TRIM(input),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(input1),status='replace')

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:7) == "&SEWARD") exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine add_bgcharge_to_molcas_inp: file '"//&
   TRIM(input)//"' is incomplete."
  close(fid1)
  close(fid2,status='delete')
  stop
 end if

 write(fid2,'(A,/,I0)') 'XField', n
 do i = 1, n, 1
  write(fid2,'(4F17.8,A)') charge(1:3,i)/Bohr_const, charge(4,i),' 0.0 0.0 0.0'
 end do ! for i
 write(fid2,'(A,/)') 'End of Input'

 BACKSPACE(fid1)
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(TRIM(input1), TRIM(input))

 return
end subroutine add_bgcharge_to_molcas_input

! add background charges into a Molpro .input file
subroutine add_bgcharge_to_molpro_inp(inpname, n, charge)
 implicit none
 integer :: i, fid, fid1, RENAME
 integer, intent(in) :: n
 integer, parameter :: iout = 6
 real(kind=8), intent(in) :: charge(4,n)
 character(len=240) :: buf, inpname1, chgname
 character(len=240), intent(in) :: inpname

 inpname1 = TRIM(inpname)//'.tmp'
 i = index(inpname, '.', back=.true.)
 chgname = inpname(1:i-1)//'.chg1'
 call lower(chgname)

 open(newunit=fid,file=TRIM(chgname),status='replace')
 write(fid,'(A)') 'Molpro background charge file generated by AutoMR of MOKIT'
 write(fid,'(I0)') n
 do i = 1, n, 1
  write(fid,'(4(1X,F18.8),A4)') charge(1:4,i), '   0'
 end do ! for i
 close(fid)

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit

  call lower(buf(1:5))
  if(buf(1:5) == 'basis') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine add_bgcharge_to_molpro_inp:'
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 write(fid1,'(A)') 'Lattice,infile='//TRIM(chgname)
 BACKSPACE(fid)

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
 return
end subroutine add_bgcharge_to_molpro_inp

! generate the .extcharge file for a given BDF .inp file, and
! add 'Extcharge' keyword into the given BDF .inp file
subroutine add_bgcharge_to_bdf_inp(inpname, n, charge)
 implicit none
 integer :: i, fid, fid1, RENAME
 integer, intent(in) :: n
 integer, parameter :: iout = 6
 real(kind=8), intent(in) :: charge(4,n)
 character(len=4) :: str
 character(len=240) :: buf, chgname, inpname1
 character(len=240), intent(in) :: inpname

 str = ' '
 i = index(inpname, '.', back=.true.)
 chgname = inpname(1:i-1)//'.extcharge'
 inpname1 = inpname(1:i-1)//'.tmp'

 open(newunit=fid,file=TRIM(chgname),status='replace')
 write(fid,'(A)') 'BDF background charge file generated by AutoMR of MOKIT'
 write(fid,'(I0)') n
 do i = 1, n, 1
  write(fid,'(A2,4(1X,F18.8))') 'H ', charge(4,i), charge(1:3,i)
 end do ! for i
 close(fid)

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  str = buf(1:4)
  call lower(str)
  if(str == '$end') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(i /= 0) then
  close(fid1,status='delete')
  close(fid)
  write(iout,'(A)') 'ERROR in subroutine add_bgcharge_to_bdf_inp: '
  stop
 end if

 write(fid1,'(A)') 'Extcharge'
 write(fid1,'(A)') ' point'
 write(fid1,'(A)') '$END'

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
 return
end subroutine add_bgcharge_to_bdf_inp

subroutine add_bgcharge_to_psi4_inp(inpname, n, charge)
 implicit none
 integer :: i, fid, fid1, RENAME
 integer, intent(in) :: n
 real(kind=8), intent(in) :: charge(4,n)
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 inpname1 = TRIM(inpname)//'.t'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:3) == 'set') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 write(fid1,'(A)') 'Chrgfield = QMMM()'
 do i = 1, n, 1
  write(fid1,'(A,F11.6,A,3(F16.8,A))') 'Chrgfield.extern.addCharge(',charge(4,i),&
   ',',charge(1,i),',',charge(2,i),',',charge(3,i),')'
 end do ! for i

 write(fid1,'(A)') "psi4.set_global_option_python('EXTERN', Chrgfield.extern)"
 write(fid1,'(/,A)') 'set {'

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
 return
end subroutine add_bgcharge_to_psi4_inp

