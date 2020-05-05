! written by jxzou at 20180814
! This is a subroutine used to transform the GAMESS format of basis sets to the
!  NWChem format of basis sets, which can be used as input file for PySCF

! Note: 1) Only basis set data and ECP data in GAMESS format are supported!
!       2) Isotopes are not supported so far!

! updated by jxzou at 20180815: Effective Core Potential (ECP) supported
! updated by jxzou at 20181212: fix the bug in ECP
! updated by jxzou at 20200324: simplify code

program main
 implicit none
 integer :: i
 integer, parameter :: iout = 6
 character(len=240) :: fname = ' '

 i = iargc()
 if(i /= 1) then
  write(iout,'(/,1X,A)') 'ERROR in subroutine bas_gms2py: wrong command line arguments!'
  write(iout,'(/,1X,A,/)') 'Example: ./bas_gms2py a.inp (a.py generated)'
  stop
 end if

 call getarg(1,fname)
 call bas_gms2py(fname)
 stop
end program main

! transform the GAMESS format of basis sets to the NWChem format,
! which can be used as input file for PySCF
subroutine bas_gms2py(inpname)
 implicit none
 integer :: i, j, k, m, lmax, charge, mult
 integer :: nline, ncol, natom
 integer, parameter :: iout = 6
 integer, parameter :: inpid = 11, pyid = 12
 real :: rtmp
 real(kind=8), allocatable :: coor(:,:), new_coor(:,:)
 real(kind=8), allocatable :: prim_gau(:,:)
 real(kind=8), parameter :: Bohr_const = 0.52917721092d0
 character(len=1) :: bastype
 character(len=2) :: tmp_elem, new_bastype
 character(len=2), allocatable :: elem(:), new_elem(:)
 character(len=240), intent(in) :: inpname
 character(len=240) buffer, buf, pyname
 ! six types of angular momentum
 character(len=1), parameter :: am_type(0:5) = ['S','P','D','F','G','H']
 logical :: alive, bohrs, ecp

 inquire(file=TRIM(inpname), exist=alive)
 if(.not. alive) then
  write(iout,'(A)') 'ERROR in subroutine bas_gms2py: input file '//&
                    TRIM(inpname)//' does not exist!'
  stop
 end if

 buf = ' '
 buffer = ' '
 pyname = ' '
 tmp_elem = ' '
 bastype = ' '
 new_bastype = ' '
 i = INDEX(inpname, '.', back=.true.)
 pyname = inpname(1:i-1)//'.py'

 ! find whether the coordinates are in au or Bohr
 bohrs = .false.
 open(unit=inpid,file=TRIM(inpname),status='old',position='rewind')

 ! find whether there exists charge and mult
 read(inpid,'(A)',iostat=i) buf
 charge = 0
 i = index(buf,'ICHAR')
 if(i /= 0) then
  read(buf(i+7:),*) charge
 end if
 mult = 1
 i = index(buf,'MULT')
 if(i /= 0) then
  read(buf(i+5:),*) mult
 end if
 mult = mult - 1 ! mult is PySCF is na - nb

 BACKSPACE(inpid)
 do while(.true.)
  read(inpid,'(A)',iostat=i) buffer
  if(i /= 0) exit
  call upper(buffer)
  if(index(buffer,'UNITS=BOHR') /= 0) then
   bohrs = .true.
   exit
  end if
 end do
 ! found

 rewind(inpid)
 do while(.true.)
  read(inpid,'(A)',iostat=i) buffer
  if(i /= 0) exit
  call upper(buffer)
  if(buffer(2:6)=='$DATA') exit
 end do
 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine bas_gms2py: No $DATA section found&
                   & in file '//TRIM(inpname)//'.'
  close(inpid)
  stop
 end if

 ! skip 2 lines: the Title Card line and the Point Group line
 read(inpid,'(A)') buffer
 read(inpid,'(A)') buffer

 natom = 0
 allocate(new_coor(3,1), source=0.0d0)
 allocate(new_elem(1))
 new_elem = ' '

 ! read the elements and coordinates
 do while(.true.)
  read(inpid,'(A)') buffer
  ! 'buffer' contains the element, relative atomic mass and coordinates
  buf = buffer
  call upper(buf)
  if(buf(2:5) == '$END') exit

  ! deal with the coordinates
  natom = natom + 1
  if(natom > 1) then
   allocate(new_coor(3,natom), source=0.0d0)
   allocate(new_elem(natom))
   new_elem = ' '
   new_coor(:,1:natom-1) = coor
   new_elem(1:natom-1) = elem
   deallocate(coor, elem)
  end if
  read(buffer,*) new_elem(natom), rtmp, new_coor(1:3,natom)
  coor = new_coor ! update coor
  elem = new_elem ! update elem
  deallocate(new_coor, new_elem)

  do while(.true.)
   read(inpid,'(A)') buffer
   if(LEN_TRIM(buffer) == 0) exit
   read(buffer,*) bastype, nline
   do i = 1, nline, 1
    read(inpid,'(A)') buffer
   end do ! for i
  end do ! for while

 end do ! for while
 ! read done

 if(natom == 0) then
  write(iout,'(A)') 'ERROR in subroutine bas_gms2py: zero atom found!'
  close(inpid)
  stop
 end if

 ! print the coordinates (optional) into .nwc file
 open(unit=pyid,file=TRIM(pyname),status='replace')
 write(pyid,'(A)') 'from pyscf import gto, scf'
 if(bohrs) coor = coor*Bohr_const
 write(pyid,'(/,A)') 'mol = gto.M()'
 write(pyid,'(A,I4,A)') '#', natom, ' atom(s)'
 write(pyid,'(A)') "mol.atom = '''"
 do i = 1, natom, 1
  write(pyid,'(A2,2X,3(1X,F17.8))') elem(i), coor(1:3,i)
 end do
 write(pyid,'(A)') "'''"
 deallocate(coor)
 ! print coordinates done

 rewind(inpid)
 do while(.true.)
  read(inpid,'(A)') buffer
  call upper(buffer)
  if(buffer(2:6)=='$DATA') exit
 end do

 ! skip 2 lines: the Title Card line and the Point Group line
 read(inpid,'(A)') buffer
 read(inpid,'(A)') buffer

 ! read and print basis sets
 write(pyid,'(/,A)') 'mol.basis = {'
 do k = 1, natom, 1
  read(inpid,'(A)') buffer
  ! 'buffer' contains the element, relative atomic mass and coordinates
  buf = buffer
  call upper(buf)
  if(buf(2:5)=='$END') exit

  ! deal with the element, i.e., atomic symbol
  read(buffer,*) tmp_elem
  tmp_elem = ADJUSTL(tmp_elem)

  ! deal with the basis data
  write(pyid,'(A)') "'"//TRIM(tmp_elem)//"': gto.basis.parse('''"
  do while(.true.)
   read(inpid,'(A)') buffer
   if(LEN_TRIM(buffer) == 0) exit

   read(buffer,*) bastype, nline
   if(bastype == 'L') then
    ncol = 3
    new_bastype = 'SP'
   else
    ncol = 2
    new_bastype(1:2) = bastype//' '
   end if
   allocate(prim_gau(ncol,nline), source=0.0d0)
   write(pyid,'(A2,4X,A2)') tmp_elem, new_bastype
   do i = 1, nline, 1
    read(inpid,*) j, prim_gau(1:ncol,i)
    if(nline==1 .and. prim_gau(2,1)<1.0d-7) prim_gau(2,1) = 1.0d0
    write(pyid,'(3(1X,E18.9))') prim_gau(:,i)
   end do ! for i
   deallocate(prim_gau)
  end do ! for while

  write(pyid,'(A)') "'''),"
 end do ! k

 BACKSPACE(pyid)
 write(pyid,'(A)') "''')}"

 ! check if any ECP data
 ! Note: it is required in inp file that ECP is defined for every atom.
 !  For no ECP atom, it must be assigned as 'ECP-NONE'
 rewind(inpid)
 ecp = .false.
 do while(.true.)
  read(inpid,'(A)',iostat=i) buffer
  if(i /= 0) exit
  call upper(buffer)
  if(buffer(2:5) == '$ECP') then
   ecp = .true.
   exit
  end if
 end do

 if(ecp) then
  write(pyid,'(/,A)') 'mol.ecp = {'

  do m = 1, natom, 1
   read(inpid,'(A)') buffer
   buf = buffer
   call upper(buf)
   if(buf(2:5) == '$END') exit
   if(index(buf,'NONE') /= 0) cycle

   write(pyid,'(A)') "'"//TRIM(elem(m))//"': gto.basis.parse_ecp('''"
   i = index(buf,'GEN')
   read(buf(i+3:),*) k, lmax
   write(pyid,'(A,I3)') TRIM(elem(m))//' nelec ', k

   allocate(prim_gau(2,1), source=0.0d0)
   do k = 1, lmax+1, 1
    read(inpid,'(A)') buffer
    if(k == 1) then
     new_bastype = 'ul'
    else
     new_bastype = am_type(k-2)
    end if
    write(pyid,'(A)') TRIM(elem(m))//' '//TRIM(new_bastype)
    read(buffer,*) nline
    do i = 1, nline, 1
     read(inpid,*) prim_gau(1,1), j, prim_gau(2,1)
     write(pyid,'(I1,1X,F15.8,3X,F15.8)') j, prim_gau(2,1), prim_gau(1,1)
    end do ! for i
   end do ! for k

   deallocate(prim_gau)
   write(pyid,'(A)') "'''),"
  end do ! for m

  BACKSPACE(pyid)
  write(pyid,'(A)') "''')}"
 end if ! for ecp

 close(inpid)
 deallocate(elem)

 write(pyid,'(/,A)') '# Remember to check the charge and spin'
 write(pyid,'(A,I0)') 'mol.charge = ', charge
 write(pyid,'(A,I0)') 'mol.spin = ', mult
 write(pyid,'(A)') 'mol.verbose = 4'
 write(pyid,'(A)') 'mol.cart = True'
 write(pyid,'(A,/)') 'mol.build()'
 close(pyid)
 return
end subroutine bas_gms2py

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

