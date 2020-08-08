! written by jxzou at 20191030
! This is a subroutine used to transform the basis sets in GAMESS format to those in Molcas/OpenMolcas format

! updated by jxzou at 20200325: add variable rtmp
! updated by jxzou at 20200326: move lower angular momentum (occurred 2nd time) forward
!  This special case happens at, e.g. def2TZVP for Co. The MOs moved forward have been
!  considered in fch2inporb. Here the basis sets should be adjusted correspondingly.
! updated by jxzou at 20200402: ECP/PP supported

! Note: Currently isotopes are not tested.

module pg
 implicit none
 integer :: natom     ! the number of atoms
 integer :: highest   ! highest angular momentum
 integer, parameter :: iout = 6
 integer, allocatable :: ram(:), ntimes(:) ! ram: relative atomic mass
 real(kind=8), parameter :: zero = 1.0d-7
 real(kind=8), allocatable :: coor(:,:)    ! Cartesian coordinates
 character(len=2), allocatable :: elem(:)  ! elements

 ! 'L' will be divided into two parts: 'S' and 'P'
 type primitive_gaussian
  character(len=1) :: stype = ' ' ! 'S', 'P', 'D', 'F', 'G', 'H', 'I'
  integer :: nline = 0
  integer :: ncol  = 0
  real(kind=8), allocatable :: coeff(:,:)
 end type primitive_gaussian

 ! 7 for 'S', 'P', 'D', 'F', 'G', 'H', 'I'
 !        1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7
 type(primitive_gaussian) :: prim_gau(7)

 ! --- below are arrays for ECP/PP ---
 type ecp_potential
  integer :: n = 0   ! size for col1, 2, and 3
  real(kind=8), allocatable :: col1(:)
  integer     , allocatable :: col2(:)
  real(kind=8), allocatable :: col3(:)
 end type ecp_potential

 type ecp4atom
  logical :: ecp = .false. ! whether this atom has ECP/PP
  integer :: core_e  = 0   ! number of core electrons
  integer :: highest = 0   ! highest angular momentum
  type(ecp_potential), allocatable :: potential(:) ! size highest+1
 end type ecp4atom

 type(ecp4atom), allocatable :: all_ecp(:) ! size natom
 logical :: ecp_exist = .false.
end module pg

program main
 use pg, only: iout
 implicit none
 integer :: i
 character(len=4) :: str
 character(len=240) :: fname
 ! fname: input file contains basis sets and Cartesian coordinates in GAMESS format
 logical :: spherical

 i = iargc()
 if(i<1 .or. i>2) then
  write(iout,'(/,A,/)') 'Example1: ./bas_gms2molcas a.inp (generating an a.input file)'
  write(iout,'(A,/)') "Example2: ./bas_gms2molcas a.inp -sph (do not write 'Cartesian all')"
  stop
 end if

 str = ' '
 fname = ' '
 spherical = .false.
 call getarg(1,fname)

 if(i == 2) then
  call getarg(2,str)

  if(str == '-sph') then
   spherical = .true.
  else
   write(iout,'(A)') 'ERROR in subroutine bas_gms2molcas: wrong command line parameters.'
   write(iout,'(A)') 'str='//str
   stop
  end if

 end if

 call bas_gms2molcas(fname, spherical)
 stop
end program main

! Transform the basis sets in GAMESS format to those in Molcas/OpenMolcas format
subroutine bas_gms2molcas(fort7, spherical)
 use pg, only: iout, natom, ram, ntimes, coor, elem, prim_gau, all_ecp, ecp_exist
 implicit none
 integer :: i, j, k, m, n, nline, rc
 integer :: fid1, fid2
 integer :: charge, mult
 real(kind=8) :: rtmp
 real(kind=8) :: exp1   ! the first exponent
 real(kind=8), parameter :: Bohr_const = 0.52917721092d0
 real(kind=8), allocatable :: coor2(:,:)
 character(len=240), intent(in) :: fort7
 character(len=240) :: buf1, buf2, input ! input is the Molcas .input file
 character(len=1) :: stype0, stype
 logical, intent(in) :: spherical
 logical :: bohrs, uhf

 ! initialization
 buf1 = ' '
 buf2 = ' '
 input = ' '

 k = INDEX(fort7, '.', .true.)
 input = fort7(1:k-1)//'.input'

 open(newunit=fid2,file=TRIM(input),status='replace')
 write(fid2,'(A)') "&GATEWAY"
 write(fid2,'(A)') 'Expert'

 open(newunit=fid1,file=TRIM(fort7),status='old',position='rewind')

 uhf = .false.
 read(fid1,'(A)') buf1
 call upper(buf1)
 i = index(buf1,'UHF')
 if(i > 0) uhf = .true.

 charge = 0; mult = 1
 i = index(buf1,'ICHAR')
 if(i > 0) read(buf1(i+7:),*) charge
 i = index(buf1,'MULT')
 if(i > 0) read(buf1(i+5:),*) mult
 BACKSPACE(fid1)

 ! find in the first 3 lines whether the coordinates are in Angstrom or Bohr
 bohrs = .false.
 do i = 1, 3
  read(fid1,'(A)') buf1
  call upper(buf1)
  if(INDEX(buf1,'UNITS=BOHR') /= 0) then
   bohrs = .true.
   exit
  end if
 end do
 ! Angstrom/Bohr determined

 ! rewind and check if any ECP data exists
 rewind(fid1)
 ecp_exist = .false.
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf1
  if(i /= 0) exit
  call upper(buf1)
  if(buf1(2:5) == '$ECP') then
   ecp_exist = .true.
   exit
  end if
 end do

 natom = 0
 if(ecp_exist) then ! if exists, find natom
  do while(.true.)
   read(fid1,'(A)') buf1
   call upper(buf1)
   if(buf1(2:5) == '$END') exit
   if(index(buf1,'ECP') /= 0) natom = natom + 1
  end do ! for while

  allocate(all_ecp(natom))
  rewind(fid1)
  do while(.true.)
   read(fid1,'(A)') buf1
   call upper(buf1)
   if(buf1(2:5) == '$ECP') exit
  end do

  do i = 1, natom, 1
   read(fid1,'(A)') buf1
   call upper(buf1)
   if(index(buf1,'NONE') /= 0) cycle

   all_ecp(i)%ecp = .true.
   k = index(buf1,'GEN')
   read(buf1(k+3:),*) all_ecp(i)%core_e, m
   all_ecp(i)%highest = m
   allocate(all_ecp(i)%potential(m+1))

   do j = 1, m+1, 1
    read(fid1,'(A)') buf1
    if(j == 1) then
     k = index(buf1,'ul')
     if(k == 0) then
      write(iout,'(A)') "ERROR in subroutine bas_gms2molcas: ECP/PP not starts with '-ul potential'."
      write(iout,'(A)') 'You should check the format of ECP/PP data in file '//TRIM(fort7)//'.'
      stop
     end if
    end if

    read(buf1,*) n
    all_ecp(i)%potential(j)%n = n
    allocate(all_ecp(i)%potential(j)%col1(n), source=0.0d0)
    allocate(all_ecp(i)%potential(j)%col2(n), source=0)
    allocate(all_ecp(i)%potential(j)%col3(n), source=0.0d0)
    do k = 1, n, 1
     read(fid1,*) all_ecp(i)%potential(j)%col1(k), all_ecp(i)%potential(j)%col2(k), &
     all_ecp(i)%potential(j)%col3(k)
    end do ! for k
   end do ! for j
  end do ! for i

 end if
 ! read ECP/PP done

 ! rewind and find the $DATA section
 rewind(fid1)
 do while(.true.)
  read(fid1,'(A)',iostat=rc) buf1
  if(rc /= 0) exit
  call upper(buf1)
  if(buf1(2:6) == '$DATA') exit
 end do
 if(rc /= 0) then
  write(iout,'(A)') 'ERROR in subroutine bas_gms2molcas: No $DATA section found&
                   & in file '//TRIM(fort7)//'.'
  close(fid1)
  stop
 end if

 ! skip 2 lines: the Title line and the Point Group line
 read(fid1,'(A)') buf1
 read(fid1,'(A)') buf1

 ! initialization: clear all primitive gaussians
 call clear_prim_gau()

 ! read the element, relative atomic mass (ram) and coordinates
 natom = 0
 do while(.true.) ! while 1
  read(fid1,'(A)',iostat=rc) buf1
  if(rc /= 0) exit
  ! 'buf1' contains the element, ram and coordinates
  buf2 = buf1
  call upper(buf2)
  if(buf2(2:5) == '$END') exit

  ! deal with the coordinates
  natom = natom + 1
  if(natom > 1) then ! enlarge the arrays
   allocate(coor2(3,natom), source=0.0d0)
   coor2(:,1:natom-1) = coor
   coor = coor2   ! auto-reallocation
   deallocate(coor2)
   elem = [elem, ' ']
   ram = [ram, 0]
   ntimes = [ntimes, 1]
  else ! natom = 1
   allocate(elem(1))
   elem(1) = ' '
   allocate(ram(1), source=0)
   allocate(ntimes(1), source=1)
   allocate(coor(3,1), source=0.0d0)
  end if
  read(buf1,*) elem(natom), rtmp, coor(1:3,natom)
  ram(natom) = INT(rtmp)
  if(bohrs) coor(1:3,natom) = coor(1:3,natom)*Bohr_const
  write(fid2,'(/,A)') 'Basis set'
  if(ecp_exist) then
   if(all_ecp(natom)%ecp) write(fid2,'(A)') TRIM(elem(natom))//'.ECP....      / inline'
  else
   write(fid2,'(A)') TRIM(elem(natom))//'.....      / inline'
  end if

  ! deal with primitive gaussians
  stype0 = ' '
  do while(.true.) ! while 2
   read(fid1,'(A)') buf1
   if(LEN_TRIM(buf1) == 0) exit

   read(buf1,*) stype, nline
   if(stype0 == ' ') then
    !------------------------------------------------------
    ! the following 5 lines are added to determine whether
    ! a angular momentum occurs more than once
    call stype2itype(stype, k)
    if( allocated(prim_gau(k)%coeff) ) then
     stype0 = stype
     BACKSPACE(fid1)
    else
    !------------------------------------------------------
     call read_prim_gau1(stype, nline, fid1)
     stype0 = stype
    end if
   else ! stype0 /= ' '
    if(stype == stype0) then
     exp1 = 0.0d0
     read(fid1,*) i, exp1
     BACKSPACE(fid1)
     call read_prim_gau2(stype, nline, fid1, exp1)
    else ! stype /= stype0
     stype0 = ' ' ! reset stype0
     BACKSPACE(fid1)
    end if
   end if
  end do ! for while 2

  ! print basis sets and ECP/PP (if any) of this atom in Molcas format
  call prt_prim_gau(fid2)

  ! print elem and coor
  call update_ntimes()
  write(fid2,'(A,I0,3(2X,F15.8),3X,A)') TRIM(elem(natom)), ntimes(natom), &
                                        coor(1:3,natom),'Angstrom'

  if(.not. spherical) write(fid2,'(A)') 'Cartesian all'
  write(fid2,'(A)') 'End of basis set'

  ! clear all primitive gaussians for next cycle
  call clear_prim_gau()
 end do ! for while 1

 if(allocated(all_ecp)) deallocate(all_ecp)
 close(fid1)
 if(natom==0 .or. rc/=0) then
  if(natom == 0) write(iout,'(A)') 'ERROR in subroutine bas_gms2molcas: zero atom found!'
  if(rc /= 0) write(iout,'(A)') "ERROR in subroutine bas_gms2molcas: it seems the '$DATA'&
                               & has no corresponding '$END'. EOF."
  close(fid2,status='delete')
  stop
 end if

 write(fid2,'(/,A)') "&SEWARD"
 write(fid2,'(/,A)') "&SCF"
 if(uhf) write(fid2,'(A3)') 'UHF'
 write(fid2,'(A,I0)') 'Charge = ', charge
 write(fid2,'(A,I0)') 'Spin = ', mult
 k = INDEX(fort7, '.', .true.)
 input = fort7(1:k-1)//'.INPORB'
 write(fid2,'(A,/)') 'FILEORB = '//TRIM(input)
 close(fid2)
 return
end subroutine bas_gms2molcas

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

! convert a (character) stype to (integer) itype
subroutine stype2itype(stype, itype)
 use pg, only: iout
 implicit none
 integer, intent(out) :: itype
 character(len=1), intent(in) :: stype

 itype = 0
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
 case('L')
 case default
  write(iout,'(A)') 'ERROR in subroutine stype2itype: stype out of range.'
  write(iout,'(A)') 'stype= '//TRIM(stype)
  stop
 end select

 return
end subroutine stype2itype

! deallocate the allocatable arrays in array prim_gau
subroutine clear_prim_gau()
 use pg, only: prim_gau
 implicit none
 integer :: i

 do i = 1, 7, 1
  prim_gau(i)%stype = ' '
  prim_gau(i)%nline = 0
  prim_gau(i)%ncol = 0
  if(allocated(prim_gau(i)%coeff)) deallocate(prim_gau(i)%coeff)
 end do
 return
end subroutine clear_prim_gau

! read this type of primitive gaussians in the 1st time
! i.e., the 1st 'S', the 1st 'L', etc.
subroutine read_prim_gau1(stype, nline, fid1)
 use pg, only: zero, prim_gau
 implicit none
 integer :: i, k, itmp, ncol, nline0
 integer, intent(in) :: nline, fid1
 real(kind=8) :: rtmp
 real(kind=8), allocatable :: coeff(:,:)
 character(len=1), intent(in) :: stype

 call stype2itype(stype, k)

 if(k /= 0) then
  prim_gau(k)%stype = stype
  prim_gau(k)%nline = nline
  prim_gau(k)%ncol = 2
  allocate(prim_gau(k)%coeff(nline,2), source=0.0d0)
  do i = 1, nline, 1
   read(fid1,*) itmp, prim_gau(k)%coeff(i,1), rtmp
   if(nline==1 .and. DABS(rtmp)<zero) rtmp = 1.0d0
   prim_gau(k)%coeff(i,2) = rtmp
  end do ! for i

 else ! k == 0, this is 'L' for some Pople-type basis sets
  ncol = prim_gau(1)%ncol
  nline0 = prim_gau(1)%nline
  coeff = prim_gau(1)%coeff ! auto-allocation
  deallocate(prim_gau(1)%coeff)
  allocate(prim_gau(1)%coeff(nline0+nline,ncol+1), source=0.0d0) ! 'S'
  prim_gau(1)%coeff(1:nline0,1:ncol) = coeff
  deallocate(coeff)
  ncol = ncol + 1
  prim_gau(1)%ncol = ncol

  prim_gau(2)%stype = stype
  prim_gau(2)%nline = nline
  prim_gau(2)%ncol = 2
  allocate(prim_gau(2)%coeff(nline,2), source=0.0d0) ! 'P'
  do i = 1, nline, 1
   read(fid1,*) itmp, prim_gau(1)%coeff(nline0+i,1), prim_gau(1)%coeff(nline0+i,ncol), &
                prim_gau(2)%coeff(i,2)
  end do ! for i
  forall(i = 1:nline) prim_gau(2)%coeff(i,1) = prim_gau(1)%coeff(nline0+i,1)
 end if

 return
end subroutine read_prim_gau1

! read this type of primitive gaussians after the 1st time
! i.e., the 2nd 'S', the 3rd 'S', the 2nd 'L', etc.
subroutine read_prim_gau2(stype, nline, fid1, exp1)
 use pg, only: zero, iout, prim_gau
 implicit none
 integer :: i, k, itmp, ncol, nline0
 integer :: ncol2, nline2, dim1
 integer, intent(in) :: nline, fid1
 real(kind=8) :: rtmp
 real(kind=8), intent(in) :: exp1
 real(kind=8), allocatable :: coeff(:,:), coeff2(:,:)
 character(len=1), intent(in) :: stype

 call stype2itype(stype, k)

 if(k /= 0) then
  nline0 = size(prim_gau(k)%coeff,1)
  ncol = size(prim_gau(k)%coeff,2)
  coeff = prim_gau(k)%coeff   ! auto-allocation
  dim1 = size(coeff,dim=1)
  deallocate(prim_gau(k)%coeff)

  if(DABS(coeff(1,1) - exp1) < zero) then
   ! For some strange basis set(e.g. LANL2DZ for Cu), nline>size(coeff,dim=1)
   ! Normally nline = size(coeff,dim=1)
   allocate(prim_gau(k)%coeff(MAX(nline,dim1),ncol+1), source=0.0d0)
   prim_gau(k)%coeff(1:dim1,1:ncol) = coeff
   deallocate(coeff)
   ncol = ncol + 1
   prim_gau(k)%ncol = ncol
   do i = 1, MIN(nline,dim1), 1
    read(fid1,*) itmp, rtmp, prim_gau(k)%coeff(i,ncol)
   end do ! for i
   if(nline > dim1) then
    do i = dim1+1, nline, 1
     read(fid1,*) itmp, prim_gau(k)%coeff(i,1), prim_gau(k)%coeff(i,ncol)
    end do ! for i
   end if
  else ! > zero, enlarge the array prim_gau(k)%coeff
   allocate(prim_gau(k)%coeff(nline0+nline,ncol+1), source=0.0d0)
   prim_gau(k)%coeff(1:nline0,1:ncol) = coeff
   deallocate(coeff)
   ncol = ncol + 1
   prim_gau(k)%ncol = ncol
   do i = 1, nline, 1
    read(fid1,*) itmp, prim_gau(k)%coeff(nline0+i,1), rtmp
    if(nline==1 .and. DABS(rtmp)<zero) rtmp = 1.0d0
    prim_gau(k)%coeff(nline0+i,ncol) = rtmp
   end do ! for i
   nline0 = nline0 + nline
   prim_gau(k)%nline = nline0
  end if

 else ! k == 0, this is 'L' for some Pople-type basis sets
  if(DABS(prim_gau(2)%coeff(1,1) - exp1) < zero) then
   write(iout,'(A)') "ERROR in subroutine read_prim_gau2: two sections of 'L'&
                    & have identical gaussian exponents. Impossible."
   write(iout,'(A)') 'Please check the input file.'
   stop
  else ! > zero, enlarge the array prim_gau(1)%coeff and prim_gau(2)%coeff
   nline0 = size(prim_gau(1)%coeff,1)
   ncol = size(prim_gau(1)%coeff,2)
   coeff = prim_gau(1)%coeff   ! auto-allocation
   deallocate(prim_gau(1)%coeff)
   allocate(prim_gau(1)%coeff(nline0+nline,ncol+1), source=0.0d0)
   prim_gau(1)%coeff(1:nline0,1:ncol) = coeff
   deallocate(coeff)
   ncol = ncol + 1
   prim_gau(1)%ncol = ncol

   nline2 = size(prim_gau(2)%coeff,1)
   ncol2 = size(prim_gau(2)%coeff,2)
   coeff2 = prim_gau(2)%coeff   ! auto-allocation
   deallocate(prim_gau(2)%coeff)
   allocate(prim_gau(2)%coeff(nline2+nline,ncol2+1), source=0.0d0)
   prim_gau(2)%coeff(1:nline2,1:ncol2) = coeff2
   deallocate(coeff2)
   ncol2 = ncol2 + 1
   prim_gau(2)%ncol = ncol2

   do i = 1, nline, 1
    read(fid1,*) itmp, rtmp, prim_gau(1)%coeff(nline0+i,ncol), prim_gau(2)%coeff(nline2+i,ncol2)
    prim_gau(1)%coeff(nline0+i,1) = rtmp
    prim_gau(2)%coeff(nline2+i,1) = rtmp
   end do ! for i
   nline0 = nline0 + nline
   nline2 = nline2 + nline
   prim_gau(1)%nline = nline0
   prim_gau(2)%nline = nline2
  end if
 end if

 return
end subroutine read_prim_gau2

! determine the highest angular momentum quantum number
subroutine get_highest_am()
 use pg, only: highest, prim_gau
 implicit none

 if(allocated(prim_gau(7)%coeff)) then
  highest = 6
 else if(allocated(prim_gau(6)%coeff)) then
  highest = 5
 else if(allocated(prim_gau(5)%coeff)) then
  highest = 4
 else if(allocated(prim_gau(4)%coeff)) then
  highest = 3
 else if(allocated(prim_gau(3)%coeff)) then
  highest = 2
 else if(allocated(prim_gau(2)%coeff)) then
  highest = 1
 else
  highest = 0
 end if
 return
end subroutine get_highest_am

! print primitive gaussians
subroutine prt_prim_gau(fid)
 use pg, only: zero, prim_gau, natom, ram, highest, all_ecp, elem, ecp_exist
 implicit none
 integer :: i, j, k, m, n, nline, ncol
 integer, intent(in) :: fid
 integer, allocatable :: list(:)
 character(len=1), parameter :: am(0:5) = ['S','P','D','F','G','H']

 call get_highest_am()
 if(ecp_exist) then
  write(fid,'(5X,I0,A1,3X,I1)') ram(natom)-all_ecp(natom)%core_e, '.', highest
 else
  write(fid,'(5X,I0,A1,3X,I1)') ram(natom), '.', highest
 end if

 do i = 1, 7, 1
  if(.not. allocated(prim_gau(i)%coeff)) cycle
  write(fid,'(A)') '* '//prim_gau(i)%stype//'-type functions'
  nline = prim_gau(i)%nline
  ncol = prim_gau(i)%ncol
  write(fid,'(2(1X,I4))') nline, ncol-1
  do j = 1, nline, 1
   write(fid,'(3X,E16.9)') prim_gau(i)%coeff(j,1)
  end do ! for j
  do j = 1, nline, 1
   write(fid,'(10(E16.9,3X))') (prim_gau(i)%coeff(j,k), k=2,ncol)
  end do ! for j
 end do ! for i

 if(.not. ecp_exist) return

 if(all_ecp(natom)%ecp) then
  m = all_ecp(natom)%highest
  write(fid,'(A2,1X,A2,1X,I0,1X,I0)') 'PP', elem(natom), all_ecp(natom)%core_e, m
  allocate(list(m+1))
  list(1) = m
  forall(i=2:m+1) list(i) = i-2
  do i = 1, m+1, 1
   n = all_ecp(natom)%potential(i)%n
   if(i == 1) then
    write(fid,'(I0,A)') n,';!'//am(list(1))//' POTENTIAL'
   else ! i > 1
    write(fid,'(I0,A)') n,';!'//am(list(i))//'-'//am(list(1))//' POTENTIAL'
   end if

   do j = 1, n, 1
    write(fid,'(I0,2(1X,F16.8))') all_ecp(natom)%potential(i)%col2(j), all_ecp(natom)%potential(i)%col3(j), &
                                  all_ecp(natom)%potential(i)%col1(j)
   end do ! for j
  end do ! for i
  write(fid,'(A1,/,A,/,A)') '*', 'Spectral', 'End of Spectral'
 end if

 return
end subroutine prt_prim_gau

! update the number of times each atom occurred
subroutine update_ntimes()
 use pg, only: natom, elem, ntimes
 implicit none
 integer :: i
 character(len=2) :: ctemp

 ctemp = ' '
 ctemp = elem(natom)

 do i = 1, natom-1, 1
  if(ctemp == elem(i)) ntimes = ntimes + 1
 end do
 return
end subroutine update_ntimes

