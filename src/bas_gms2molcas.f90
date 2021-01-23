! written by jxzou at 20191030: transform basis sets in GAMESS format to (Open)Molcas format
! updated by jxzou at 20200325: add variable rtmp
! updated by jxzou at 20200326: move lower angular momentum (occurred 2nd time) forward
!   This special case happens at, e.g. def2TZVP for Co. The MOs moved forward have been
!   considered in fch2inporb. Here the basis sets should be adjusted correspondingly.
! updated by jxzou at 20200402: ECP/PP supported
! updated by jxzou at 20201110: dealing with the extra primitive duplicate, e.g.
!   cc-pVDZ for Co. The original way is to both extend column and row. Now the
!   trick is to only extend column, thus will not generate the primitive duplicate.
!   The primitive duplicate is not allowed in relativistic calculations.
! updated by jxzou at 20201209: move some subroutines to read_gms_inp.f90

! Note: Currently isotopes are not tested.
program main
 use pg, only: iout
 implicit none
 integer :: i
 character(len=4) :: str
 character(len=240) :: fname
 ! fname: input file contains basis sets and Cartesian coordinates in GAMESS format
 logical :: alive, spherical

 i = iargc()
 if(i<1 .or. i>2) then
  write(iout,'(/,A,/)') 'Example1: bas_gms2molcas a.inp (generating an a.input file)'
  write(iout,'(A,/)') "Example2: bas_gms2molcas a.inp -sph (without 'Cartesian all')"
  stop
 end if

 str = ' '
 fname = ' '
 spherical = .false.
 call getarg(1,fname)

 inquire(file=TRIM(fname),exist=alive)
 if(.not. alive) then
  write(iout,'(A)') 'ERROR in subroutine bas_gms2molcas: file does not exist!'
  write(iout,'(A)') 'Filename='//TRIM(fname)
  stop
 end if

 if(i == 2) then
  call getarg(2,str)

  if(str == '-sph') then
   spherical = .true.
  else
   write(iout,'(A)') 'ERROR in subroutine bas_gms2molcas: wrong command line arguments!'
   write(iout,'(A)') 'str='//str
   stop
  end if

 end if

 call bas_gms2molcas(fname, spherical)
 stop
end program main

! Transform the basis sets in GAMESS format to those in (Open)Molcas format
subroutine bas_gms2molcas(fort7, spherical)
 use pg, only: iout, natom, ram, ntimes, coor, elem, prim_gau, all_ecp, ecp_exist
 implicit none
 integer :: i, j, k, m, n, nline, rc, rel
 integer :: fid1, fid2
 integer :: charge, mult
 real(kind=8) :: rtmp
 real(kind=8) :: exp1   ! the first exponent
 real(kind=8), parameter :: Bohr_const = 0.52917721092d0
 real(kind=8), allocatable :: coor2(:,:)
 character(len=240), intent(in) :: fort7
 character(len=240) :: buf1, buf2, input
 ! input is the (Open)Molcas .input file
 character(len=1) :: stype0, stype
 logical, intent(in) :: spherical
 logical :: bohrs, uhf, X2C

 ! initialization
 buf1 = ' '
 buf2 = ' '
 input = ' '

 k = INDEX(fort7, '.', .true.)
 input = fort7(1:k-1)//'.input'

 open(newunit=fid2,file=TRIM(input),status='replace')
 write(fid2,'(A)') "&GATEWAY"
! write(fid2,'(A)') 'Expert'

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
  all_ecp(:)%ecp = .false.
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
      write(iout,'(A)') "ERROR in subroutine bas_gms2molcas: ECP/PP does not&
                       & starts with '-ul potential'."
      write(iout,'(A)') 'You should check the format of ECP/PP data in file '//TRIM(fort7)
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
  write(fid2,'(A)') 'Basis set'
  if(ecp_exist) then
   if(all_ecp(natom)%ecp) then
    write(fid2,'(A)') TRIM(elem(natom))//'.ECP....      / inline'
   else
    write(fid2,'(A)') TRIM(elem(natom))//'.....      / inline'
   end if
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
    ! an angular momentum occurs more than once
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

 call check_X2C_in_gms_inp(fort7, X2C)
 if(X2C) write(fid2,'(A)') 'RX2C'
 write(fid2,'(/,A)') "&SEWARD"

 call check_DKH_in_gms_inp(fort7, rel)
 select case(rel)
 case(-2) ! nothing
 case(-1) ! RESC
  write(iout,'(A)') 'ERROR in subroutine bas_gms2molcas: RESC keywords detected.'
  write(iout,'(A)') 'But RESC is not supported in (Open)Molcas.'
  stop
 case(0,1,2,4)  ! DKH0/1/2/4
  if(.not. X2C) write(fid2,'(A,I2.2,A)') 'Relativistic = R',rel,'O'
  !if(.not. X2C) write(fid2,'(A,I2.2,A)') 'R',rel,'O'
 case default
  write(iout,'(A)') 'ERROR in subroutine bas_gms2molcas: rel out of range!'
  write(iout,'(A,I0)') 'rel=', rel
  stop
 end select

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

