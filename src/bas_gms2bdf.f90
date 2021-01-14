! written by jxzou at 20210104: copied from file bas_gms2molcas.f90 and modified

! Note: Currently isotopes are not tested.
program main
 use pg, only: iout
 implicit none
 integer :: i
 character(len=4) :: str
 character(len=240) :: fname
 ! fname: input file contains basis sets and Cartesian coordinates in GAMESS format
 logical :: alive, uhf

 i = iargc()
 if(i<1 .or. i>2) then
  write(iout,'(/,A)') ' ERROR in subroutine bas_gms2bdf: wrong command line arguments!'
  write(iout,'(A)') ' Example 1: bas_gms2bdf a.inp (generating an a_bdf.inp file)'
  write(iout,'(A,/)') ' Example 2: bas_gms2bdf a.inp -uhf'
  stop
 end if

 fname = ' '
 call getarg(1,fname)
 fname = ADJUSTL(fname)

 inquire(file=TRIM(fname),exist=alive)
 if(.not. alive) then
  write(iout,'(A)') 'ERROR in subroutine bas_gms2bdf: file does not exist!'
  write(iout,'(A)') 'Filename='//TRIM(fname)
  stop
 end if

 uhf = .false.
 if(i == 2) then
  str = ' '
  call getarg(2,str)
  if(str == '-uhf') then
   uhf = .true.
  else
   write(iout,'(A)') 'ERROR in subroutine bas_gms2bdf: wrong command line arguments!'
   write(iout,'(A)') "The 2nd argument can only be '-uhf'."
   stop
  end if
 end if

 call bas_gms2bdf(fname, uhf)
 stop
end program main

! Transform the basis sets in GAMESS format to those in BDF format
subroutine bas_gms2bdf(fort7, uhf)
 use pg, only: iout, natom, ram, ntimes, coor, elem, prim_gau, all_ecp, ecp_exist
 implicit none
 integer :: i, j, k, m, n, nline, rc, rel
 integer :: fid1, fid2
 integer :: charge, mult
 real(kind=8) :: rtmp
 real(kind=8) :: exp1   ! the first exponent
 character(len=7) :: str
 character(len=240), intent(in) :: fort7
 character(len=240) :: buf, input, basfile
 ! input is the BDF input file
 ! if you do not like the suffix .bdf, you can change it into .inp
 character(len=1) :: stype0, stype
 logical :: X2C
 logical, intent(in) :: uhf

 ! initialization
 buf = ' '
 input = ' '
 basfile = ' '

 k = index(fort7, '.', back=.true.)
 input = fort7(1:k-1)//'_bdf.inp'
 call read_natom_from_gms_inp(fort7, natom)

 open(newunit=fid1,file=TRIM(fort7),status='old',position='rewind')
 read(fid1,'(A)') buf
 call upper(buf)

 charge = 0; mult = 1
 i = index(buf,'ICHAR')
 if(i > 0) read(buf(i+7:),*) charge
 i = index(buf,'MULT')
 if(i > 0) read(buf(i+5:),*) mult
 BACKSPACE(fid1)

 ! check if any ECP data exists
 ecp_exist = .false.
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  call upper(buf(3:5))
  if(buf(2:5) == '$ECP') then
   ecp_exist = .true.
   exit
  end if
 end do ! for while

 if(ecp_exist) then ! if exists, read ECP
  basfile = 'ECP.'//fort7(1:k-1)//'.BAS'
  allocate(all_ecp(natom))

  do i = 1, natom, 1
   read(fid1,'(A)') buf
   call upper(buf)
   if(index(buf,'NONE') /= 0) cycle

   all_ecp(i)%ecp = .true.
   k = index(buf,'GEN')
   read(buf(k+3:),*) all_ecp(i)%core_e, m
   all_ecp(i)%highest = m
   allocate(all_ecp(i)%potential(m+1))

   do j = 1, m+1, 1
    read(fid1,'(A)') buf
    if(j == 1) then
     k = index(buf,'ul')
     if(k == 0) then
      write(iout,'(A)') "ERROR in subroutine bas_gms2bdf: ECP/PP not starts with '-ul potential'."
      write(iout,'(A)') 'You should check the format of ECP/PP data in file '//TRIM(fort7)//'.'
      stop
     end if
    end if

    read(buf,*) n
    all_ecp(i)%potential(j)%n = n
    allocate(all_ecp(i)%potential(j)%col1(n), source=0d0)
    allocate(all_ecp(i)%potential(j)%col2(n), source=0)
    allocate(all_ecp(i)%potential(j)%col3(n), source=0d0)
    do k = 1, n, 1
     read(fid1,*) all_ecp(i)%potential(j)%col1(k), all_ecp(i)%potential(j)%col2(k), &
     all_ecp(i)%potential(j)%col3(k)
    end do ! for k
   end do ! for j
  end do ! for i
 else
  basfile = fort7(1:k-1)//'.BAS'
 end if
 ! read ECP/PP done

 close(fid1)
 call upper(basfile)
 open(newunit=fid2,file=TRIM(input),status='replace')
 write(fid2,'(A)') '$COMPASS'
 write(fid2,'(A)') 'Title'
 write(fid2,'(A)') 'auto-generated file by the bas_gms2bdf utility of MOKIT'
 write(fid2,'(A)') 'Geometry'

 allocate(elem(natom), ram(natom), coor(3,natom), ntimes(natom))
 call read_elem_nuc_coor_from_gms_inp(fort7, natom, elem, ram, coor)
 call calc_ntimes(natom, elem, ntimes)
 call check_X2C_in_gms_inp(fort7, X2C)
 call check_DKH_in_gms_inp(fort7, rel)

 do i = 1, natom, 1
  str = ' '
  write(str,'(A,I0)') TRIM(elem(i)), ntimes(i)
  write(fid2,'(A7,1X,3F18.8)') str, coor(1:3,i)
 end do ! for i

 deallocate(coor)
 write(fid2,'(A)') 'End Geometry'
 write(fid2,'(A)') 'Basis'
 write(fid2,'(A)') TRIM(basfile)
 write(fid2,'(A)') 'Nosymm'
 write(fid2,'(A)') '$END'
 write(fid2,'(/,A)') '$XUANYUAN'

 if(rel>-1 .or. X2C) then
  write(fid2,'(A)') 'Scalar'

  if(rel>-1 .and. (.not.X2C)) then
   write(iout,'(A)') 'Warning in subroutine bas_gms2bdf: the BDF program does&
                    & not support DKH Hamiltonian.'
   write(iout,'(A)') "Spin-free X2C keyword 'Scalar' are written in $XUANYUAN instead."
  end if
 end if

 write(fid2,'(A)') '$END'
 write(fid2,'(/,A)') '$SCF'
 write(fid2,'(A)') 'CheckLin'
 write(fid2,'(A)') 'TolLin'
 write(fid2,'(A)') '1.D-6'
 if(uhf) then
  write(fid2,'(A3)') 'UHF'
 else
  if(mult /= 1) then
   write(fid2,'(A4)') 'ROHF'
  else
   write(fid2,'(A3)') 'RHF'
  end if
 end if
 write(fid2,'(A,/,I0)') 'Charge', charge
 write(fid2,'(A,/,I0)') 'Spin', mult
 write(fid2,'(A)') 'Guess'
 write(fid2,'(A)') 'read'
 write(fid2,'(A)') '$END'
 close(fid2)

 ! reopen the .inp/.dat file and find the $DATA section
 open(newunit=fid1,file=TRIM(fort7),status='old',position='rewind')
 do while(.true.)
  read(fid1,'(A)',iostat=rc) buf
  if(rc /= 0) exit
  call upper(buf)
  if(buf(2:6) == '$DATA') exit
 end do ! for while
 if(rc /= 0) then
  write(iout,'(A)') 'ERROR in subroutine bas_gms2bdf: No $DATA section found&
                   & in file '//TRIM(fort7)//'.'
  close(fid1)
  stop
 end if

 ! skip 2 lines: the Title line and the Point Group line
 read(fid1,'(A)') buf
 read(fid1,'(A)') buf

 ! initialization: clear all primitive gaussians
 call clear_prim_gau()
 open(newunit=fid2,file=TRIM(basfile),status='replace')

 ! read the element, relative atomic mass (ram) and coordinates
 do i = 1, natom, 1
  read(fid1,'(A)',iostat=rc) buf
  ! 'buf' contains the element, ram and coordinates
  if(rc /= 0) exit
  if(buf(2:2) == '$') call upper(buf(3:5))
  if(buf(2:5) == '$END') exit

  ! deal with primitive gaussians
  stype0 = ' '
  do while(.true.)
   read(fid1,'(A)') buf
   if(LEN_TRIM(buf) == 0) exit

   read(buf,*) stype, nline
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
     read(fid1,*) j, exp1
     BACKSPACE(fid1)
     call read_prim_gau2(stype, nline, fid1, exp1)
    else ! stype /= stype0
     stype0 = ' ' ! reset stype0
     BACKSPACE(fid1)
    end if
   end if
  end do ! for while

  ! print basis sets and ECP/PP (if any) of this atom in BDF format
  call prt_prim_gau_bdf(i, fid2)

  ! clear all primitive gaussians for next cycle
  call clear_prim_gau()
 end do ! for i

 close(fid1)
 close(fid2)
 deallocate(ram, ntimes, elem)
 if(allocated(all_ecp)) deallocate(all_ecp)
 return
end subroutine bas_gms2bdf

! print primitive gaussians
subroutine prt_prim_gau_bdf(iatom, fid)
 use pg, only: prim_gau, ram, highest, all_ecp, elem, ntimes, ecp_exist
 implicit none
 integer :: i, j, k, m, n, nline, ncol
 integer, intent(in) :: iatom, fid
 integer, allocatable :: list(:)
 character(len=1), parameter :: am(0:5) = ['s','p','d','f','g','h']

 call get_highest_am()
 write(fid,'(A)') '****'
 write(fid,'(A,I0,3X,I0,3X,I1)') TRIM(elem(iatom)),ntimes(iatom),ram(iatom),highest

 do i = 1, 7, 1
  if(.not. allocated(prim_gau(i)%coeff)) cycle
  write(fid,'(A)',advance='no') prim_gau(i)%stype
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

 if(all_ecp(iatom)%ecp) then
  write(fid,'(A)') 'ECP'
  m = all_ecp(iatom)%highest
  write(fid,'(A,2(I0,1X),I0)') TRIM(elem(iatom)),ntimes(iatom),all_ecp(iatom)%core_e,m
  allocate(list(m+1))
  list(1) = m
  forall(i=2:m+1) list(i) = i-2
  do i = 1, m+1, 1
   n = all_ecp(iatom)%potential(i)%n
   if(i == 1) then
    write(fid,'(A,I0)') am(list(i))//' potential ', n
   else ! i > 1
    write(fid,'(A,I0)') am(list(i))//' potential ', n
   end if

   do j = 1, n, 1
    write(fid,'(I0,2(1X,F16.8))') all_ecp(iatom)%potential(i)%col2(j), all_ecp(iatom)%potential(i)%col3(j), &
                                  all_ecp(iatom)%potential(i)%col1(j)
   end do ! for j
  end do ! for i
 end if

 return
end subroutine prt_prim_gau_bdf

