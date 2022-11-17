! written by jxzou at 20210831: do rigid/relaxed scan

! currently only rigid scan is supported
subroutine do_PES_scan()
 use mr_keyword, only: rigid_scan, relaxed_scan, scan_nstep, scan_val
 implicit none
 integer :: i
 real(kind=8), external :: calc_an_int_coor
 character(len=24) :: data_string = ' '

 if(.not. (rigid_scan .or. relaxed_scan)) return
 write(6,'(//,A)') 'Enter subroutine do_PES_scan...'

 call read_scan_var_from_gjf()
 if(rigid_scan) then
  write(6,'(A)') 'Rigid scan values:'
 else
  write(6,'(A)') 'Relaxed scan values:'
 end if
 write(6,'(10F6.3)') (scan_val(i),i=1,scan_nstep)

 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine do_PES_scan at '//TRIM(data_string)
 return
end subroutine do_PES_scan

! read scan variables/coordinates from gjf
subroutine read_scan_var_from_gjf()
 use mol, only: scan_itype, scan_atoms
 use mr_keyword, only: gjfname, hf_fch, skiphf, scan_nstep, scan_val
 implicit none
 integer :: i, j, k, m, nblank0, nblank, fid
 integer, external :: detect_ncol_in_buf
 real(kind=8) :: rtmp0, stepsize
 real(kind=8), external :: calc_an_int_coor
 real(kind=8), allocatable :: coor(:,:), rtmp(:,:)
 character(len=240) :: buf
 character(len=43), parameter :: error_warn='ERROR in subroutine read_scan_var_from_gjf:'

 scan_atoms = 0     ! initialization
 buf = ' '
 nblank = 0
 if(skiphf) then
  nblank0 = 2
 else
  nblank0 = 3
 end if

 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) nblank = nblank + 1
  if(nblank == nblank0) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') error_warn//' wrong format of scan coordinate.'
  close(fid)
  stop
 end if

 read(fid,'(A)') buf
 close(fid)

 select case(buf(1:1))
 case('B')
  i = 5
 case('A')
  i = 6
 case('D')
  i = 7
 case default
  write(6,'(A)') error_warn//' invalid scan variable: '//buf(1:1)
  stop
 end select

 j = detect_ncol_in_buf(buf)
 if(j < i) then
  write(6,'(A)') error_warn//' invalid rigid scan syntax.'
  write(6,'(A)') 'buf='//TRIM(buf)
  stop
 end if

 select case(buf(1:1))
 case('B') ! bond
  scan_itype = 1
  read(buf(3:),*) scan_atoms(1:2), scan_nstep
 case('A') ! angle
  scan_itype = 2
  read(buf(3:),*) scan_atoms(1:3), scan_nstep
 case('D') ! dihedral
  scan_itype = 3
  read(buf(3:),*) scan_atoms(1:4), scan_nstep
 case default
  write(6,'(A)') error_warn//" invalid scan variable '"//buf(1:1)//"'"
  stop
 end select

 if(scan_nstep < 1) then
  write(6,'(A,I0)') error_warn//' invalid scan step=', scan_nstep
  stop
 end if

 allocate(scan_val(scan_nstep), source=0d0)
 j = LEN_TRIM(buf)

 if(buf(j:j) == '}') then ! given a set of values {}
  i = index(buf, '{')
  read(buf(i+1:j-1),fmt=*,iostat=m) (scan_val(k),k=1,scan_nstep)
  if(m /= 0) then
   write(6,'(/,A)') error_warn//' wrong scan syntax.'
   write(6,'(A)') 'buf='//TRIM(buf)
   stop
  end if

 else ! given the step length/interval
  ! calculate the current value of the target scan variable
  call read_natom_from_fch(hf_fch, k)
  allocate(coor(3,k))
  call read_coor_from_fch(hf_fch, k, coor)
  k = scan_itype + 1
  allocate(rtmp(3,k))
  forall(i = 1:k) rtmp(:,i) = coor(:,scan_atoms(i))
  deallocate(coor)
  rtmp0 = calc_an_int_coor(k, rtmp(:,1:k))
  deallocate(rtmp)
  i = index(buf(1:j), ' ', back=.true.)
  read(buf(i+1:j),*) stepsize
  forall(i = 1:scan_nstep) scan_val(i) = rtmp0 + DBLE(i)*stepsize
 end if

 call check_scan_val(scan_itype, scan_nstep, scan_val)
 return
end subroutine read_scan_var_from_gjf

! check whether the array scan_val is valid/reasonable
subroutine check_scan_val(scan_itype, n, scan_val)
 implicit none
 integer :: i, j
 integer, intent(in) :: scan_itype, n
 real(kind=8) :: rtmp
 real(kind=8), intent(in) :: scan_val(n)
 character(len=35), parameter :: error_warn='ERROR in subroutine check_scan_var:'

 do i = 1, n-1, 1
  rtmp = scan_val(i) - scan_val(i+1)

  select case(scan_itype)
  case(1) ! scan a bond
   if(rtmp < 0d0) then
    write(6,'(/,A)') error_warn//' scan values must be in descending order.'
    write(6,'(A,10F6.3)') 'scan_val=',(scan_val(j),j=1,n)
    stop
   end if
   if(rtmp<1d-3 .or. rtmp>10d0) then
    write(6,'(/,A)') error_warn//' intervals of the scanned bond are'
    write(6,'(A)') 'too large or too small. It should be >0.001 and <10.0.'
    stop
   end if
  case(2,3) ! scan an angle or a dihedral
   rtmp = DABS(rtmp)
   if(rtmp<1d0 .or. rtmp>100d0) then
    write(6,'(/,A)') error_warn//' intervals of the scanned angle/dihedral'
    write(6,'(A)') 'are too large or too small. It should be >1 and <100.0&
                     & degree.'
    stop
   end if
  case default
   write(6,'(/,A,I0)') error_warn//' invalid scan_itype=',scan_itype
   stop
  end select
 end do ! for i

 return
end subroutine check_scan_val

