! written by jxzou at 20201209: move some subroutines from bas_gms2molcas.f90 to this file

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

  prim_gau(2)%stype = 'P'
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
   if(DABS(coeff(nline0,1) - exp1) < zero) then
    allocate(prim_gau(k)%coeff(nline0,ncol+1), source=0.0d0)
    prim_gau(k)%coeff(1:nline0,1:ncol) = coeff
    deallocate(coeff)
    ncol = ncol + 1
    prim_gau(k)%ncol = ncol
    read(fid1,*)   ! skip one line
    prim_gau(k)%coeff(nline0,1) = exp1
    prim_gau(k)%coeff(nline0,ncol) = 1.0d0
   else
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
 use pg, only: prim_gau, natom, ram, highest, all_ecp, elem, ecp_exist
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

! update the number of times each atom occurred
subroutine calc_ntimes(natom, elem, ntimes)
 implicit none
 integer :: i, j
 integer, intent(in) :: natom
 integer, intent(out) :: ntimes(natom)
 character(len=2) :: tmp
 character(len=2), intent(in) :: elem(natom)

 ntimes = 1
 do i = 1, natom, 1
  tmp = elem(i)

  do j = 1, i-1, 1
   if(tmp == elem(j)) ntimes(i) = ntimes(i) + 1
  end do ! for j
 end do ! for i
 return
end subroutine calc_ntimes

