! written by jxzou at 20210404: generate BDF .BAS, .inp, .scforb/.inporb files
!  from Gaussian .fch(k) file
! updated by jxzou at 20210114: enlarge coeff(nbf,nif) to coeff(nbf,nbf) for
!  linear dependence usage
! updated by jxzou at 20210407: remove '-uhf' and add automatic determination

program main
 use util_wrapper, only: formchk, fch2inp_wrap
 implicit none
 integer :: i, k, system
 character(len=3) :: str
 character(len=240) :: fchname, inpname
 logical :: prt_no

 i = iargc()
 if(i<1 .or. i>2) then
  write(6,'(/,A)') ' ERROR in subroutine fch2bdf: wrong command line arguments!'
  write(6,'(A)')   ' Example 1 (R(O)HF, UHF): fch2bdf a.fch     (-> a_bdf.inp a.scforb)'
  write(6,'(A,/)') ' Example 2 (CAS NO)     : fch2bdf a.fch -no (-> a_bdf.inp a.inporb)'
  stop
 end if

 fchname = ' '
 inpname = ' '
 call getarg(1, fchname)
 call require_file_exist(fchname)

 ! if .chk file provided, convert into .fch file automatically
 k = LEN_TRIM(fchname)
 if(fchname(k-3:k) == '.chk') then
  call formchk(fchname)
  fchname = fchname(1:k-3)//'fch'
 end if

 if(i == 2) then
  call getarg(2, str)
  if(str /= '-no') then
   write(6,'(A)') 'ERROR in subroutine fch2bdf: wrong command line arguments!'
   write(6,'(A)') "The 2nd argument can only be '-no'."
   stop
  else ! str = '-no'
   prt_no = .true.
  end if
 else ! i = 1
  prt_no = .false.
 end if

 call fch2bdf(fchname, prt_no)
 call fch2inp_wrap(fchname, .false., 0, 0)

 i = index(fchname, '.fch', back=.true.)
 inpname = fchname(1:i-1)//'.inp'
 i = system('bas_gms2bdf '//TRIM(inpname))

 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine fch2bdf: failed to call utility bas_gms2bdf.'
  write(6,'(A)') 'The file '//TRIM(fchname)//' may be incomplete.'
  stop
 end if

 call delete_file(inpname)
end program main

! read the MOs in .fch(k) file and adjust its p,d,f,g, etc. functions order
!  of Gaussian to that of BDF
subroutine fch2bdf(fchname, prt_no)
 use fch_content, only: Bohr_const, check_uhf_in_fch
 implicit none
 integer :: i, j, k, m, length, natom, orbid
 integer :: na, nb, nif, nbf, nbf0, nbf1, mult
 integer :: n3pmark, n5dmark, n7fmark, n9gmark, n11hmark
 integer, allocatable :: shell_type(:), shell2atom_map(:)
 ! mark the index where d, f, g, h functions begin
 integer, allocatable :: p_mark(:), d_mark(:), f_mark(:), g_mark(:), h_mark(:)
 integer, allocatable :: nuc(:), ntimes(:)
 real(kind=8), allocatable :: coeff(:,:), occ_num(:), coor(:,:)
 character(len=2), allocatable :: elem(:)
 character(len=240) :: orbfile
 character(len=240), intent(in) :: fchname
 logical, intent(in) :: prt_no
 logical :: uhf

 orbfile = ' '
 call check_nobasistransform_in_fch(fchname)
 call check_nosymm_in_fch(fchname)
 call read_ncontr_from_fch(fchname, k)

 allocate(shell_type(2*k), source=0)
 allocate(shell2atom_map(2*k), source=0)
 call read_shltyp_and_shl2atm_from_fch(fchname, k, shell_type, shell2atom_map)

 if( ANY(shell_type>1) ) then ! whether Cartesian/spherical harmonic
  deallocate(shell_type, shell2atom_map)
  write(6,'(A)') 'ERROR in subroutine fch2bdf: Cartesian-type basis functions&
                 & are not supported in BDF.'
  write(6,'(A)') 'You can use spherical harmonic basis functions.'
  stop
 end if

 call read_na_and_nb_from_fch(fchname, na, nb)
 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 nbf0 = nbf   ! save a copy of the original nbf

 ! read MO Coefficients
 ! The array size of occ_num is actually 2*nif or nif, and coeff is (nbf,nif)
 ! or (nbf,2*nif). But BDF requires them to use nbf. If linear dependence
 ! occurs, the extra elements are zero and must be printed into .scforb/.inporb.
 call check_uhf_in_fch(fchname, uhf) ! determine whether UHF
 if(uhf) then
  allocate(occ_num(2*nbf), coeff(nbf,2*nbf))
  call read_eigenvalues_from_fch(fchname,nif,'a',occ_num(1:nif))
  call read_eigenvalues_from_fch(fchname,nif,'b',occ_num(nif+1:2*nif))
  call read_mo_from_fch(fchname, nbf, nif, 'a', coeff(:,1:nif))
  call read_mo_from_fch(fchname, nbf, nif, 'b', coeff(:,nif+1:2*nif))
  nif = 2*nif   ! double the size
 else
  allocate(occ_num(nbf), coeff(nbf,nbf))
  call read_eigenvalues_from_fch(fchname, nif, 'a', occ_num)
  call read_mo_from_fch(fchname, nbf, nif, 'a', coeff(1:nbf,1:nif))
 end if

! first we adjust the basis functions in each MO according to the Shell to atom map
 ! 1) split the 'L' into 'S' and 'P', this is to ensure that D comes after L functions
 call split_L_func(k, shell_type, shell2atom_map, length)

 ! 2) sort the shell_type and shell2atom_map by ascending order
 ! MOs will be adjusted accordingly
 call sort_shell_and_mo(length, shell_type, shell2atom_map, nbf, nif, coeff)
! adjust done

! then we adjust the basis functions in each MO according to the type of basis functions
 k = length  ! update k
 n3pmark = 0
 n5dmark = 0
 n7fmark = 0
 n9gmark = 0
 n11hmark = 0
 allocate(p_mark(k), d_mark(k), f_mark(k), g_mark(k), h_mark(k))
 p_mark = 0
 d_mark = 0
 f_mark = 0
 g_mark = 0
 h_mark = 0
 nbf = 0
 do i = 1, k, 1
  select case(shell_type(i))
  case( 0)   ! S
   nbf = nbf + 1
  case( 1)   ! 3P
   n3pmark = n3pmark + 1
   p_mark(n3pmark) = nbf + 1
   nbf = nbf + 3
  case(-1)   ! SP or L
   n3pmark = n3pmark + 1
   p_mark(n3pmark) = nbf + 2
   nbf = nbf + 4
  case(-2)   ! 5D
   n5dmark = n5dmark + 1
   d_mark(n5dmark) = nbf + 1
   nbf = nbf + 5
  case(-3)   ! 7F
   n7fmark = n7fmark + 1
   f_mark(n7fmark) = nbf + 1
   nbf = nbf + 7
  case(-4)   ! 9G
   n9gmark = n9gmark + 1
   g_mark(n9gmark) = nbf + 1
   nbf = nbf + 9
  case(-5)   ! 11H
   n11hmark = n11hmark + 1
   h_mark(n11hmark) = nbf + 1
   nbf = nbf + 11
  case default
   write(6,'(A)') 'ERROR in subroutine fch2bdf: shell_type(i) out of range.'
   write(6,'(2(A,I0))') 'i=', i, ', k=', k
   stop
  end select
 end do ! for i

 ! adjust the order of d, f, etc. functions
 do i = 1, n3pmark, 1
  call fch2bdf_permute_3p(nif,coeff(p_mark(i):p_mark(i)+2,:))
 end do
 do i = 1, n5dmark, 1
  call fch2bdf_permute_5d(nif,coeff(d_mark(i):d_mark(i)+4,:))
 end do
 do i = 1, n7fmark, 1
  call fch2bdf_permute_7f(nif,coeff(f_mark(i):f_mark(i)+6,:))
 end do
 do i = 1, n9gmark, 1
  call fch2bdf_permute_9g(nif,coeff(g_mark(i):g_mark(i)+8,:))
 end do
 do i = 1, n11hmark, 1
  call fch2bdf_permute_11h(nif,coeff(h_mark(i):h_mark(i)+10,:))
 end do
! adjustment finished

 deallocate(p_mark, d_mark, f_mark, g_mark, h_mark)

! move the 2nd, 3rd, ... Zeta basis functions forward
 i = 0
 nbf = 0
 do while(i < k)
  i = i + 1
  j = shell2atom_map(i)
  m = shell_type(i)
  nbf1 = nbf
  select case(m)
  case( 0)   ! S
   nbf = nbf + 1
  case( 1)   ! 3P
   nbf = nbf + 3
  case(-1)   ! SP or L
   nbf = nbf + 4
  case(-2)   ! 5D
   nbf = nbf + 5
  case(-3)   ! 7F
   nbf = nbf + 7
  case(-4)   ! 9G
   nbf = nbf + 9
  case(-5)   ! 11H
   nbf = nbf + 11
  end select
  if(m == 0) cycle

  length = 1
  do while(i < k)
   i = i + 1
   if(shell_type(i) /= m) exit
   if(shell2atom_map(i) /= j) exit
   length = length + 1
  end do
  if(i < k) i = i - 1
  if(length > 1) then
   call zeta_mv_forwd(nbf1, m, length, nbf0, nif, coeff)
   nbf = nbf1 + length*(nbf-nbf1)
  end if
 end do
 deallocate(shell_type, shell2atom_map)
! move done

! print MOs into BDF .scforb or .inporb
 i = index(fchname, '.fch', back=.true.)
 if(prt_no) then
  orbfile = fchname(1:i-1)//'_bdf.inporb'
 else
  orbfile = fchname(1:i-1)//'_bdf.scforb'
 end if

 open(newunit=orbid,file=TRIM(orbfile),status='replace')
 write(orbid,'(A)') 'TITLE - transferred MOs'

 call read_charge_and_mult_from_fch(fchname, i, mult)

 if((mult==1 .and. (.not.uhf)) .or. prt_no) then
  write(orbid,'(A)') '$MOCOEF  1'  ! RHF/CAS
 else
  write(orbid,'(A)') '$MOCOEF  2'  ! ROHF, UHF
 end if

 if(uhf) nif = nif/2
 write(orbid,'(A,I9,A)') 'SYM=  1 NORB=', nbf0, ' ALPHA'
 do i = 1, nbf0, 1
  write(orbid,'(5(2X,E23.16))') (coeff(j,i),j=1,nbf0)
 end do ! for i

 if(uhf) then    ! UHF
  write(orbid,'(A,I9,A)') 'SYM=  1 NORB=', nbf0, ' BETA'
  do i = 1, nbf0, 1
   write(orbid,'(5(2X,E23.16))') (coeff(j,i+nbf0),j=1,nbf0)
  end do ! for i
 else if(mult/=1 .and. (.not.prt_no)) then ! ROHF
  write(orbid,'(A,I9,A)') 'SYM=  1 NORB=', nbf0, ' BETA'
  do i = 1, nbf0, 1
   write(orbid,'(5(2X,E23.16))') (coeff(j,i),j=1,nbf0)
  end do ! for i
 end if
 deallocate(coeff)

 write(orbid,'(A)') 'ORBITAL ENERGY'
 write(orbid,'(5(2X,E23.16))') (occ_num(i),i=1,nbf0)

 if(uhf) then
  write(orbid,'(5(2X,E23.16))') (occ_num(i),i=nbf0+1,2*nbf0)
 else if(mult/=1 .and. (.not.prt_no)) then
  write(orbid,'(5(2X,E23.16))') (occ_num(i),i=1,nbf0)
 end if

 write(orbid,'(A)') 'OCCUPATION'
 if(prt_no) then
  write(orbid,'(10(1X,E11.5))') (occ_num(i),i=1,nbf0)
  write(orbid,'(A)') '$END'
  close(orbid)
  return
 end if

 occ_num = 0d0
 occ_num(1:na) = 1d0
 write(orbid,'(10(1X,E11.5))') (occ_num(i),i=1,nbf0)

 if(mult/=1 .or. uhf) then
  if(nb < na) occ_num(nb+1:na) = 0d0
  write(orbid,'(A)') 'OCCUPATION'
  write(orbid,'(10(1X,E11.5))') (occ_num(i),i=1,nbf0)
 end if
 deallocate(occ_num)
 write(orbid,'(A)') '$END'

 call read_natom_from_fch(fchname, natom)
 write(orbid,'(A,1X,I8)') '$COORD', natom

 allocate(elem(natom), coor(3,natom), nuc(natom), ntimes(natom))
 call read_elem_and_coor_from_fch(fchname, natom, elem, nuc, coor, i, j)
 deallocate(nuc)
 call calc_ntimes(natom, elem, ntimes)
 coor = coor/Bohr_const

 do i = 1, natom, 1
  write(orbid,'(A,I0,3(1X,F19.12))') TRIM(elem(i)),ntimes(i),coor(1:3,i)
 end do ! for i
 deallocate(elem, ntimes, coor)

 write(orbid,'(A)') '$END'

 if(mult /= 1) then
  write(orbid,'(A,2(1X,I8),A)') '$NBF', nbf0, nbf0, '   2'
 else
  write(orbid,'(A,2(1X,I8),A)') '$NBF', nbf0, nbf0, '   1'
 end if

 write(orbid,'(A,2X,I7,1X,I7)') 'NOCC', na, nb
 write(orbid,'(A)') '$END'
 close(orbid)
! print done
end subroutine fch2bdf

! move the 2nd, 3rd, ... Zeta basis functions forward
subroutine zeta_mv_forwd(i0, shell_type, length, nbf, nif, coeff2)
 implicit none
 integer i, j, k
 integer, intent(in) :: i0, shell_type, length, nbf, nif
 integer, parameter :: num0(-5:5) = [11, 9, 7, 5, 0, 0, 3, 6, 10, 15, 21]
 !                                   11H 9G 7F 5D L  S 3P 6D 10F 15G 21H
 real(kind=8), intent(inout) :: coeff2(nbf,nif)
 real(kind=8), allocatable :: coeff(:,:)

 if(length == 1) return

 if(shell_type==0 .or. shell_type==-1) then
  write(6,'(A)') 'ERROR in subroutine zeta_mv_forwd: this element of&
                 & shell_type is 0 or -1. Impossible.'
  write(6,'(2(A,I0))') 'shell_type=', shell_type, ', length=', length
  stop
 end if

 coeff = coeff2
 k = num0(shell_type)
 do i = 1, k, 1
  do j = 1, length, 1
   coeff(i0+j+(i-1)*length,:) = coeff2(i0+i+(j-1)*k,:)
  end do ! for j
 end do ! for i

 coeff2 = coeff
 deallocate(coeff)
end subroutine zeta_mv_forwd

subroutine fch2bdf_permute_3p(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(3) = [2,3,1]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(3,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical p functions in Gaussian
! To: the order of spherical p functions in BDF
! 1    2    3
! Px,  Py,  Pz
! Py,  Pz,  Px

 allocate(coeff2(3,nif), source=0d0)
 forall(i = 1:3) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
end subroutine fch2bdf_permute_3p

subroutine fch2bdf_permute_5d(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(5) = [5, 3, 1, 2, 4]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(5,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical d functions in Gaussian
! To: the order of spherical d functions in BDF
! 1    2    3    4    5
! d0 , d+1, d-1, d+2, d-2
! d-2, d-1, d0 , d+1, d+2

 allocate(coeff2(5,nif), source=0d0)
 forall(i = 1:5) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
end subroutine fch2bdf_permute_5d

subroutine fch2bdf_permute_7f(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(7) = [7, 5, 3, 1, 2, 4, 6]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(7,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical f functions in Gaussian
! To: the order of spherical f functions in BDF
! 1    2    3    4    5    6    7
! f0 , f+1, f-1, f+2, f-2, f+3, f-3
! f-3, f-2, f-1, f0 , f+1, f+2, f+3

 allocate(coeff2(7,nif), source=0d0)
 forall(i = 1:7) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
end subroutine fch2bdf_permute_7f

subroutine fch2bdf_permute_9g(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(9) = [9, 7, 5, 3, 1, 2, 4, 6, 8]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(9,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical g functions in Gaussian
! To: the order of spherical g functions in BDF
! 1    2    3    4    5    6    7    8    9
! g0 , g+1, g-1, g+2, g-2, g+3, g-3, g+4, g-4
! g-4, g-3, g-2, g-1, g0 , g+1, g+2, g+3, g+4

 allocate(coeff2(9,nif), source=0d0)
 forall(i = 1:9) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
end subroutine fch2bdf_permute_9g

subroutine fch2bdf_permute_11h(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(11) = [11, 9, 7, 5, 3, 1, 2, 4, 6, 8, 10]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(11,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of spherical h functions in Gaussian
! To: the order of spherical h functions in BDF
! 1    2    3    4    5    6    7    8    9    10   11
! h0 , h+1, h-1, h+2, h-2, h+3, h-3, h+4, h-4, h+5, h-5
! h-5, h-4, h-3, h-2, h-1, h0 , h+1, h+2, h+3, h+4, h+5

 allocate(coeff2(11,nif), source=0d0)
 forall(i = 1:11) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
end subroutine fch2bdf_permute_11h

