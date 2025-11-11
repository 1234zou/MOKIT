! written by jxzou at 20210107: adjust the orders of p,d,f,g, etc functions in
!  BDF, to that in Gaussian .fch(k) file
! Originally copied from bdf2fch.f90, some modifications are made
! updated by jxzou at 20210111: transfer orbital energies when transferring HF/DFT orbitals
! updated by jxzou at 20230316: add an optional argument new_fch

program main
 implicit none
 integer :: i
 character(len=3) :: str = ' '
 character(len=240) :: orbname, fchname, new_fch
 logical :: prt_no

 i = iargc()
 if(i<2 .or. i>4) then
  write(6,'(/,A)') ' ERROR in program bdf2fch: wrong command line arguments!'
  write(6,'(A)')   ' Example 1 (R(O)HF/UHF): bdf2fch a.scforb a.fch'
  write(6,'(A)')   ' Example 2 (R(O)HF/UHF): bdf2fch a.scforb a.fch a_new.fch'
  write(6,'(A)')   ' Example 3 (CAS)       : bdf2fch a.casorb a.fch'
  write(6,'(A)')   ' Example 4 (CAS)       : bdf2fch a.casorb a.fch a_new.fch'
  write(6,'(A)')   ' Example 5 (CAS NO)    : bdf2fch a.casorb a.fch -no'
  write(6,'(A,/)') ' Example 6 (CAS NO)    : bdf2fch a.casorb a.fch a_new.fch -no'
  stop
 end if

 str = ' '; orbname = ' '; fchname = ' '; new_fch = ' '
 call getarg(1, orbname)
 call require_file_exist(orbname)

 call getarg(2, fchname)
 call require_file_exist(fchname)

 prt_no = .false.
 if(i == 3) then
  call getarg(3, new_fch)
  if(TRIM(new_fch) == '-no') then
   prt_no = .true.
   new_fch = ' ' ! update new_fch
  end if
 end if

 if(i == 4) then
  call getarg(4, str)
  if(str /= '-no') then
   write(6,'(/,A)') "ERROR in subroutine bdf2fch: the 4th argument is wrong! O&
                    &nly '-no' is accepted."
   write(6,'(A)') "But you specify '"//str//"'."
   stop
  end if
  prt_no = .true.
 end if

 call bdf2fch(orbname, fchname, new_fch, prt_no)
end program main

! read the MOs in orbital file of BDF and adjust its p,d,f,g,h functions
!  order to that of Gaussian
subroutine bdf2fch(orbname, fchname, new_fch, prt_no)
 implicit none
 integer :: i, j, k, m, length
 integer :: na, nb, nbf, nif, nbf0, nbf1
 integer :: n3pmark, n5dmark, n7fmark, n9gmark, n11hmark
 integer, allocatable :: shell_type(:), shell2atom_map(:)
 integer, allocatable :: idx(:), idx2(:)
 ! mark the index where p, d, f, g, h functions begin
 integer, allocatable :: p_mark(:), d_mark(:), f_mark(:), g_mark(:), h_mark(:)
 character(len=240), intent(in) :: orbname, fchname
 ! orbname: one of .scforb, .casorb, .casorb1, canorb, file of BDF
 character(len=240), intent(inout) :: new_fch
 ! If the user provides the filename new_fch, use the filename; if not, replace
 ! orbitals in fchname
 real(kind=8), allocatable :: coeff(:,:), coeff2(:,:)
 real(kind=8), allocatable :: occ_num(:), ev(:)
 logical :: uhf
 logical, intent(in) :: prt_no

 call check_uhf_in_fch(fchname, uhf) ! determine whether UHF or not
 call read_na_and_nb_from_fch(fchname, na, nb)
 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 nbf0 = nbf ! make a copy of nbf

 ! read MO Coefficients from orbital file of BDF
 if(uhf) then
  allocate(coeff(nbf,2*nif), ev(2*nif))
  call read_ev_from_bdf_orb(orbname, nif, 'a', ev(1:nif))
  call read_ev_from_bdf_orb(orbname, nif, 'b', ev(nif+1:2*nif))
  call read_mo_from_bdf_orb(orbname, nbf, nif, 'a', coeff(:,1:nif))
  call read_mo_from_bdf_orb(orbname, nbf, nif, 'b', coeff(:,nif+1:2*nif))
  nif = 2*nif   ! double the size
 else
  if(prt_no) then
   allocate(coeff(nbf,nif), occ_num(nif))
   call read_on_from_bdf_orb(orbname, nif, 'a', occ_num)
   call read_mo_from_bdf_orb(orbname, nbf, nif, 'a', coeff)
  else
   allocate(coeff(nbf,nif), ev(nif))
   call read_ev_from_bdf_orb(orbname, nif, 'a', ev)
   call read_mo_from_bdf_orb(orbname, nbf, nif, 'a', coeff)
  end if
 end if

 call read_ncontr_from_fch(fchname, k)
 allocate(shell_type(2*k), source=0)
 allocate(shell2atom_map(2*k), source=0)
 call read_shltyp_and_shl2atm_from_fch(fchname, k, shell_type, shell2atom_map)

 if(ANY(shell_type>1)) then
  write(6,'(A)') 'ERROR in subroutine bdf2fch: Cartesian-type basis functions&
                 & detected. Cannot deal with'
  write(6,'(A)') 'that. BDF supports only spherical harmonic basis functions.'
  write(6,'(A)') "You should add keywords '5D 7F' in Gaussian input file&
                 & to obtain a new .fch(k) file."
  stop
 end if

! first we adjust the basis functions in each MO according to the Shell to atom map
 ! 1) split the 'L' into 'S' and 'P', this is to ensure that D comes after L functions
 call split_L_func(k, shell_type, shell2atom_map, length)

 allocate(idx(nbf))
 forall(i = 1:nbf) idx(i) = i

 ! 2) sort the shell_type and shell2atom_map by ascending order,
 !  the indices of MOs will be adjusted accordingly
 call sort_shell_and_mo_idx(length, shell_type, shell2atom_map, nbf, idx)
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
  end select
 end do ! for i

 ! adjust the order of p, d, f, etc. functions
 do i = 1, n3pmark, 1
  call bdf2fch_permute_3p(idx(p_mark(i):p_mark(i)+2))
 end do
 do i = 1, n5dmark, 1
  call bdf2fch_permute_5d(idx(d_mark(i):d_mark(i)+4))
 end do
 do i = 1, n7fmark, 1
  call bdf2fch_permute_7f(idx(f_mark(i):f_mark(i)+6))
 end do
 do i = 1, n9gmark, 1
  call bdf2fch_permute_9g(idx(g_mark(i):g_mark(i)+8))
 end do
 do i = 1, n11hmark, 1
  call bdf2fch_permute_11h(idx(h_mark(i):h_mark(i)+10))
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
  end do ! for while

  if(i < k) i = i - 1
  if(length > 1) then
   call zeta_mv_forwd(nbf1, m, length, nbf0, idx)
   nbf = nbf1 + length*(nbf-nbf1)
  end if
 end do ! for while

 deallocate(shell_type, shell2atom_map)
! move done

 nbf = nbf0
 allocate(idx2(nbf), coeff2(nbf,nif))
 forall(i = 1:nbf) idx2(idx(i)) = i
 forall(i=1:nbf, j=1:nif) coeff2(i,j) = coeff(idx2(i),j)
 deallocate(idx, idx2, coeff)

 if(LEN_TRIM(new_fch) == 0) then
  new_fch = fchname
 else
  call copy_file(fchname, new_fch, .false.)
 end if

 ! print MOs and orbital energies/occupation numbers into .fch(k) file
 if(uhf) then
  nif = nif/2
  call write_eigenvalues_to_fch(new_fch, nif, 'a', ev(1:nif), .true.)
  call write_eigenvalues_to_fch(new_fch, nif, 'b', ev(nif+1:2*nif), .true.)
  call write_mo_into_fch(new_fch, nbf, nif, 'a', coeff2(:,1:nif))
  call write_mo_into_fch(new_fch, nbf, nif, 'b', coeff2(:,nif+1:2*nif))
 else
  if(prt_no) then
   call write_eigenvalues_to_fch(new_fch, nif, 'a', occ_num, .true.)
   call write_mo_into_fch(new_fch, nbf, nif, 'a', coeff2)
  else
   call write_eigenvalues_to_fch(new_fch, nif, 'a', ev, .true.)
   call write_mo_into_fch(new_fch, nbf, nif, 'a', coeff2)
  end if
 end if

 deallocate(coeff2)
 if(allocated(occ_num)) deallocate(occ_num)
end subroutine bdf2fch

! move the 2nd, 3rd, ... Zeta basis functions forward
subroutine zeta_mv_forwd(i0, shell_type, length, nbf, idx2)
 implicit none
 integer :: i, j, k
 integer, intent(in) :: i0, shell_type, length, nbf
 integer, parameter :: num0(-5:5) = [11, 9, 7, 5, 0, 0, 3, 6, 10, 15, 21]
 !                                   11H 9G 7F 5D L  S 3P 6D 10F 15G 21H
 integer, intent(inout) :: idx2(nbf)
 integer, allocatable :: idx(:)

 if(length == 1) return

 if(shell_type==0 .or. shell_type==-1) then
  write(6,'(A)') 'ERROR in subroutine zeta_mv_forwd: this element of&
                   & shell_type is 0 or -1. Impossible.'
  write(6,'(2(A,I0))') 'shell_type=', shell_type, ', length=', length
  stop
 end if

 idx = idx2
 k = num0(shell_type)

 do i = 1, k, 1
  do j = 1, length, 1
   idx(i0+j+(i-1)*length) = idx2(i0+i+(j-1)*k)
  end do ! for j
 end do ! for i

 idx2 = idx
 deallocate(idx)
end subroutine zeta_mv_forwd

subroutine bdf2fch_permute_3p(idx)
 implicit none
 integer :: i, idx0(3)
 integer, parameter :: order(3) = [2,3,1]
 integer, intent(inout) :: idx(3)
! From: the order of spherical p functions in Gaussian
! To: the order of spherical p functions in BDF
! 1    2    3
! Px,  Py,  Pz
! Py,  Pz,  Px

 idx0 = idx
 forall(i = 1:3) idx(i) = idx0(order(i))
end subroutine bdf2fch_permute_3p

subroutine bdf2fch_permute_5d(idx)
 implicit none
 integer :: i, idx0(5)
 integer, parameter :: order(5) = [5, 3, 1, 2, 4]
 integer, intent(inout) :: idx(5)
! From: the order of spherical d functions in Gaussian
! To: the order of spherical d functions in BDF
! 1    2    3    4    5
! d0 , d+1, d-1, d+2, d-2
! d-2, d-1, d0 , d+1, d+2

 idx0 = idx
 forall(i = 1:5) idx(i) = idx0(order(i))
end subroutine bdf2fch_permute_5d

subroutine bdf2fch_permute_7f(idx)
 implicit none
 integer :: i, idx0(7)
 integer, parameter :: order(7) = [7, 5, 3, 1, 2, 4, 6]
 integer, intent(inout) :: idx(7)
! From: the order of spherical f functions in Gaussian
! To: the order of spherical f functions in BDF
! 1    2    3    4    5    6    7
! f0 , f+1, f-1, f+2, f-2, f+3, f-3
! f-3, f-2, f-1, f0 , f+1, f+2, f+3

 idx0 = idx
 forall(i = 1:7) idx(i) = idx0(order(i))
end subroutine bdf2fch_permute_7f

subroutine bdf2fch_permute_9g(idx)
 implicit none
 integer :: i, idx0(9)
 integer, parameter :: order(9) = [9, 7, 5, 3, 1, 2, 4, 6, 8]
 integer, intent(inout) :: idx(9)
! From: the order of spherical g functions in Gaussian
! To: the order of spherical g functions in BDF
! 1    2    3    4    5    6    7    8    9
! g0 , g+1, g-1, g+2, g-2, g+3, g-3, g+4, g-4
! g-4, g-3, g-2, g-1, g0 , g+1, g+2, g+3, g+4

 idx0 = idx
 forall(i = 1:9) idx(i) = idx0(order(i))
end subroutine bdf2fch_permute_9g

subroutine bdf2fch_permute_11h(idx)
 implicit none
 integer :: i, idx0(11)
 integer, parameter :: order(11) = [11, 9, 7, 5, 3, 1, 2, 4, 6, 8, 10]
 integer, intent(inout) :: idx(11)
! From: the order of spherical h functions in Gaussian
! To: the order of spherical h functions in BDF
! 1    2    3    4    5    6    7    8    9    10   11
! h0 , h+1, h-1, h+2, h-2, h+3, h-3, h+4, h-4, h+5, h-5
! h-5, h-4, h-3, h-2, h-1, h0 , h+1, h+2, h+3, h+4, h+5

 idx0 = idx
 forall(i = 1:11) idx(i) = idx0(order(i))
end subroutine bdf2fch_permute_11h

