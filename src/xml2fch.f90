! written by jxzou at 20201214: adjust the orders of d,f,g, etc functions in
!  Molpro, to that in Gaussian .fch(k) file
! Originally copied from orb2fch.f90, some modifications are made
! updated by jxzou at 20210412: remove '-uhf', add automatic determination

program main
 implicit none
 integer :: i
 character(len=3) :: str
 character(len=240) :: fchname, xmlname
 logical :: prt_no

 i = iargc()
 if(.not. (i==2 .or. i==3)) then
  write(6,'(/,A)') ' ERROR in subroutine xml2fch: wrong command line arguments!'
  write(6,'(A)')   ' Example 1 (for R(O)HF, UHF): xml2fch a.xml a.fch'
  write(6,'(A,/)') ' Example 2 (for CAS NO)     : xml2fch a.xml a.fch -no'
  stop
 end if

 fchname = ' '
 call getarg(1,xmlname)
 call require_file_exist(xmlname)
 call getarg(2,fchname)
 call require_file_exist(fchname)
 prt_no = .false.

 if(i == 3) then
  call getarg(3, str)
  if(str /= '-no') then
   write(6,'(/,1X,A)') "ERROR in subroutine xml2fch: the 3rd argument is&
                       & wrong! Only '-no' is accepted."
   stop
  else
   prt_no = .true.
  end if
 end if

 call xml2fch(xmlname, fchname, prt_no)
end program main

! nbf: the number of basis functions
! nif: the number of independent functions, i.e., the number of MOs

! read the MOs in orbital file of Molpro and adjust its d,f,g,h functions
!  order to that of Gaussian
subroutine xml2fch(xmlname, fchname, prt_no)
 implicit none
 integer :: i, j, k, length, na, nb, nbf, nif, nbf0
 integer :: n6dmark, n10fmark, n15gmark, n21hmark, n28imark
 integer :: n5dmark, n7fmark, n9gmark, n11hmark, n13imark
 integer, allocatable :: shell_type(:), shell2atom_map(:)
 integer, allocatable :: idx(:), idx2(:), d_mark(:), f_mark(:), g_mark(:), &
  h_mark(:), i_mark(:)
 character(len=240), intent(in) :: xmlname, fchname
 real(kind=8), allocatable :: coeff(:,:), coeff2(:,:), occ_num(:)
 logical :: uhf, sph
 logical, intent(in) :: prt_no

 call check_uhf_in_fch(fchname, uhf)
 call read_na_and_nb_from_fch(fchname, na, nb)
 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 nbf0 = nbf ! make a copy of nbf

 ! read MO Coefficients from Molpro .xml file
 if(uhf) then
  allocate(coeff(nbf,2*nif))
  call read_mo_from_xml(xmlname, nbf, nif, 'a', coeff(:,1:nif))
  call read_mo_from_xml(xmlname, nbf, nif, 'b', coeff(:,nif+1:2*nif))
  nif = 2*nif   ! double the size
 else ! not UHF
  if(prt_no) then
   allocate(occ_num(nif), coeff(nbf,nif))
   call read_on_from_xml(xmlname, nif, 'a', occ_num)
   call read_mo_from_xml(xmlname, nbf, nif, 'a', coeff)
  else
   allocate(coeff(nbf,nif))
   call read_mo_from_xml(xmlname, nbf, nif, 'a', coeff)
  end if
 end if

 call read_ncontr_from_fch(fchname, k)
 allocate(shell_type(2*k), source=0)
 allocate(shell2atom_map(2*k), source=0)
 call read_shltyp_and_shl2atm_from_fch(fchname, k, shell_type, shell2atom_map)

 if(ANY(shell_type<-1) .and. ANY(shell_type>1)) then
  write(6,'(A)') 'ERROR in subroutine xml2fch: mixed spherical harmonic/&
                 &Cartesian functions detected.'
  write(6,'(A)') 'You probably used a basis set like 6-31G(d) in Gaussian. Its&
                 & default setting is (6D,7F).'
  write(6,'(A)') "You need to add '5D 7F' or '6D 10F' keywords in Gaussian."
  stop
 else if( ANY(shell_type>1) ) then
  sph = .false.
 else
  sph = .true.
 end if

! first we adjust the basis functions in each MO according to the Shell to atom map
 ! 1) split the 'L' into 'S' and 'P', this is to ensure that D comes after L functions
 call split_L_func(k, shell_type, shell2atom_map, length)

 allocate(idx(nbf))
 forall(i = 1:nbf) idx(i) = i

 ! 2) sort the shell_type and shell2atom_map by ascending order,
 !  the indices of MOs will be adjusted accordingly
 call sort_shell_and_mo_idx(length, shell_type, shell2atom_map, nbf, idx)
 deallocate(shell2atom_map)
! adjust done

! then we adjust the basis functions in each MO according to the type of basis functions
 k = length  ! update k
 allocate(d_mark(k), f_mark(k), g_mark(k), h_mark(k))

 ! adjust the order of d, f, etc. functions
 if(sph) then ! spherical harmonic
  call read_mark_from_shltyp_sph(k, shell_type, n5dmark, n7fmark, n9gmark, &
                 n11hmark, n13imark, d_mark, f_mark, g_mark, h_mark, i_mark)
  do i = 1, n5dmark, 1
   call xml2fch_permute_5d(idx(d_mark(i):d_mark(i)+4))
  end do
  do i = 1, n7fmark, 1
   call xml2fch_permute_7f(idx(f_mark(i):f_mark(i)+6))
  end do
  do i = 1, n9gmark, 1
   call xml2fch_permute_9g(idx(g_mark(i):g_mark(i)+8))
  end do
  do i = 1, n11hmark, 1
   call xml2fch_permute_11h(idx(h_mark(i):h_mark(i)+10))
  end do
 else  ! Cartesian-type basis
  call read_mark_from_shltyp_cart(k, shell_type, n6dmark, n10fmark, n15gmark, &
                    n21hmark, n28imark, d_mark, f_mark, g_mark, h_mark, i_mark)
  do i = 1, n10fmark, 1
   call xml2fch_permute_10f(idx(f_mark(i):f_mark(i)+9))
  end do
  do i = 1, n15gmark, 1
   call xml2fch_permute_15g(idx(g_mark(i):g_mark(i)+14))
  end do
  do i = 1, n21hmark, 1
   call xml2fch_permute_21h(idx(h_mark(i):h_mark(i)+20))
  end do
 end if
! adjustment finished

 deallocate(shell_type, d_mark, f_mark, g_mark, h_mark, i_mark)

 nbf = nbf0
 allocate(idx2(nbf), coeff2(nbf,nif))
 forall(i = 1:nbf) idx2(idx(i)) = i
 forall(i=1:nbf, j=1:nif) coeff2(i,j) = coeff(idx2(i),j)
 deallocate(idx, idx2, coeff)

! print MOs into .fch(k) file
 if(uhf) then
  nif = nif/2
  call write_mo_into_fch(fchname, nbf, nif, 'a', coeff2(:,1:nif))
  call write_mo_into_fch(fchname, nbf, nif, 'b', coeff2(:,nif+1:2*nif))
 else
  if(prt_no) then
   call write_eigenvalues_to_fch(fchname, nif, 'a', occ_num, .true.)
   call write_mo_into_fch(fchname, nbf, nif, 'a', coeff2)
  else
   call write_mo_into_fch(fchname, nbf, nif, 'a', coeff2)
  end if
 end if
! print done

 deallocate(coeff2)
 if(allocated(occ_num)) deallocate(occ_num)
end subroutine xml2fch

subroutine xml2fch_permute_5d(idx)
 implicit none
 integer :: i, idx0(5)
 integer, parameter :: order(5) = [1,5,2,4,3]
 integer, intent(inout) :: idx(5)

 idx0 = idx
 forall(i = 1:5) idx(i) = idx0(order(i))
end subroutine xml2fch_permute_5d

subroutine xml2fch_permute_7f(idx)
 implicit none
 integer :: i, idx0(7)
 integer, parameter :: order(7) = [2,3,1,6,5,7,4]
 integer, intent(inout) :: idx(7)

 idx0 = idx
 forall(i = 1:7) idx(i) = idx0(order(i))
end subroutine xml2fch_permute_7f

subroutine xml2fch_permute_10f(idx)
 implicit none
 integer :: i, idx0(10)
 integer, parameter :: order(10) = [1,2,3,5,6,4,9,7,8,10]
 integer, intent(inout) :: idx(10)

 idx0 = idx
 forall(i = 1:10) idx(i) = idx0(order(i))
end subroutine xml2fch_permute_10f

subroutine xml2fch_permute_9g(idx)
 implicit none
 integer :: i, idx0(9)
 integer, parameter :: order(9) = [1,5,2,8,3,4,9,6,7]
 integer, intent(inout) :: idx(9)

 idx0 = idx
 forall(i = 1:9) idx(i) = idx0(order(i))
end subroutine xml2fch_permute_9g

subroutine xml2fch_permute_15g(idx)
 implicit none
 integer :: i, idx0(15)
 integer, parameter :: order(15) = [15,5,1,14,13,9,4,6,2,12,10,3,11,8,7]
 integer, intent(inout) :: idx(15)

 idx0 = idx
 forall(i = 1:15) idx(i) = idx0(order(i))
end subroutine xml2fch_permute_15g

subroutine xml2fch_permute_11h(idx)
 implicit none
 integer :: i, idx0(11)
 integer, parameter :: order(11) = [2,3,4,6,9,7,8,11,1,10,5]
 integer, intent(inout) :: idx(11)

 idx0 = idx
 forall(i = 1:11) idx(i) = idx0(order(i))
end subroutine xml2fch_permute_11h

subroutine xml2fch_permute_21h(idx)
 implicit none
 integer :: i, idx0(21)
 integer, intent(inout) :: idx(21)

 idx0 = idx
 forall(i = 1:21) idx(i) = idx0(22-i)
end subroutine xml2fch_permute_21h

