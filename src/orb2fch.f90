! written by jxzou at 20200907: adjust the orders of d,f,g, etc functions in
!  OpenMolcas, to that in Gaussian .fch(k) file
! Originally copied from fch2inporb.f90, some modifications are made
! updated by jxzou at 20210413: remove '-uhf', add automatic determination

program main
 implicit none
 integer :: i
 character(len=3) :: str = ' '
 character(len=240) :: fchname, orbname
 logical :: prt_no

 i = iargc()
 if(.not. (i==2 .or. i==3)) then
  write(6,'(/,A)') ' ERROR in subroutine orb2fch: wrong command line arguments!'
  write(6,'(A)')   ' Example 1 (for RHF)   : orb2fch a.ScfOrb a.fch'
  write(6,'(A)')   ' Example 2 (for UHF)   : orb2fch a.UhfOrb a.fch'
  write(6,'(A)')   ' Example 3 (for CAS)   : orb2fch a.RasOrb a.fch'
  write(6,'(A)')   ' Example 4 (for UNO)   : orb2fch a.UnaOrb a.fch -no'
  write(6,'(A,/)') ' Example 5 (for CAS NO): orb2fch a.RasOrb.1 a.fch -no'
  stop
 end if

 fchname = ' '; orbname = ' '
 call getarg(1,orbname)
 call require_file_exist(orbname)

 call getarg(2,fchname)
 call require_file_exist(fchname)

 prt_no = .false.
 if(i == 3) then
  call getarg(3, str)
  if(str /= '-no') then
   write(6,'(/,A)') "ERROR in subroutine orb2fch: the 3rd argument is wrong! O&
                    &nly '-no' is accepted."
   stop
  else
   prt_no = .true.
  end if
 end if

 call orb2fch(orbname, fchname, prt_no)
end program main

! read the MOs in orbital file of OpenMolcas and adjust its d,f,g,h functions
!  order to that of Gaussian
subroutine orb2fch(orbname, fchname, prt_no)
 use fch_content, only: check_uhf_in_fch
 implicit none
 integer :: i, j, k, m, length
 integer :: na, nb, nbf, nif, nbf0, nbf1
 integer :: n6dmark, n10fmark, n15gmark, n21hmark
 integer :: n5dmark, n7fmark, n9gmark, n11hmark
 integer, allocatable :: shell_type(:), shell2atom_map(:)
 integer, allocatable :: idx(:), idx2(:)
 ! mark the index where d, f, g, h functions begin
 integer, allocatable :: d_mark(:), f_mark(:), g_mark(:), h_mark(:)
 character(len=240), intent(in) :: orbname, fchname
 ! orbname is one of .ScfOrb, .RasOrb, .RasOrb.1, .UnaOrb, .UhfOrb file of OpenMolcas
 real(kind=8), allocatable :: coeff(:,:), coeff2(:,:), occ_num(:), norm(:)
 logical :: uhf, sph
 logical, intent(in) :: prt_no

 call check_uhf_in_fch(fchname, uhf)
 call read_na_and_nb_from_fch(fchname, na, nb)
 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 nbf0 = nbf ! make a copy of nbf

 ! read MO Coefficients from .ScfOrb, .RasOrb, .RasOrb.1, .UnaOrb, or .UhfOrb
 ! file of OpenMolcas
 if(uhf) then
  allocate(coeff(nbf,2*nif))
  call read_mo_from_orb(orbname, nbf, nif, 'a', coeff(:,1:nif))
  call read_mo_from_orb(orbname, nbf, nif, 'b', coeff(:,nif+1:2*nif))
  nif = 2*nif   ! double the size
 else
  if(prt_no) then
   allocate(occ_num(nif), coeff(nbf,nif))
   call read_on_from_orb(orbname, nif, 'a', occ_num)
   call read_mo_from_orb(orbname, nbf, nif, 'a', coeff)
  else
   allocate(coeff(nbf,nif))
   call read_mo_from_orb(orbname, nbf, nif, 'a', coeff)
  end if
 end if

 call read_ncontr_from_fch(fchname, k)
 allocate(shell_type(2*k), source=0)
 allocate(shell2atom_map(2*k), source=0)
 call read_shltyp_and_shl2atm_from_fch(fchname, k, shell_type, shell2atom_map)

 ! check if any spherical functions
 if(ANY(shell_type<-1) .and. ANY(shell_type>1)) then
  write(6,'(A)') 'ERROR in subroutine orb2fch: mixed spherical harmonic/Cartes&
                 &ian functions detected.'
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
 allocate(norm(nbf), source=1d0)

 ! 2) sort the shell_type and shell2atom_map by ascending order,
 !  the indices of MOs will be adjusted accordingly
 call sort_shell_and_mo_idx(length, shell_type, shell2atom_map, nbf, idx)
! adjust done

! then we adjust the basis functions in each MO according to the type of basis functions
 k = length  ! update k
 allocate(d_mark(k), f_mark(k), g_mark(k), h_mark(k))

 ! adjust the order of d, f, etc. functions
 if(sph) then ! spherical harmonic
  call read_mark_from_shltyp_sph(k, shell_type, n5dmark, n7fmark, n9gmark, &
                                 n11hmark, d_mark, f_mark, g_mark, h_mark)
  call orb2fch_permute_sph(n5dmark, n7fmark, n9gmark, n11hmark, k, d_mark, &
                           f_mark, g_mark, h_mark, nbf, idx)
 else ! Cartesian-type basis
  call read_mark_from_shltyp_cart(k, shell_type, n6dmark, n10fmark, n15gmark,&
                                  n21hmark, d_mark, f_mark, g_mark, h_mark)
  call orb2fch_permute_cart(n6dmark, n10fmark, n15gmark, n21hmark, k, d_mark, &
                            f_mark, g_mark, h_mark, nbf, idx, norm)
 end if
! adjustment finished

 deallocate(d_mark, f_mark, g_mark, h_mark)

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
  case(-2)   ! 5D
   nbf = nbf + 5
  case( 2)   ! 6D
   nbf = nbf + 6
  case(-3)   ! 7F
   nbf = nbf + 7
  case( 3)   ! 10F
   nbf = nbf + 10
  case(-4)   ! 9G
   nbf = nbf + 9
  case( 4)   ! 15G
   nbf = nbf + 15
  case(-5)   ! 11H
   nbf = nbf + 11
  case( 5)   ! 21H
   nbf = nbf + 21
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
   call zeta_mv_forwd_idx(nbf1, m, length, nbf0, idx, norm)
   nbf = nbf1 + length*(nbf-nbf1)
  end if
 end do ! for while

 deallocate(shell_type, shell2atom_map)
! move done

 nbf = nbf0
 allocate(idx2(nbf), coeff2(nbf,nif))
 forall(i = 1:nbf) idx2(idx(i)) = i
 forall(i=1:nbf, j=1:nif) coeff2(i,j) = coeff(idx2(i),j)*norm(idx2(i))
 deallocate(idx, idx2, norm, coeff)

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
end subroutine orb2fch

