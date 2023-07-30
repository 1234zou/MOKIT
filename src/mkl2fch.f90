! written by jxzou at 20200627: transfer MOs from .mkl to .fch(k)
! updated by jxzou at 20201213: read NOONs from .mkl and save to .fch
! updated by jxzou at 20210413: remove '-uhf', add automatic determination

! The 'Shell types' array in Gaussian .fch file:
!
!   Spherical     |     Cartesian
! -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5
!  H  G  F  D  L  S  P  D  F  G  H
!
! 'L' is 'SP' in Pople-type basis sets

program main
 implicit none
 integer :: i, j, no_type
 character(len=4) :: str = ' '
 character(len=240) :: mklname, fchname
 logical :: alive

 i = iargc()
 if(i<1 .or. i>3) then
  write(6,'(/,A)') ' ERROR in subroutine mkl2fch: wrong command line arguments!'
  write(6,'(A)')   ' Example 1 (R(O)HF, UHF, CAS)    : mkl2fch a.mkl'
  write(6,'(A)')   ' Example 2 (R(O)HF, UHF, CAS)    : mkl2fch a.mkl a.fch'
  write(6,'(A)')   ' Example 3 (UNO, CAS NO)         : mkl2fch a.mkl a.fch -no'
  write(6,'(A,/)') ' Example 4 (UHF, UMP2, UCCSD NSO): mkl2fch a.mkl a.fch -nso'
  stop
 end if

 mklname = ' '; fchname = ' '
 no_type = 0   ! default, canonical orbitals assumed

 call getarg(1, mklname)
 call require_file_exist(mklname)

 if(i == 1) then
  call find_specified_suffix(mklname, '.mkl', j)
  fchname = mklname(1:j-1)//'.fch'
 else
  call getarg(2, fchname)
  call find_specified_suffix(fchname, '.fch', j)
 end if

 if(i > 2) then
  call getarg(3, str)

  select case(TRIM(str))
  case('-no')  ! (spatial) natural orbtials
   no_type = 1
  case('-nso') ! natural spin orbitals
   no_type = 2
  case default
   write(6,'(/,A)') "ERROR in subroutine mkl2fch: the third argument can only&
                   & be '-no' or '-sno'."
   write(6,'(A)') "But you specify '"//TRIM(str)//"'."
   stop
  end select
 end if

 inquire(file=TRIM(fchname),exist=alive)

 if(alive) then
  write(6,'(/,A)') 'Remark from program mkl2fch: fchname '//TRIM(fchname)//' ex&
                   &ists.'
  write(6,'(A)') 'It will be used directly.'
  call mkl2fch(mklname, fchname, no_type)
 else
  write(6,'(/,A)') 'Warning from program mkl2fch: fchname '//TRIM(fchname)//' d&
                   &oes not'
  write(6,'(A)') 'exist, the program mkl2fch is trying to generate one from scr&
                 &atch...'
  call mkl2fch_direct(mklname, fchname, no_type)
  write(6,'(/,A,/)') 'Done generation.'
 end if
end program main

! convert .fch(k) file (Gaussian) to .mkl file (Molekel, ORCA)
subroutine mkl2fch(mklname, fchname, no_type)
 use fch_content
 use mkl_content, only: read_mo_from_mkl, read_on_from_mkl, read_ev_from_mkl
 implicit none
 integer :: i, k, nf3mark, ng3mark, nh3mark
 integer, intent(in) :: no_type
 integer, allocatable :: f3_mark(:), g3_mark(:), h3_mark(:)
 real(kind=8), allocatable :: noon(:), coeff(:,:)
 character(len=240), intent(in) :: mklname, fchname
 logical :: uhf

 call find_specified_suffix(fchname, '.fch', i)
 call check_uhf_in_fch(fchname, uhf) ! determine whether UHF
 call read_fch(fchname, uhf)         ! read content in .fch(k) file

 if(no_type==2 .and. (.not.uhf)) then
  write(6,'(/,A)') 'ERROR in subroutine mkl2fch: Natural Spin Orbitals requeste&
                   &d. But this is'
  write(6,'(A)') 'not a UHF-type .fch file.'
  stop
 end if

 ! check if any Cartesian functions
 if( ANY(shell_type > 1) ) then
  write(6,'(/,A)') 'ERROR in subroutine mkl2fch: Cartesian-type basis functions&
                   & detected in file'
  write(6,'(A)') TRIM(fchname)//'.'
  write(6,'(A)') "ORCA supports spherical functions only. You need to add '5D 7&
                 &F' keywords in"
  write(6,'(A)') 'Gaussian input file.'
  stop
 else if( ANY(shell_type < -5) ) then
  write(6,'(/,A)') 'ERROR in subroutine mkl2fch: angular momentum too high! not&
                   & supported.'
  stop
 end if
 ! check done

 ! read Alpha and/or Beta MOs from .fch(k) file
 call read_mo_from_mkl(mklname, nbf, nif, 'a', alpha_coeff)
 if(uhf) then ! UHF
  call read_mo_from_mkl(mklname, nbf, nif, 'b', beta_coeff)
  k = 2*nif
  allocate(coeff(nbf,k))
  coeff(:,1:nif) = alpha_coeff
  coeff(:,nif+1:) = beta_coeff
 else         ! R(O)HF
  k = nif
  allocate(coeff(nbf,k), source=alpha_coeff)
 end if

 ! find F+3, G+3 and H+3 functions, multiply them by -1
 allocate(f3_mark(nbf), g3_mark(nbf), h3_mark(nbf))
 call get_bas_mark_from_shltyp(ncontr, shell_type, nbf, nf3mark, ng3mark, &
                               nh3mark, f3_mark, g3_mark, h3_mark)
 deallocate(shell_type)
 call update_mo_using_bas_mark(nbf, k, nf3mark, ng3mark, nh3mark, f3_mark, &
                               g3_mark, h3_mark, coeff)
 deallocate(f3_mark, g3_mark, h3_mark)

 if(uhf) then ! UHF
  alpha_coeff = coeff(:,1:nif)
  beta_coeff = coeff(:,nif+1:)
 else         ! R(O)HF
  alpha_coeff = coeff
 end if
 deallocate(coeff)

 allocate(noon(nif))
 select case(no_type)
 case(0) ! canonical orbitals
  call read_ev_from_mkl(mklname, nif, 'a', noon)
  call write_eigenvalues_to_fch(fchname, nif, 'a', noon, .true.)
  if(uhf) then
   call read_ev_from_mkl(mklname, nif, 'b', noon)
   call write_eigenvalues_to_fch(fchname, nif, 'b', noon, .true.)
  end if
 case(1) ! natural orbitals
  call read_on_from_mkl(mklname, nif, 'a', noon)
  call write_eigenvalues_to_fch(fchname, nif, 'a', noon, .true.)
 case(2) ! natural spin orbitals
  call read_on_from_mkl(mklname, nif, 'a', noon)
  call write_eigenvalues_to_fch(fchname, nif, 'a', noon, .true.)
  call read_on_from_mkl(mklname, nif, 'b', noon)
  call write_eigenvalues_to_fch(fchname, nif, 'b', noon, .true.)
 case default
  write(6,'(/,A)') 'ERROR in subroutine mkl2fch: no_type out of range!'
  write(6,'(A,I0)') 'no_type=', no_type
  stop
 end select

 call write_mo_into_fch(fchname, nbf, nif, 'a', alpha_coeff)
 deallocate(alpha_coeff)

 if(uhf) then
  call write_mo_into_fch(fchname, nbf, nif, 'b', beta_coeff)
  deallocate(beta_coeff)
 end if

 deallocate(noon)
end subroutine mkl2fch

! Convert a ORCA .mkl file to a Gaussian fch file. The fch file will be generated
! from scratch
subroutine mkl2fch_direct(mklname, fchname, no_type)
 use fch_content
 use mkl_content, only: check_uhf_in_mkl, read_mkl, read_on_from_mkl, nuc, &
  charge0=>charge, mult0=>mult, natom0=>natom, ncontr0=>ncontr, nbf0=>nbf, &
  nif0=>nif, shell_type0=>shell_type, shl2atm, alpha_coeff0=>alpha_coeff, &
  beta_coeff0=>beta_coeff, elem0=>elem, coor0=>coor, all_pg, ev_a, ev_b
 implicit none
 integer :: i, j, k, nc, ne, nline, ncol, nf3mark, ng3mark, nh3mark
 integer, intent(in) :: no_type
 integer, allocatable :: f3_mark(:), g3_mark(:), h3_mark(:)
 real(kind=8), allocatable :: coeff(:,:)
 character(len=240), intent(in) :: mklname, fchname
 logical :: has_sp = .false.

 call check_uhf_in_mkl(mklname, is_uhf)
 call read_mkl(mklname, is_uhf, .true.)
 deallocate(nuc)

 charge = charge0
 mult = mult0
 natom = natom0
 ncontr = ncontr0
 nbf = nbf0
 nif = nif0

 allocate(ielem(natom))
 forall(i = 1:natom) ielem(i) = elem2nuc(elem0(i))
 deallocate(elem0)
 ne = SUM(ielem) - charge
 nopen = mult - 1
 na = (ne + nopen)/2
 nb = ne - na

 allocate(coor(3,natom), source=coor0)
 deallocate(coor0)

 allocate(shell_type(ncontr), shell2atom_map(ncontr), prim_per_shell(ncontr))
 shell_type = shell_type0
 shell2atom_map = shl2atm
 deallocate(shell_type0, shl2atm)

 ! find nprim from type all_pg
 k = 0; nprim = 0
 do i = 1, natom, 1
  nc = all_pg(i)%nc
  do j = 1, nc, 1
   nline = all_pg(i)%prim_gau(j)%nline
   prim_per_shell(k+1) = nline
   k = k + 1
   nprim = nprim + nline
   if(.not. has_sp) then
    if(all_pg(i)%prim_gau(j)%ncol == 3) has_sp = .true.
   end if
  end do ! for j
 end do ! for i

 allocate(prim_exp(nprim), source=0d0)
 allocate(contr_coeff(nprim), source=0d0)
 if(has_sp) allocate(contr_coeff_sp(nprim), source=0d0)

 ! convert type all_pg to arrays prim_exp, contr_coeff, and contr_coeff_sp
 k = 0
 do i = 1, natom, 1
  nc = all_pg(i)%nc
  do j = 1, nc, 1
   nline = all_pg(i)%prim_gau(j)%nline
   ncol = all_pg(i)%prim_gau(j)%ncol
   prim_exp(k+1:k+nline) = all_pg(i)%prim_gau(j)%coeff(:,1)
   contr_coeff(k+1:k+nline) = all_pg(i)%prim_gau(j)%coeff(:,2)
   if(ncol==3) contr_coeff_sp(k+1:k+nline) = all_pg(i)%prim_gau(j)%coeff(:,3)
   k = k + nline
  end do ! for j
 end do ! for i

 deallocate(all_pg)

 ! update MO coefficients
 if(is_uhf) then ! UHF
  k = 2*nif
  allocate(coeff(nbf,k))
  coeff(:,1:nif) = alpha_coeff0
  coeff(:,nif+1:) = beta_coeff0
  deallocate(alpha_coeff0, beta_coeff0)
 else            ! R(O)HF
  k = nif
  allocate(coeff(nbf,k), source=alpha_coeff0)
  deallocate(alpha_coeff0)
 end if

 ! find F+3, G+3 and H+3 functions, multiply them by -1
 allocate(f3_mark(nbf), g3_mark(nbf), h3_mark(nbf))
 call get_bas_mark_from_shltyp(ncontr, shell_type, nbf, nf3mark, ng3mark, &
                               nh3mark, f3_mark, g3_mark, h3_mark)
 call update_mo_using_bas_mark(nbf, k, nf3mark, ng3mark, nh3mark, f3_mark, &
                               g3_mark, h3_mark, coeff)
 deallocate(f3_mark, g3_mark, h3_mark)

 if(is_uhf) then ! UHF
  allocate(alpha_coeff(nbf,nif), source=coeff(:,1:nif))
  allocate(beta_coeff(nbf,nif), source=coeff(:,nif+1:))
 else            ! R(O)HF
  allocate(alpha_coeff(nbf,nif), source=coeff)
 end if
 deallocate(coeff)
 ! update MO coefficients done

 ! generate arrays tot_dm (and spin_dm)
 allocate(eigen_e_a(nif), tot_dm(nbf,nbf))
 select case(no_type)
 case(0) ! canonical orbitals
  call check_na_nb_ecp_in_mkl(mklname, is_uhf, nif, ne, na, nb)
  eigen_e_a = ev_a
  ev_a = 0d0
  if(is_uhf) then ! UHF
   ev_a(1:na) = 1d0
   call calc_dm_using_mo_and_on(nbf, nif, alpha_coeff, ev_a, tot_dm)
   deallocate(ev_a)
   allocate(eigen_e_b(nif), source=ev_b)
   ev_b = 0d0; ev_b(1:nb) = 1d0
   allocate(spin_dm(nbf,nbf))
   call calc_dm_using_mo_and_on(nbf, nif, beta_coeff, ev_b, spin_dm)
   deallocate(ev_b)
   allocate(coeff(nbf,nbf), source=tot_dm)
   tot_dm = coeff + spin_dm
   spin_dm = coeff - spin_dm
   deallocate(coeff)
  else            ! R(O)HF
   ev_a(1:nb) = 2d0
   if(na > nb) ev_a(nb+1:na) = 1d0
   call calc_dm_using_mo_and_on(nbf, nif, alpha_coeff, ev_a, tot_dm)
   deallocate(ev_a)
  end if
 case(1) ! natural orbitals
  deallocate(ev_a)
  call read_on_from_mkl(mklname, nif, 'a', eigen_e_a)
  call calc_dm_using_mo_and_on(nbf, nif, alpha_coeff, eigen_e_a, tot_dm)
 case(2) ! natural spin orbitals
  deallocate(ev_a, ev_b)
  allocate(eigen_e_b(nif))
  call read_on_from_mkl(mklname, nif, 'a', eigen_e_a)
  call read_on_from_mkl(mklname, nif, 'b', eigen_e_b)
  call calc_dm_using_mo_and_on(nbf, nif, alpha_coeff, eigen_e_a, tot_dm)
  allocate(spin_dm(nbf,nbf))
  call calc_dm_using_mo_and_on(nbf, nif, beta_coeff, eigen_e_b, spin_dm)
  allocate(coeff(nbf,nbf), source=tot_dm)
  tot_dm = coeff + spin_dm
  spin_dm = coeff - spin_dm
  deallocate(coeff)
 case default
  write(6,'(/,A)') 'ERROR in subroutine mkl2fch: no_type out of range!'
  write(6,'(A,I0)') 'no_type=', no_type
  stop
 end select

 call write_fch(fchname)
 deallocate(ielem, coor, iatom_type, shell_type, prim_per_shell, shell2atom_map,&
            prim_exp, contr_coeff, eigen_e_a, alpha_coeff, tot_dm)
 if(allocated(contr_coeff_sp)) deallocate(contr_coeff_sp)
 if(allocated(spin_dm)) deallocate(spin_dm)
 if(is_uhf) deallocate(eigen_e_b, beta_coeff)
end subroutine mkl2fch_direct

! check na, nb with those calculated from $OCC_ALPHA(and $OCC_BETA) in .mkl file
subroutine check_na_nb_ecp_in_mkl(mklname, uhf, nif, ne, na, nb)
 use mkl_content, only: read_on_from_mkl
 implicit none
 integer :: i, ne1, na1, nb1
 integer, intent(in) :: nif
 integer, intent(inout) :: ne, na, nb
 real(kind=8), parameter :: diff = 1d-4
 real(kind=8), allocatable :: on_a(:), on_b(:)
 character(len=240), intent(in) :: mklname
 logical, intent(in) :: uhf

 if(uhf) then ! UHF
  allocate(on_a(nif), on_b(nif))
  call read_on_from_mkl(mklname, nif, 'a', on_a)
  call read_on_from_mkl(mklname, nif, 'b', on_b)
  na1 = INT(SUM(on_a))
  nb1 = INT(SUM(on_b))
  deallocate(on_a, on_b)
 else         ! R(O)HF
  allocate(on_a(nif))
  call read_on_from_mkl(mklname, nif, 'a', on_a)
  na1 = 0; nb1 = 0
  do i = 1, nif, 1
   if(DABS(on_a(i) - 2d0) < diff) then
    na1 = na1 + 1
    nb1 = nb1 + 1
   else if(DABS(on_a(i) - 1d0) < diff) then
    na1 = na1 + 1
   else if(on_a(i) < diff) then
    exit
   end if
  end do ! for i
  deallocate(on_a)
 end if

 ne1 = na1 + nb1

 if((na==na1 .and. nb/=nb1) .or. (na/=na1 .and. nb==nb1) .or. na1>na .or. nb1>nb) then
  write(6,'(/,A)') 'ERROR in subroutine check_na_nb_ecp_in_mkl: internal error.'
  write(6,'(A,4I4)') 'na, nb, na1, nb1=', na, nb, na1, nb1 
  stop
 end if

 if(na1<na .and. nb1<nb) then
  write(6,'(/,A)') 'Warning from subroutine check_na_nb_ecp_in_mkl: ECP/PP dete&
                   &cted.'
  write(6,'(A)') 'NOTE: the .mkl file does not contain ECP/PP information. If y&
                 &ou use ECP/PP'
  write(6,'(A)') '(in ORCA .inp file), there would be no ECP in the generated i&
                 &nput file. You'
  write(6,'(A)') 'should manually add ECP data into the generated input file. I&
                 &f you are using'
  write(6,'(A,/)') 'an all-electron basis set, there is no problem.'
  ! update the number of electrons
  na = na1; nb = nb1; ne = ne1
 end if
end subroutine check_na_nb_ecp_in_mkl

