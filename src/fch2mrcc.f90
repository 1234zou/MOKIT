! written by jxzou at 20250327: transfer MOs from Gaussian to MRCC
! limitation: one element with different basis sets is not supported

program main
 use util_wrapper, only: formchk
 implicit none
 integer :: i
 character(len=240) :: fchname

 i = iargc()
 if(i /= 1) then
  write(6,'(/,A)') ' ERROR in subroutine fch2mrcc: wrong command line argument!'
  write(6,'(A,/)') ' Example: fch2mrcc water.fch'
  stop
 end if

 fchname = ' '
 call getarg(1, fchname)
 call require_file_exist(fchname)

 ! if .chk file provided, convert into .fch file automatically
 i = LEN_TRIM(fchname)
 if(fchname(i-3:i) == '.chk') then
  call formchk(fchname)
  fchname = fchname(1:i-3)//'fch'
 end if

 call fch2mrcc(fchname)
end program main

subroutine fch2mrcc(fchname)
 use fch_content
 implicit none
 integer :: i, j, k, m, length, icart, nif1
 integer :: n3pmark, n5dmark, n7fmark, n9gmark, n11hmark, n13imark, n6dmark, &
  n10fmark, n15gmark, n21hmark, n28imark
 integer, allocatable :: idx(:), p_mark(:), d_mark(:), f_mark(:), g_mark(:), &
  h_mark(:), i_mark(:)
 character(len=240), intent(in) :: fchname
 real(kind=8), allocatable :: coeff(:,:), coeff2(:,:), norm(:)
 logical :: uhf, sph, ecp, lin_dep

 call check_nobasistransform_in_fch(fchname)
 call check_nosymm_in_fch(fchname)

 sph = .true.
 call find_icart_in_fch(fchname, .false., icart)
 if(icart == 2) sph = .false.

 call find_irel_in_fch(fchname, irel)
 if(irel /= -1) then
  write(6,'(/,A)') REPEAT('-',79)
  write(6,'(A)') 'Warning from subroutine fch2mrcc: relativistic Hamiltonian de&
                 &tected.'
  write(6,'(A)') 'MRCC itself does not support relativistic Hamiltonian. It mus&
                 &t be used with'
  write(6,'(A)') 'CFOUR or Dirac. Anyway, the utility will continue...'
  write(6,'(A)') REPEAT('-',79)
 end if

 uhf = .false.; ecp = .false.; lin_dep = .false.
 call check_uhf_in_fch(fchname, uhf) ! determine whether UHF
 call read_fch(fchname, uhf)
 if(LenNCZ > 0) ecp = .true.
 if(nbf > nif) lin_dep = .true.

 ! generate the file MINP
 call prt_mrcc_minp(natom, charge, mult, uhf, sph, ecp, lin_dep, elem, coor, &
                    LPSkip)
 deallocate(coor)

 ! generate files GENBAS and ECPDATA(if needed)
 call prt_cfour_genbas(ecp, .true.)

 if(uhf) then ! UHF
  allocate(coeff(nbf,2*nif))
  coeff(:,1:nif) = alpha_coeff
  coeff(:,nif+1:2*nif) = beta_coeff
  deallocate(alpha_coeff, beta_coeff)
  nif1 = 2*nif
 else         ! R(O) HF
  allocate(coeff(nbf,nif), source=alpha_coeff)
  deallocate(alpha_coeff)
  nif1 = nif
 end if

 ! enlarge arrays shell_type and shell2atom_map, using d_mark as a tmp array
 k = ncontr
 allocate(d_mark(k), source=shell_type)
 deallocate(shell_type)
 allocate(shell_type(2*k), source=0)
 shell_type(1:k) = d_mark

 d_mark = shell2atom_map
 deallocate(shell2atom_map)
 allocate(shell2atom_map(2*k), source=0)
 shell2atom_map(1:k) = d_mark
 deallocate(d_mark)

! first we adjust the basis functions in each MO according to the Shell to atom map
! this is to ensure that D comes after L functions
 ! split the 'L' into 'S' and 'P'
 call split_L_func(k, shell_type, shell2atom_map, length)
 allocate(idx(nbf))
 forall(i = 1:nbf) idx(i) = i

 ! sort the shell_type, shell_to_atom_map by ascending order
 ! MOs will be adjusted accordingly
 call sort_shell_and_mo_idx(length, shell_type, shell2atom_map, nbf, idx)
 deallocate(shell2atom_map)
! adjust done

 ! record the indices of d, f, g and h functions
 k = length  ! update k
 allocate(p_mark(k), d_mark(k), f_mark(k), g_mark(k), h_mark(k), i_mark(k))

 ! adjust the order of 3p functions
 call read_3pmark_from_shltyp(k, shell_type, n3pmark, p_mark)
 do i = 1, n3pmark, 1
  j = p_mark(i)
  m = idx(j+2); idx(j+2) = idx(j); idx(j) = m
 end do ! for i
 deallocate(p_mark)

 if(sph) then
  call read_mark_from_shltyp_sph(k, shell_type, n5dmark, n7fmark, n9gmark, &
                 n11hmark, n13imark, d_mark, f_mark, g_mark, h_mark, i_mark)
  ! adjust the order of 5d/7f/9g/11h functions
  call fch2inporb_permute_sph(n5dmark, n7fmark, n9gmark, n11hmark, k, d_mark, &
                              f_mark, g_mark, h_mark, nbf, idx)
  allocate(coeff2(nbf,nif1), source=coeff)
  do i = 1, nif1, 1
   do j = 1, nbf, 1
    coeff(j,i) = coeff2(idx(j),i)
   end do ! for j
  end do ! for i
  deallocate(coeff2)
 else
  allocate(norm(nbf), source=1d0)
  call read_mark_from_shltyp_cart(k, shell_type, n6dmark, n10fmark, n15gmark, &
                    n21hmark, n28imark, d_mark, f_mark, g_mark, h_mark, i_mark)
  ! adjust the order of 6d/10f/15g/21h functions
  call fch2mrcc_permute_cart(n6dmark, n10fmark, n15gmark, n21hmark, k, d_mark,&
                             f_mark, g_mark, h_mark, nbf, idx, norm)
  allocate(coeff2(nbf,nif1), source=coeff)
  do i = 1, nif1, 1
   do j = 1, nbf, 1
    coeff(j,i) = coeff2(idx(j),i)*norm(j)
   end do ! for j
  end do ! for i
  deallocate(coeff2, norm)
 end if

 deallocate(shell_type, d_mark, f_mark, g_mark, h_mark, i_mark, idx)

 ! create/print MRCC orbital file MOCOEF
 if(uhf) then
  call prt_mrcc_mocoef(nbf, nif, coeff(:,1:nif), .false.)
  call prt_mrcc_mocoef(nbf, nif, coeff(:,nif+1:nif1), .true.)
 else
  call prt_mrcc_mocoef(nbf, nif, coeff, .false.)
 end if
 deallocate(coeff)
end subroutine fch2mrcc

! print MRCC MINP file
subroutine prt_mrcc_minp(natom, charge, mult, uhf, sph, ecp, lin_dep, elem, &
                         coor, LPSkip)
 use fch_content, only: period_nelem, period_elem, elem2nuc, nuc2elem
 use basis_data, only: ghost_elem
 implicit none
 integer :: i, k, max_ielem, fid
 integer, intent(in) :: natom, charge, mult
 integer, intent(in) :: LPSkip(natom)
 integer, allocatable :: ielem(:)
 real(kind=8), intent(in) :: coor(3,natom)
 character(len=2), intent(in) :: elem(natom)
 character(len=2), allocatable :: new_elem(:)
 logical, intent(in) :: uhf, sph, ecp, lin_dep
 logical :: has_ghost

 allocate(ielem(natom))
 max_ielem = 2

 do i = 1, natom, 1
  k = elem2nuc(elem(i))
  ielem(i) = k
  if(k > max_ielem) max_ielem = k
 end do ! for i

 if(ANY(ielem == 0)) then
  has_ghost = .true.
  ghost_elem = nuc2elem(max_ielem+1)
  allocate(new_elem(natom), source=elem)
  do i = 1, natom, 1
   if(ielem(i) == 0) new_elem(i) = ghost_elem
  end do ! for i
  write(6,'(/,A)') "Remark from subroutine prt_mrcc_minp: ghost atoms detected.&
                   & Using '"//ghost_elem//"' as"
  write(6,'(A)') 'the ghost atom symbols. The case that all ghost atoms have th&
                 &e same basis set'
  write(6,'(A,/)') 'is supported. Please DO NOT use different basis sets for di&
                   &fferent ghost atoms.'
 else
  has_ghost = .false.
  deallocate(ielem)
  allocate(new_elem(natom), source=elem)
 end if

 open(newunit=fid,file='MINP',status='replace')
 write(fid,'(A)') '# generated by fch2mrcc of MOKIT'
 write(fid,'(A)') 'mem=8GB'
 write(fid,'(A,I0)') 'charge=', charge
 write(fid,'(A,I0)') 'mult=', mult

 if(uhf) then
  write(fid,'(A)') 'calc=UHF'
 else if(mult > 1) then
  write(fid,'(A)') 'calc=ROHF'
 else
  write(fid,'(A)') 'calc=RHF'
 end if
 write(fid,'(A)') 'basis=PVTZ'

 if(ecp) then ! ECP/PP
  write(fid,'(A)') 'ecp=special'
  do i = 1, natom, 1
   if(LPSkip(i) == 0) then
    write(fid,'(A)') 'ECP-10-MDF'
   else
    write(fid,'(A)') 'none'
   end if
  end do ! for i
 else         ! no ECP/PP
  write(fid,'(A)') 'ecp=none'
 end if

 if(.not. sph) write(fid,'(A)') 'gauss=cart'
 if(lin_dep) write(fid,'(A)') 'ovltol=1e-6'
 write(fid,'(A)') 'scfiguess=mo'
 write(fid,'(A)') 'symm=0'
 ! The molecule will still be re-oriented when symmetry is turned off. Using
 ! zero background point charge keeps it fixed.
 write(fid,'(A)') 'qmmm=Amber'
 write(fid,'(/,A)') 'geom=xyz'
 write(fid,'(I0,/)') natom

 do i = 1, natom, 1
  write(fid,'(A2,3(1X,F18.8))') new_elem(i), coor(:,i)
 end do ! for i

 deallocate(new_elem)
 if(has_ghost) then
  write(fid,'(/,A)') 'ghost=serialno'
  do i = 1, natom, 1
   if(ielem(i) == 0) write(fid,'(I0,1X)',advance='no') i
  end do ! for i
  deallocate(ielem)
 end if

 write(fid,'(/,A)') 'pointcharges'
 write(fid,'(A)') '0'
 close(fid)
end subroutine prt_mrcc_minp

! read MRCC orbital file MOCOEF
subroutine read_mrcc_mocoef(nbf, nif, coeff)
 implicit none
 integer :: i, k, ii(2), fid
 integer, intent(in) :: nbf, nif
 integer(kind=8) :: k1
 real(kind=8), intent(out) :: coeff(nbf,nif)
 real(kind=8), allocatable :: mo(:)

 coeff = 0d0
 open(newunit=fid,file='MOCOEF',form='unformatted',status='old',position='rewind')
 read(fid) k1
 allocate(mo(2_8*k1))

 rewind(fid)
 read(fid) k1, mo
 close(fid)
 k = INT(k1, kind=4)

 do i = 1, k, 1
  ii = TRANSFER(mo(2*i), 1_4, 2)
  coeff(ii(1),ii(2)) = mo(2*i-1)
 end do ! for i

 deallocate(mo)
end subroutine read_mrcc_mocoef

! create/print MRCC orbital file MOCOEF
subroutine prt_mrcc_mocoef(nbf, nif, coeff, append)
 implicit none
 integer :: i, j, k, ncoeff, fid
 integer, intent(in) :: nbf, nif
 real(kind=8), intent(in) :: coeff(nbf,nif)
 real(kind=8), allocatable :: mo(:)
 logical, intent(in) :: append

 ncoeff = nbf*nif
 allocate(mo(2*ncoeff))
 k = 0

 do i = 1, nif, 1
  do j = 1, nbf, 1
   k = k + 1
   mo(2*k-1) = coeff(j,i)
   mo(2*k) = TRANSFER([j,i], 0d0)
  end do ! for j
 end do ! for i

 if(append) then
  open(newunit=fid,file='MOCOEF',form='unformatted',status='old',position='append')
 else
  open(newunit=fid,file='MOCOEF',form='unformatted',status='replace')
 end if

 write(fid) INT(ncoeff,kind=8), mo
 close(fid)
 deallocate(mo)
end subroutine prt_mrcc_mocoef

