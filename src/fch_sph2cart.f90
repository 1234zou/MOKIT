! written by jxzou at 20250923: convert a .fch(k) file with spherical harmonic
! functions into a new .fch with pure Cartesian functions. Sometimes this is
! called 5D to 6D, but note that higher angular momenta like F,G,H are also
! converted.

!TODO: support ECP/PP
!TODO: expand AO overlap integrals from sph -> cart, instead of calculating

program main
 use util_wrapper, only: formchk
 implicit none
 integer :: i, narg, RENAME
 character(len=240) :: sph_fch, cart_fch, def_fch

 narg = iargc()
 if(narg<1 .or. narg>2) then
  write(6,'(/,A)')' ERROR in program fch_sph2cart: wrong command line argument!'
  write(6,'(A)')  ' Example 1 (generate h2o_c.fch): fch_sph2cart h2o.fch'
  write(6,'(A,/)')' Example 2 (specify filename)  : fch_sph2cart h2o.fch h2o_ne&
                   &w.fch'
  stop
 end if

 sph_fch = ' '
 call getarg(1, sph_fch)
 call require_file_exist(sph_fch)

 ! if .chk file provided, convert into .fch file automatically
 i = LEN_TRIM(sph_fch)
 if(sph_fch(i-3:i) == '.chk') then
  call formchk(sph_fch)
  sph_fch = sph_fch(1:i-3)//'fch'
 end if

 call find_specified_suffix(sph_fch, '.fch', i)
 def_fch = sph_fch(1:i-1)//'_c.fch'

 call fch_sph2cart(sph_fch)

 if(narg == 2) then
  call getarg(2, cart_fch)
  if(TRIM(def_fch) /= TRIM(cart_fch)) then
   i = RENAME(TRIM(def_fch), TRIM(cart_fch))
  end if
 end if
end program main

subroutine fch_sph2cart(fchname)
 use fch_content
 implicit none
 integer :: i, nif0, nbf1, nif1, icart
 real(kind=8), allocatable :: ev(:), sph_coeff(:,:), cart_coeff(:,:), &
  cart_ovlp(:,:), new_mo(:,:)
 character(len=240) :: cart_fch
 character(len=240), intent(in) :: fchname
 logical :: uhf

 call find_specified_suffix(fchname, '.fch', i)
 cart_fch = fchname(1:i-1)//'_c.fch'
 call find_icart_in_fch(fchname, .false., icart)

 write(6,'(A)') REPEAT('-',79)

 if(icart==0 .or. icart==2) then
  if(icart == 0) then
   write(6,'(A)') 'Remark from subroutine fch_sph2cart: angular momenta are low&
                  &er than D. No'
   write(6,'(A)') 'conversion is needed. Now copy to '//TRIM(cart_fch)//' ...'
  else
   write(6,'(A)') 'Remark from subroutine fch_sph2cart: this is already a .fc&
                  &h(k) file with pure'
   write(6,'(A)') 'Cartesian functions. No conversion is needed. Now copy to '&
                  //TRIM(cart_fch)//' ...'
  end if
  write(6,'(A)') REPEAT('-',79)
  call sys_copy_file(fchname, cart_fch, .false.)
  return
 end if

 write(6,'(A)') 'Spherical harmonic functions detected. Now convert to pure Car&
                &tesian functions...'
 write(6,'(A)') REPEAT('-',79)
#ifdef _WIN32
 write(6,'(/,A)') 'ERROR in subroutine fch_sph2cart: this utility cannot be use&
                  &d on Windows'
 write(6,'(A)') 'since it needs to call PySCF (and PySCF cannot be used on Wind&
                &ows currently).'
 write(6,'(A)') 'Please use it on Linux/MacOS.'
 stop
#endif

 call check_uhf_in_fch(fchname, uhf)
 call read_fch(fchname, uhf)
 is_uhf = uhf

 if(uhf) then ! UHF
  nif1 = 2*nif
  allocate(sph_coeff(nbf,nif1))
  sph_coeff(:,1:nif) = alpha_coeff
  sph_coeff(:,nif+1:nif1) = beta_coeff
  deallocate(alpha_coeff, beta_coeff)
 else            ! R(O)HF
  nif1 = nif
  allocate(sph_coeff(nbf,nif1), source=alpha_coeff)
  deallocate(alpha_coeff)
 end if

 nbf1 = nbf + COUNT(shell_type==-2) + 3*COUNT(shell_type==-3) + &
         & 6*COUNT(shell_type==-4) + 10*COUNT(shell_type==-5) + &
           15*COUNT(shell_type==-6)
 ! [6D,10F,15G,21H,28I] - [5D,7F,9G,11H,13I] = [1,3,6,10,15]

! allocate(sph_ovlp(nbf,nbf), cart_ovlp(nbf1,nbf1))
! call get_gau_ao_ovlp_from_pyscf(fchname, nbf, sph_ovlp)
! call ovlp_sph2cart(ncontr, shell_type, nbf, nbf1, sph_ovlp, cart_ovlp)
 allocate(cart_coeff(nbf1,nif1))
 call mo_sph2cart(ncontr, shell_type, nbf, nbf1, nif1, sph_coeff, cart_coeff)
 deallocate(sph_coeff)

 ! The array dimensions of MO coefficients are (nbf1,nif2), where nif2 <= nbf1.
 ! But we do not know whether there is basis set linear dependency using pure
 ! Cartesian functions, currently we might as well set nif2 = nbf1. And set
 ! mo_coeff(:,nif+1:nbf1) to all zero. In the future, we need to compute the AO
 ! integral matrix at pure Cartesian basis sets, and see if there is basis set
 ! linear dependency.
 if(uhf) then
  allocate(alpha_coeff(nbf1,nbf1), source=0d0)
  alpha_coeff(:,1:nif) = cart_coeff(:,1:nif)
  allocate(beta_coeff(nbf1,nbf1), source=0d0)
  beta_coeff(:,1:nif) = cart_coeff(:,nif+1:nif1)

  allocate(ev(nif), source=eigen_e_b)
  deallocate(eigen_e_b)
  allocate(eigen_e_b(nbf1), source=0d0)
  eigen_e_b(1:nif) = ev
  deallocate(ev)

  if(allocated(spin_dm)) then
   deallocate(spin_dm)
   allocate(spin_dm(nbf1,nbf1), source=0d0)
  end if
 else
  allocate(alpha_coeff(nbf1,nbf1), source=0d0)
  alpha_coeff(:,1:nif) = cart_coeff
 end if

 deallocate(cart_coeff)
 allocate(ev(nif), source=eigen_e_a)
 deallocate(eigen_e_a)
 allocate(eigen_e_a(nbf1), source=0d0)
 eigen_e_a(1:nif) = ev
 deallocate(ev)

 if(allocated(tot_dm)) then
  deallocate(tot_dm)
  allocate(tot_dm(nbf1,nbf1), source=0d0)
 end if

 where (shell_type < -1)
  shell_type = -shell_type
 end where

 nbf = nbf1; nif0 = nif; nif = nbf1
 call write_fch(cart_fch)

 allocate(cart_ovlp(nbf1,nbf1))
 call get_gau_ao_ovlp_from_pyscf(cart_fch, nbf1, cart_ovlp)
 call get_nmo_from_ao_ovlp(nbf1, cart_ovlp, i)
 if(i < nbf1) then
  write(6,'(/,A)') 'ERROR in subroutine fch_sph2cart: basis set linear dependen&
                   &cy detected'
  write(6,'(A)') 'for current pure Cartesian basis. Not supported yet.'
  stop
 end if

 allocate(new_mo(nbf1,nbf1))
 call construct_vir(nbf1, nbf1, nif0+1, alpha_coeff, cart_ovlp, new_mo)

 if(uhf) then
  call write_mo_into_fch(cart_fch, nbf1, nbf1, 'a', new_mo)
  call construct_vir(nbf1, nbf1, nif0+1, beta_coeff, cart_ovlp, new_mo)
  deallocate(cart_ovlp)
  call write_mo_into_fch(cart_fch, nbf1, nbf1, 'b', new_mo)
 else
  deallocate(cart_ovlp)
  call write_mo_into_fch(cart_fch, nbf1, nbf1, 'a', new_mo)
 end if

 deallocate(new_mo)
 call free_arrays_in_fch_content()
end subroutine fch_sph2cart

