! written by jxzou at 20251026: convert Gaussian .fch(k) -> NWChem .nw

program main
 use util_wrapper, only: formchk
 implicit none
 integer :: i, k
 character(len=240) :: fchname

 i = iargc()
 if(i /= 1) then
  write(6,'(/,A)')' ERROR in program fch2nwc: wrong command line arguments!'
  write(6,'(A,/)')' Example: fch2nwc h2o.fch'
  stop
 end if

 fchname = ' '
 call getarg(1, fchname)
 call require_file_exist(fchname)

 ! if .chk file provided, convert into .fch file automatically
 k = LEN_TRIM(fchname)
 if(fchname(k-3:k) == '.chk') then
  call formchk(fchname)
  fchname = fchname(1:k-4)//'.fch'
 end if

 call fch2nwc(fchname)
end program main

! convert Gaussian .fch(k) -> NWChem .nw
subroutine fch2nwc(fchname)
 use fch_content
 implicit none
 integer :: i, k, icart
 real(kind=8), allocatable :: coeff(:,:)
 character(len=240) :: inpname
 character(len=240), intent(in) :: fchname
 logical :: uhf, ecp, sph

 call find_specified_suffix(fchname, '.fch', i)
 inpname = fchname(1:i-1)//'.nw'

 call check_nobasistransform_in_fch(fchname)
 call check_nosymm_in_fch(fchname)
 call find_irel_in_fch(fchname, irel)
 call check_uhf_in_fch(fchname, uhf) ! determine whether UHF
 call read_fch(fchname, uhf) ! read content in .fch(k) file
 is_uhf = uhf
 ecp = .false.
 if(LenNCZ > 0) ecp = .true.

 sph = .true.
 call find_icart_from_shell_type(.false., ncontr, shell_type, icart)
 if(icart == 2) sph = .false.

 call prt_nwchem_inp(inpname, sph)

 ! update MO coefficients
 if(uhf) then ! UHF
  k = 2*nif
  allocate(coeff(nbf,k))
  coeff(:,1:nif) = alpha_coeff
  coeff(:,nif+1:) = beta_coeff
 else         ! R(O)HF
  k = nif
  allocate(coeff(nbf,k), source=alpha_coeff)
 end if

 ! find F+3, G+3 and H+3 functions, multiply them by -1
! allocate(d_mark(ncontr), f_mark(ncontr), g_mark(ncontr), h_mark(ncontr), &
!          i_mark(ncontr))
! call read_mark_from_shltyp_sph(ncontr, shell_type, ndmark, nfmark, ngmark, &
!                      nhmark, nimark, d_mark, f_mark, g_mark, h_mark, i_mark)
! deallocate(d_mark)
! call update_mo_using_bas_mark(nbf, k, nfmark, ngmark, nhmark, nimark, ncontr,&
!                               f_mark, g_mark, h_mark, i_mark, coeff)
! deallocate(f_mark, g_mark, h_mark, i_mark)

 if(uhf) then ! UHF
  alpha_coeff = coeff(:,1:nif)
  beta_coeff = coeff(:,nif+1:)
 else         ! R(O)HF
  alpha_coeff = coeff
 end if
 deallocate(coeff)
 ! update MO coefficients done

 call free_arrays_in_fch_content()
end subroutine fch2nwc

! print/create .nw file
subroutine prt_nwchem_inp(inpname, sph)
 use fch_content
 implicit none
 integer :: i, j, k, m, n, fid
 integer, allocatable :: ntimes(:)
 character(len=1), parameter :: am_type(0:6) = ['S','P','D','F','G','H','I']
 character(len=10) :: str10
 character(len=240), intent(in) :: inpname
 logical, intent(in) :: sph

 ! find the number of times of each atom occurs
 allocate(ntimes(natom))
 call calc_ntimes(natom, elem, ntimes)

 open(newunit=fid,file=TRIM(inpname),status='replace')
 write(fid,'(A)') 'title "NWChem input file produced by fch2nwc of MOKIT"'
 write(fid,'(A)') 'memory total 4 gb'
 write(fid,'(A)') 'geometry'
 write(fid,'(A)') 'symmetry c1'
 do i = 1, natom, 1
  write(fid,'(A,I0,1X,3F17.9)') TRIM(elem(i)), ntimes(i), (coor(j,i),j=1,3)
 end do ! for i
 write(fid,'(A)') 'end'

 if(sph) then
  write(fid,'(A)') 'basis spherical'
 else
  write(fid,'(A)') 'basis cartesian'
 end if

 k = 0
 do i = 1, ncontr, 1
  j = shell2atom_map(i)
  write(str10,'(A,I0)') TRIM(elem(j)), ntimes(j)
  m = shell_type(i); n = prim_per_shell(i)

  if(m == -1) then ! 'L' or 'SP' in Pople-type basis sets
   write(fid,'(A)') TRIM(str10)//' S'
   do j = k+1, k+n, 1
    write(fid,'(2(2X,ES15.8))') prim_exp(j), contr_coeff(j)
   end do ! for j
   write(fid,'(A)') TRIM(str10)//' P'
   do j = k+1, k+n, 1
    write(fid,'(2(2X,ES15.8))') prim_exp(j), contr_coeff_sp(j)
   end do ! for j
  else             ! m<-1 or m=0,1
   write(fid,'(A)') TRIM(str10)//' '//am_type(IABS(m))
   do j = k+1, k+n, 1
    write(fid,'(2(2X,ES15.8))') prim_exp(j), contr_coeff(j)
   end do ! for j
  end if

  k = k + n
 end do ! for i

 deallocate(ntimes)
 write(fid,'(A)') 'end'
 write(fid,'(A,I0)') 'charge ', charge
 write(fid,'(A)') 'set lindep:tol 1e-6'

 write(fid,'(A)') 'scf'
 write(fid,'(A,I0)') ' nopen ', mult-1
 if(is_uhf) then
  write(fid,'(A)') ' uhf'
 else if(mult > 1) then
  write(fid,'(A)') ' rohf'
 else
  write(fid,'(A)') ' rhf'
 end if
 write(fid,'(A)') ' tol2e 1e-12'
 write(fid,'(A)') 'end'

 write(fid,'(A)') 'property'
 write(fid,'(A)') ' moldenfile'
 write(fid,'(A)') ' molden_norm nwchem'
 write(fid,'(A)') 'end'
 write(fid,'(A)') 'task scf property'
 close(fid)
end subroutine prt_nwchem_inp

