! written by jxzou at 20250916: transfer MOs from Gaussian to OpenQP
! TODO: ghost atoms case
! TODO: ECP/PP support
! Limitation: not support different basis set for atoms of the same element

program main
 use util_wrapper, only: formchk
 implicit none
 integer :: i, k, sf_type
 character(len=8) :: str8
 character(len=240) :: fchname

 i = iargc()
 if(.not. (i==1 .or. i==2)) then
  write(6,'(/,A)')' ERROR in program fch2openqp: wrong command line argument!'
  write(6,'(A)')  ' Example 1 (R(O)HF, UHF): fch2openqp h2o.fch'
  write(6,'(A)')  ' Example 2 (SF-CIS)     : fch2openqp O2_rohf.fch -sfcis'
  write(6,'(A)')  ' Example 3 (SF-TDDFT)   : fch2openqp O2_roks.fch -sf'
  write(6,'(A)')  ' Example 4 (MRSF-CIS)   : fch2openqp O2_rohf.fch -mrsfcis'
  write(6,'(A,/)')' Example 5 (MRSF-TDDFT) : fch2openqp O2_roks.fch -mrsf'
  stop
 end if

 fchname = ' '
 call getarg(1, fchname)
 call require_file_exist(fchname)

 ! if .chk file provided, convert into .fch file automatically
 k = LEN_TRIM(fchname)
 if(fchname(k-3:k) == '.chk') then
  call formchk(fchname)
  fchname = fchname(1:k-3)//'fch'
 end if

 sf_type = 0
 if(i > 1) then
  call getarg(2, str8)
  select case(TRIM(str8))
  case('-sfcis')
   sf_type = 1
  case('-sf')
   sf_type = 2
  case('-mrsfcis')
   sf_type = 3
  case('-mrsf')
   sf_type = 4
  case default
   write(6,'(/,A)') 'ERROR in program fch2openqp: the 2nd command line argument&
                    & can only be'
   write(6,'(A)') "'-sfcis', '-sf', '-mrsfcis, '-mrsf'. But got "//TRIM(str8)
   stop
  end select
 end if

 call fch2openqp(fchname, sf_type)
end program main

subroutine fch2openqp(fchname, sf_type)
 use fch_content
 implicit none
 integer :: i, j, k, length, nif1, hf_type, icart
 integer :: n6dmark, n10fmark, n15gmark, n21hmark, n28imark
 integer, allocatable :: idx(:), d_mark(:), f_mark(:), g_mark(:), h_mark(:), &
  i_mark(:)
 integer, intent(in) :: sf_type
 character(len=30) :: dftname
 character(len=240) :: proname, bas_json, cart_fch
 character(len=240), intent(in) :: fchname
 real(kind=8), allocatable :: coeff(:,:), coeff2(:,:), ev(:)
 logical :: uhf, lin_dep

 call find_specified_suffix(fchname, '.fch', i)
 proname = fchname(1:i-1)
 cart_fch = fchname(1:i-1)//'_c.fch'
 bas_json = fchname(1:i-1)//'_bas.json'
 call lower(bas_json) ! lower case is needed for this file

 call check_nobasistransform_in_fch(fchname)
 call check_nosymm_in_fch(fchname)

 call find_icart_in_fch(fchname, .false., icart)
 if(icart == 1) then
  write(6,'(/,A)') 'ERROR in subroutine fch2openqp: spherical harmonic function&
                    &s detected.'
  write(6,'(A)') 'OpenQP supports only pure Cartesian functions (6D 10F) curren&
                 &tly. There are'
  write(6,'(A)') 'two solutions:'
  write(6,'(/,A)') '(1) Project current MOs onto pure Cartesian basis then use &
                   &fch2openqp, e.g.'
  write(6,'(A)') '```'
  write(6,'(A)') 'fch_sph2cart '//TRIM(fchname)
  write(6,'(A)') 'fch2openqp '//TRIM(cart_fch)
  write(6,'(A)') '```'
  write(6,'(A)') 'Remember to add extra arguments -mrsf/-sf if you need them.'
  write(6,'(/,A)') '(2) Add keywords `6D 10F` into .gjf if you are using Gaussi&
                   &an to generate this'
  write(6,'(A)') '.fch(k) file, or set `mol.cart=True` in Python script if PySC&
                 &F and the `fchk()`'
  write(6,'(A,/)') 'module are used. Note that (1) is usually much faster than &
                   &(2).'
  stop
 end if

 uhf = .false.; lin_dep = .false.; dftname = ' '
 call check_uhf_in_fch(fchname, uhf) ! determine whether UHF
 call read_fch(fchname, uhf)
 if(nbf > nif) lin_dep = .true.

 hf_type = 0
 if(uhf) then
  hf_type = 3
 else if(mult > 1) then
  hf_type = 2
 else
  hf_type = 1
 end if

 if((sf_type==1 .or. sf_type==2) .and. mult<3) then
  write(6,'(/,A)') 'ERROR in subroutine fch2openqp: spin-flip methods requires &
                   &the spin multiplicity'
  write(6,'(A,I0)') 'to be >=3. But got mult=', mult
  stop
 end if

 if((sf_type==3 .or. sf_type==4) .and. mult/=3) then
  write(6,'(/,A)') 'ERROR in subroutine fch2openqp: MRSF methods requires the s&
                   &pin multiplicity'
  write(6,'(A,I0)') 'to be 3 currently. But got mult=', mult
  stop
 end if

 if(uhf .and. sf_type>0) then
  write(6,'(/,A)') 'ERROR in subroutine fch2openqp: UHF-based SF/MRSF methods a&
                   &re not supported.'
  write(6,'(A)') 'Please ROHF/ROKS as reference.'
 end if

 ! generate the OpenQP input file
 call prt_openqp_inp(proname, dftname, natom, charge, mult, sf_type, uhf, &
                     ielem, coor)

 ! generate the basis set file xxx_bas.json (ECP/PP not considered so far)
 call prt_bas_json(bas_json)

 ! copy alpha (and beta) MO coefficients into the array coeff
 if(uhf) then ! UHF
  allocate(coeff(nbf,2*nif), ev(2*nif))
  coeff(:,1:nif) = alpha_coeff
  coeff(:,nif+1:2*nif) = beta_coeff
  ev(1:nif) = eigen_e_a
  ev(nif+1:2*nif) = eigen_e_b
  nif1 = 2*nif
 else         ! R(O) HF
  allocate(coeff(nbf,nif), source=alpha_coeff)
  allocate(ev(nif), source=eigen_e_a)
  nif1 = nif
 end if

 ! enlarge arrays shell_type and shell2atom_map, using f_mark as a tmp array
 k = ncontr
 allocate(f_mark(k), source=shell_type)
 deallocate(shell_type)
 allocate(shell_type(2*k), source=0)
 shell_type(1:k) = f_mark

 f_mark = shell2atom_map
 deallocate(shell2atom_map)
 allocate(shell2atom_map(2*k), source=0)
 shell2atom_map(1:k) = f_mark
 deallocate(f_mark)

! first we adjust the basis functions in each MO according to the Shell to atom map
! this is to ensure that D comes after L functions
 ! split the 'L' into 'S' and 'P'
 call split_L_func(k, shell_type, shell2atom_map, length)
 allocate(idx(nbf))
 forall(i = 1:nbf) idx(i) = i

 ! sort the shell_type, shell_to_atom_map by ascending order
 ! MOs will be adjusted accordingly
 call sort_shell_and_mo_idx(length, shell_type, shell2atom_map, nbf, idx)
! adjust done

 ! record the indices of f, g and h functions
 k = length  ! update k
 allocate(d_mark(k), f_mark(k), g_mark(k), h_mark(k), i_mark(k))
 call read_mark_from_shltyp_cart(k, shell_type, n6dmark, n10fmark, n15gmark, &
                   n21hmark, n28imark, d_mark, f_mark, g_mark, h_mark, i_mark)
 deallocate(d_mark, i_mark)
 ! adjust the order of 10f/15g/21h functions
 call fch2openqp_permute_cart(n10fmark, n15gmark, n21hmark, k, f_mark, g_mark, &
                              h_mark, nbf, idx)
 deallocate(f_mark, g_mark, h_mark)
 allocate(coeff2(nbf,nif1), source=coeff)

 do i = 1, nif1, 1
  do j = 1, nbf, 1
   coeff(j,i) = coeff2(idx(j),i)
  end do ! for j
 end do ! for i
 deallocate(coeff2, idx)

 ! generate the wave function file xxx_wfn.json
 coor = coor/Bohr_const ! convert Anstrom to Bohr
 call prt_openqp_wfn_json(proname, hf_type, natom, nbf, nif1, na, nb, ielem, &
                          coor, ev, coeff)
 deallocate(ev, coeff)

 call free_arrays_in_fch_content()
end subroutine fch2openqp

! print/create an OpenQP input file
subroutine prt_openqp_inp(proname, dftname, natom, charge, mult, sf_type, uhf, &
                          nuc, coor)
 implicit none
 integer :: i, fid
 integer, intent(in) :: natom, charge, mult, sf_type
 integer, intent(in) :: nuc(natom)
 real(kind=8), intent(in) :: coor(3,natom)
 character(len=30) :: dftname1
 character(len=30), intent(in) :: dftname
 character(len=240) :: inpname, wfn_json, bas_json
 character(len=240), intent(in) :: proname
 logical, intent(in) :: uhf
 logical :: dft

 inpname = TRIM(proname)//'.inp'
 wfn_json = TRIM(proname)//'_wfn.json'
 bas_json = TRIM(proname)//'_bas.json'
 call lower(bas_json) ! lower case is needed for this file

 dft = .false.
 if(LEN_TRIM(dftname) > 0) then
  dft = .true.
  dftname1 = dftname
 end if
 if(sf_type==2 .or. sf_type==4) then
  dft = .true.
  dftname1 = 'bhhlyp'
 end if

 open(newunit=fid,file=TRIM(inpname),status='replace')
 write(fid,'(A,/,A)') '[input]', 'system='

 do i = 1, natom, 1
  write(fid,'(I4,3(1X,F17.8))') nuc(i), coor(:,i)
 end do ! for i
 write(fid,'(A,I0)') 'charge=', charge
 write(fid,'(A)') 'runtype=energy'
 if(dft) write(fid,'(A)') 'functional='//TRIM(dftname1)

 write(fid,'(A)') 'basis=file:'//TRIM(bas_json)
 if(sf_type > 0) then
  write(fid,'(A)') 'method=tdhf'
 else
  write(fid,'(A)') 'method=hf'
 end if

 write(fid,'(/,A)') '[guess]'
 write(fid,'(A)') 'type=json'
 write(fid,'(A)') 'file='//TRIM(wfn_json)
 write(fid,'(A)') 'save_mol=True'

 write(fid,'(/,A)') '[scf]'
 write(fid,'(A,I0)') 'multiplicity=', mult
 write(fid,'(A,I0)') 'conv=1e-7'
 write(fid,'(A)') 'maxit=50'
 write(fid,'(A)') 'save_molden=False'
 write(fid,'(A)',advance='no') 'type='
 if(uhf) then
  write(fid,'(A)') 'uhf'
 else if(mult > 1) then
  write(fid,'(A)') 'rohf'
 else
  write(fid,'(A)') 'rhf'
 end if

 if(sf_type > 0) then
  write(fid,'(/,A)') '[tdhf]'
  if(sf_type < 3) then
   write(fid,'(A)') 'type=sf'
  else
   write(fid,'(A)') 'type=mrsf'
  end if
  write(fid,'(A)') 'nstate=5'
 end if

 if(dft) then
  write(fid,'(/,A)') '[dftgrid]'
  !write(fid,'(A)') 'rad_type=log3'
  write(fid,'(A)') 'rad_npts=99'
  write(fid,'(A)') 'ang_npts=590'
 end if
 close(fid)
end subroutine prt_openqp_inp

! print the density matrix array into an OpenQP .json file
subroutine prt_dm2openqp_json(fid, nbf, dm, alpha)
 implicit none
 integer :: i, j
 integer, intent(in) :: fid, nbf
 real(kind=8), intent(in) :: dm(nbf,nbf)
 logical, intent(in) :: alpha

 if(alpha) then
  write(fid,'(2X,A)') '"OQP::DM_A": ['
 else
  write(fid,'(2X,A)') '"OQP::DM_B": ['
 end if

 do i = 1, nbf-1, 1
  do j = 1, i, 1
   write(fid,'(4X,ES15.8,A)') dm(j,i),','
  end do ! for j
 end do ! for i
 do j = 1, nbf-1, 1
  write(fid,'(4X,ES15.8,A)') dm(j,nbf),','
 end do ! for j
 write(fid,'(4X,ES15.8,A)') dm(nbf,nbf)

 write(fid,'(2X,A)') '],'
end subroutine prt_dm2openqp_json

! print alpha/beta MO eigenvalues (which are usually orbital energies) into an
! OpenQP .json file
subroutine prt_ev2openqp_json(fid, nif, ev, alpha)
 implicit none
 integer :: i
 integer, intent(in) :: fid, nif
 real(kind=8), intent(in) :: ev(nif)
 logical, intent(in) :: alpha

 if(alpha) then
  write(fid,'(2X,A)') '"OQP::E_MO_A": ['
 else
  write(fid,'(2X,A)') '"OQP::E_MO_B": ['
 end if

 do i = 1, nif-1, 1
  write(fid,'(4X,ES15.8,A)') ev(i),','
 end do ! for i
 write(fid,'(4X,ES15.8)') ev(nif)

 write(fid,'(2X,A)') '],'
end subroutine prt_ev2openqp_json

! print alpha/beta MO coefficients into an OpenQP .json file
subroutine prt_mo2openqp_json(fid, nbf, nif, mo, alpha)
 implicit none
 integer :: i, j
 integer, intent(in) :: fid, nbf, nif
 real(kind=8), intent(in) :: mo(nbf,nif)
 logical, intent(in) :: alpha

 if(alpha) then
  write(fid,'(2X,A)') '"OQP::VEC_MO_A": ['
 else
  write(fid,'(2X,A)') '"OQP::VEC_MO_B": ['
 end if

 do i = 1, nif, 1
  write(fid,'(4X,A)') '['
  do j = 1, nbf-1, 1
   write(fid,'(8X,ES15.8,A)') mo(j,i),','
  end do ! for j
  write(fid,'(8X,ES15.8)') mo(nbf,i)
  if(i == nif) then
   write(fid,'(4X,A)') ']'
  else
   write(fid,'(4X,A)') '],'
  end if
 end do ! for i

 write(fid,'(2X,A)') '],'
end subroutine prt_mo2openqp_json

! Print/create an OpenQP wave function file xxx_wfn.json.
! Note: 1) the input array mo must contain OpenQP-convention MO coefficients;
!       2) the input array coor must contain Cartesian coordinates in Bohr;
!       3) for UHF-type MOs (i.e. hf_type=3), nif is two times of No. MOs.
subroutine prt_openqp_wfn_json(proname, hf_type, natom, nbf, nif, na, nb, &
                               ielem, coor, ev, mo)
 implicit none
 integer :: i, j, fid
 integer, intent(in) :: hf_type, natom, nbf, nif, na, nb
 integer, intent(in) :: ielem(natom)
 real(kind=8), intent(in) :: coor(3,natom), ev(nif), mo(nbf,nif)
 real(kind=8), allocatable :: dm(:,:)
 character(len=240) :: wfn_json, bas_json
 character(len=240), intent(in) :: proname

 wfn_json = TRIM(proname)//'_wfn.json'
 bas_json = TRIM(proname)//'_bas.json'
 call lower(bas_json) ! lower case is needed for this file

 open(newunit=fid,file=TRIM(wfn_json),status='replace')
 write(fid,'(A)') '{'

 allocate(dm(nbf,nbf))
 call calc_cct(nbf, na, mo(:,1:na), dm)

 select case(hf_type)
 case(1) ! RHF
  dm = 2d0*dm
  call prt_dm2openqp_json(fid, nbf, dm, .true.)
  deallocate(dm)
  call prt_ev2openqp_json(fid, nif, ev, .true.)
  call prt_mo2openqp_json(fid, nbf, nif, mo, .true.)
 case(2) ! ROHF
  call prt_dm2openqp_json(fid, nbf, dm, .true.)
  call calc_cct(nbf, nb, mo(:,1:nb), dm)
  call prt_dm2openqp_json(fid, nbf, dm, .false.)
  deallocate(dm)
  call prt_ev2openqp_json(fid, nif, ev, .true.)
  call prt_ev2openqp_json(fid, nif, ev, .false.)
  call prt_mo2openqp_json(fid, nbf, nif, mo, .true.)
  call prt_mo2openqp_json(fid, nbf, nif, mo, .false.)
 case(3) ! UHF
  call prt_dm2openqp_json(fid, nbf, dm, .true.)
  call calc_cct(nbf, nb, mo(:,1:nb), dm)
  call prt_dm2openqp_json(fid, nbf, dm, .false.)
  deallocate(dm)
  i = nif/2
  call prt_ev2openqp_json(fid, i, ev(1:i), .true.)
  call prt_ev2openqp_json(fid, i, ev(i+1:nif), .false.)
  call prt_mo2openqp_json(fid, nbf, i, mo(:,1:i), .true.)
  call prt_mo2openqp_json(fid, nbf, i, mo(:,i+1:nif), .false.)
 case default
  write(6,'(/,A,I0)') 'ERROR in subroutine prt_openqp_wfn_json: invalid hf_type&
                      &=', hf_type
  stop
 end select

 write(fid,'(2X,A)') '"atoms": ['
 do i = 1, natom-1, 1
  write(fid,'(4X,I0,A)') ielem(i),','
 end do ! for i
 write(fid,'(4X,I0)') ielem(natom)
 write(fid,'(2X,A)') '],'

 write(fid,'(2X,A)') '"coord": ['
 do i = 1, natom-1, 1
  do j = 1, 3
   write(fid,'(4X,F18.8,A)') coor(j,i),','
  end do ! for j
 end do ! for i
 write(fid,'(4X,F18.8,A)') coor(1,natom),','
 write(fid,'(4X,F18.8,A)') coor(2,natom),','
 write(fid,'(4X,F18.8,A)') coor(3,natom)
 write(fid,'(2X,A)') '],'

 write(fid,'(2X,A)') '"energy": 0.0,'
 write(fid,'(2X,A)') '"td_energies": ['
 write(fid,'(4X,A)') '0.0'
 write(fid,'(2X,A)') '],'
 write(fid,'(2X,A)') '"grad": [],'
 write(fid,'(2X,A)') '"nac": [],'
 write(fid,'(2X,A)') '"soc": [],'
 write(fid,'(2X,A)') '"hess": [],'

 write(fid,'(2X,A)') '"json": {'
 select case(hf_type)
 case(1) ! RHF
  write(fid,'(4X,A)') '"scf_type": "rhf",'
 case(2) ! ROHF
  write(fid,'(4X,A)') '"scf_type": "rohf",'
 case(3) ! UHF
  write(fid,'(4X,A)') '"scf_type": "uhf",'
 end select
 write(fid,'(4X,A)') '"basis": "file:'//TRIM(bas_json)//'",'
 write(fid,'(4X,A)') '"library": ""'
 write(fid,'(2X,A)') '}'

 write(fid,'(A)') '}'
 close(fid)
end subroutine prt_openqp_wfn_json

