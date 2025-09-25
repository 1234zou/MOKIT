! written by jxzou at 20220816: transfer MOs from Gaussian -> Q-Chem
! Thanks to wsr for the previous version of fch2qchem (https://gitlab.com/jeanwsr/mokit)

! Current limitations:
! 1) not supported for GHF
! 2) not supported for different basis sets for the same element

program main
 use util_wrapper, only: formchk
 implicit none
 integer :: i, npair
 character(len=8) :: str
 character(len=240) :: fchname
 character(len=61), parameter :: error_warn = ' ERROR in subroutine fch2qchem:&
                                              & wrong command line arguments!'
 logical :: sfcis, sasfcis, sftd, sasf

 str = ' '; fchname = ' '; npair = 0
 sfcis = .false.; sasfcis = .false.; sftd = .false.; sasf = .false.

 i = iargc()
 if(i<1 .or. i>3) then
  write(6,'(/,A)') error_warn
  write(6,'(A)')  ' Example 1 (HF/DFT)   : fch2qchem h2o.fch'
  write(6,'(A)')  ' Example 2 (GVB)      : fch2qchem h2o.fch -gvb 2 (nopen is auto-detected)'
  write(6,'(A)')  ' Example 3 (SF-CIS)   : fch2qchem high_spin.fch -sfcis'
  write(6,'(A)')  ' Example 4 (SA-SF-CIS): fch2qchem high_spin.fch -sasfcis'
  write(6,'(A)')  ' Example 5 (SF-TDDFT) : fch2qchem high_spin.fch -sf'
  write(6,'(A,/)')' Example 6 (SA-SF-DFT): fch2qchem high_spin.fch -sasf'
  stop
 end if

 call getarg(1, fchname)
 call require_file_exist(fchname)

 if(i > 1) then
  call getarg(2, str)
  select case(TRIM(str))
  case('-gvb')
   if(i == 2) then
    write(6,'(/,A)') error_warn
    write(6,'(A)') 'You forget to specify the number of GVB pairs.'
    stop
   end if
   call getarg(3, str)
   read(str,*) npair
   if(npair < 0) then
    write(6,'(/,A)') error_warn
    write(6,'(A)') 'The 3rd argument npair should be >=0.'
    stop
   end if
  case('-sfcis')
   if(i == 3) then
    write(6,'(/,A)') error_warn
    write(6,'(A)') "Only two arguments are allowed when '-sfcis' is specified."
    stop
   end if
   sfcis = .true.
  case('-sasfcis')
   if(i == 3) then
    write(6,'(/,A)') error_warn
    write(6,'(A)') "Only two arguments are allowed when '-sasfcis' is specified."
    stop
   end if
   sasfcis = .true.
  case('-sf')
   if(i == 3) then
    write(6,'(/,A)') error_warn
    write(6,'(A)') "Only two arguments are allowed when '-sf' is specified."
    stop
   end if
   sftd = .true.
  case('-sasf')
   if(i == 3) then
    write(6,'(/,A)') error_warn
    write(6,'(A)') "Only two arguments are allowed when '-sasf' is specified."
    stop
   end if
   sasf = .true.
  case default
   write(6,'(/,A)') error_warn
   write(6,'(A)') "The 2nd argument can only be one of '-gvb', '-sasfcis', '-sf&
                  &cis', '-sf', '-sasf'."
   stop
  end select
 end if

 ! if .chk file provided, convert into .fch file automatically
 i = LEN_TRIM(fchname)
 if(fchname(i-3:i) == '.chk') then
  call formchk(fchname)
  fchname = fchname(1:i-3)//'fch'
 end if

 call fch2qchem(fchname, npair, sfcis, sasfcis, sftd, sasf)
end program main

subroutine fch2qchem(fchname, npair, sfcis, sasfcis, sftd, sasf)
 use fch_content
 implicit none
 integer :: i, j, k, m, n, n1, n2, nif1, icart, fid, purecart(4), SYSTEM
 integer, intent(in) :: npair
 integer, allocatable :: idx(:)
 character(len=1) :: str = ' '
 character(len=2) :: str2 = '  '
 character(len=1), parameter :: am_type(-1:6) = ['L','S','P','D','F','G','H','I']
 character(len=1), parameter :: am_type1(0:6) = ['s','p','d','f','g','h','i']
 character(len=240) :: proname, inpname, dirname
 character(len=240), intent(in) :: fchname
 real(kind=8), allocatable :: coeff0(:,:), coeff(:,:)
 logical :: spin_flip, uhf, sph, has_sp, ecp, so_ecp
 logical, intent(in) :: sfcis, sasfcis, sftd, sasf

 spin_flip = (sfcis .or. sasfcis .or. sftd .or. sasf)
 if(npair>0 .and. spin_flip) then
  write(6,'(/,A)') 'ERROR in subroutine fch2qchem: npair>0 is not allowed when &
                   &any SF-type method'
  write(6,'(A)') 'is activated. Please check your input arguments.'
  stop
 end if

 call find_specified_suffix(fchname, '.fch', i)
 proname = fchname(1:i-1)
 inpname = fchname(1:i-1)//'.in'
 call check_nobasistransform_in_fch(fchname)
 call check_nosymm_in_fch(fchname)
 sph = .true.
 call find_icart_in_fch(fchname, .false., icart)
 if(icart == 2) sph = .false.

 uhf = .false.; has_sp = .false.; ecp = .false.; so_ecp = .false.

 call check_uhf_in_fch(fchname, uhf) ! determine whether UHF
 if((sasfcis .or. sasf) .and. uhf) then
  write(6,'(/,A)') 'ERROR in subroutine fch2qchem: SA-SF-DFT must be based on a&
                   &n ROHF/ROKS reference.'
  write(6,'(A)') 'It seems that you provide a UHF/UKS .fch(k) file.'
  stop
 end if
 if(sftd .and. (.not.uhf)) then
  write(6,'(/,A)') 'ERROR in subroutine fch2qchem: SF-TDDFT must be based on a&
                   & UHF/UKS reference'
  write(6,'(A)') 'in Q-Chem. It seems that you provide an ROHF/ROKS .fch(k) file.'
  stop
 end if

 call read_fch(fchname, uhf)
 if(spin_flip .and. mult<3) then
  write(6,'(/,A)') 'ERROR in subroutine fch2qchem: SF-type methods must be base&
                   &d on a high-spin'
  write(6,'(A)') 'reference wave function, where the spin multiplicity should b&
                 &e >=3.'
  stop
 end if

 purecart = 1
 if(ANY(shell_type == 2)) purecart(4) = 2 ! 6D
 if(ANY(shell_type == 3)) purecart(3) = 2 ! 10F
 if(ANY(shell_type == 4)) purecart(2) = 2 ! 15G
 if(ANY(shell_type == 5)) purecart(1) = 2 ! 21H
 if(LenNCZ > 0) ecp = .true.

! Firstly, generate the input file (.in)
 open(newunit=fid,file=TRIM(inpname),status='replace')
 write(fid,'(A)') '$comment'
 write(fid,'(A)') ' file generated by fch2qchem utility of MOKIT'
 write(fid,'(A,/)') '$end'

 write(fid,'(A)') '$molecule'
 write(fid,'(I0,1X,I0)') charge, mult
 do i = 1, natom, 1
  write(fid,'(A2,1X,3(1X,F18.8))') elem(i), coor(:,i)
 end do ! for i
 deallocate(coor)
 write(fid,'(A,/)') '$end'

 write(fid,'(A)') '$rem'
 write(fid,'(A)') 'method hf'
 if(uhf) then ! UHF
  write(fid,'(A)') 'unrestricted true'
 else         ! R(O)HF
  write(fid,'(A)') 'unrestricted false'
 end if
 if(ecp) write(fid,'(A)') 'ecp gen'
 write(fid,'(A)') 'basis gen'
 write(fid,'(A)') 'scf_guess read'
 write(fid,'(A)') 'scf_convergence 8'
 write(fid,'(A)') 'thresh 12'
 write(fid,'(A,1X,4I0)') 'purecart', (purecart(i),i=1,4)
 write(fid,'(A)') 'sym_ignore true'
 if(npair > 0) then
  ! This is used for projection-equation PP:
  write(fid,'(A)') 'correlation pp'
  write(fid,'(A,I0)') 'gvb_n_pairs ',npair
  write(fid,'(A)') 'gvb_restart true'
  if(nopen > 0) write(fid,'(A,I0)') 'gvb_do_rohf ',npair
  ! This is used for variational PP (i.e. GVB):
  !write(fid,'(A)') 'gen_scfman false'
  !write(fid,'(A)') 'mp2_restart_no_scf true'
  !write(fid,'(A)') 'scf_algorithm diis'
  !write(fid,'(A)') 'correlation ccvb'
  !write(fid,'(A)') 'ccvb_guess 2'
  !write(fid,'(A)') 'ccvb_method 4'
  !write(fid,'(A,I0)') 'gvb_n_pairs ', npair
  !write(fid,'(A)') 'gvb_restart false'
 end if

 if(sfcis) then
  write(fid,'(A)') 'correlation ci'
  write(fid,'(A)') 'eom_corr cis'
  write(fid,'(A)') 'ccman2 false'
  write(fid,'(A)') 'sf_states 5'
 end if

 if(sftd .or. sasfcis .or. sasf) then
  if(.not. sasfcis) write(fid,'(A)') 'exchange bhhlyp'
  write(fid,'(A)') 'cis_n_roots 5'
  if(sftd) then
   write(fid,'(A)') 'spin_flip true'
  else
   write(fid,'(A)') 'sasf_rpa true'
  end if
  write(fid,'(A)') 'xc_grid 000099000590'
 end if

 write(fid,'(A)') 'gui = 2' ! generate fchk
 write(fid,'(A)') 'mem_total 4000'
 write(fid,'(A,/)') '$end'

 ! print basis sets into the .in file
 write(fid,'(A)') '$basis'
 write(fid,'(A,1X,A)') elem(1), '0'
 k = 0
 do i = 1, ncontr, 1
  m = shell2atom_map(i)
  if(m > 1) then
   if(shell2atom_map(i-1) == m-1) then
    write(fid,'(A4,/,A,1X,A)') '****', elem(m), '0'
   end if
  end if

  m = shell_type(i); n = prim_per_shell(i)
  if(m < -1) m = -m
  if(m == -1) then
   str2 = 'SP'
  else
   str2 = am_type(m)//' '
  end if
  write(fid,'(A2,1X,I2,3X,A)') str2, n, '1.00'

  has_sp = .false.
  if(allocated(contr_coeff_sp)) then
   if(ANY(contr_coeff_sp(k+1:k+n) > 1d-6)) has_sp = .true.
  end if

  if(has_sp) then
   do j = k+1, k+n, 1
    write(fid,'(3(2X,ES15.8))') prim_exp(j), contr_coeff(j), contr_coeff_sp(j)
   end do ! for j
  else ! no SP in this paragraph
   do j = k+1, k+n, 1
    write(fid,'(2(2X,ES15.8))') prim_exp(j), contr_coeff(j)
   end do ! for j
  end if

  k = k + n
 end do ! for i

 write(fid,'(A4,/,A)') '****','$end'
 deallocate(ielem, prim_per_shell, prim_exp, contr_coeff, shell2atom_map)
 if(allocated(contr_coeff_sp)) deallocate(contr_coeff_sp)

 if(ecp) then
  write(fid,'(/,A)') '$ecp'

  do i = 1, natom, 1
   if(LPSkip(i) /= 0) cycle
   write(fid,'(A)') elem(i)//'     0'
   write(fid,'(A,2X,I2,2X,I3)') TRIM(elem(i))//'-ECP', LMax(i), NINT(RNFroz(i))
   str = am_type1(LMax(i))

   do j = 1, 10, 1
    n1 = KFirst(i,j); n2 = KLast(i,j)
    if(n1 == 0) exit
    if(j == 1) then
     write(fid,'(A)') str//' potential'
    else
     write(fid,'(A)') am_type1(j-2)//'-'//str//' potential'
    end if
    write(fid,'(I2)') n2-n1+1
    do n = n1, n2, 1
     write(fid,'(I0,2(3X,ES15.8))') NLP(n), ZLP(n), CLP(n)
    end do ! for n
   end do ! for j

   write(fid,'(A)') '****'
  end do ! for i

  write(fid,'(A)') '$end'
  deallocate(KFirst, KLast, Lmax, LPSkip, NLP, RNFroz, CLP, CLP2, ZLP)
 end if

 close(fid)
 deallocate(elem)

! Secondly, permute MO coefficients and generate the orbital file 53.0
 if(uhf) then ! UHF
  allocate(coeff(nbf,2*nif))
  coeff(:,1:nif) = alpha_coeff
  coeff(:,nif+1:2*nif) = beta_coeff
  deallocate(alpha_coeff, beta_coeff)
  nif1 = 2*nif
 else         ! R(O) HF
  allocate(coeff(nbf,nif))
  coeff(:,:) = alpha_coeff
  deallocate(alpha_coeff)
  nif1 = nif
 end if

 allocate(coeff0(nbf,nif1), source=coeff)
 allocate(idx(nbf))
 call get_fch2qchem_permute_idx(sph, ncontr, shell_type, nbf, idx)

 forall(i=1:nif1, j=1:nbf) coeff(j,i) = coeff0(idx(j),i)
 deallocate(shell_type, coeff0, idx)

 allocate(alpha_coeff(nbf,nif), source=coeff(:,1:nif))
 if(uhf) allocate(beta_coeff(nbf,nif), source=coeff(:,nif+1:2*nif))
 deallocate(coeff)

 call create_dir(proname)
 !open(newunit=fid,file='53.0',form='binary')
 open(newunit=fid,file=TRIM(proname)//'/53.0',access='stream')

 write(unit=fid) alpha_coeff
 if(uhf) then
  write(unit=fid) beta_coeff
  deallocate(beta_coeff)
 else
  write(unit=fid) alpha_coeff
 end if

 write(unit=fid) eigen_e_a
 if(uhf) then
  write(unit=fid) eigen_e_b
  deallocate(eigen_e_b)
 else
  write(unit=fid) eigen_e_a
 end if

 deallocate(alpha_coeff, eigen_e_a)
 close(fid)
 ! For variational PP (i.e. GVB), there is no need for the file 169.0. This file
 ! is required when projection-equation PP is invoked.
 if(npair > 0) then
  call sys_copy_file(TRIM(proname)//'/53.0', TRIM(proname)//'/169.0', .false.)
  !if(nopen > 0) ??
 end if

 ! move the newly created directory into $QCSCRATCH/
 dirname = ' '
 call getenv('QCSCRATCH', dirname)

 if(LEN_TRIM(dirname) == 0) then
  write(6,'(/,A)') '$QCSCRATCH not found. '//TRIM(proname)//' put in the curren&
                   &t directory.'
  write(6,'(A)') 'You need to put the directory into $QCSCRATCH/ before running&
                 & qchem.'
 else
  call remove_dir(TRIM(dirname)//'/'//TRIM(proname))
  i = SYSTEM('mv '//TRIM(proname)//' '//TRIM(dirname)//'/')
  if(i == 0) then
   write(6,'(/,A)') '$QCSCRATCH found. Directory '//TRIM(proname)//' moved into &
                    &$QCSCRATCH/'
   write(6,'(A)') 'You can run:'
   write(6,'(A,/)') 'qchem '//TRIM(inpname)//' '//TRIM(proname)//'.out '//&
                     TRIM(proname)
  else
   write(6,'(/,A)') 'Warning in subroutine fch2qchem: failed to move directory&
                    & into '//TRIM(dirname)//'/'
  end if
 end if
end subroutine fch2qchem

