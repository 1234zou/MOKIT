! written by jxzou at 20230217

! Q-Chem .fch(k) -> Amesp (.aip, .amo)
subroutine qchem2amesp(fchname, aipname)
 use util_wrapper, only: fch2amo_wrap
 implicit none
 integer :: i
 character(len=240) :: std_fch
 character(len=240), intent(in) :: fchname, aipname
!f2py intent(in) :: fchname, aipname

 call find_specified_suffix(aipname, '.aip', i)
 call find_specified_suffix(fchname, '.', i)
 std_fch = fchname(1:i-1)//'_std.fch'
 call standardize_fch(fchname)
 call fch2amo_wrap(std_fch, aipname)
 call delete_file(TRIM(std_fch))
end subroutine qchem2amesp

! Q-Chem .fch(k) -> BDF (.inp, .scforb, .BAS)
subroutine qchem2bdf(fchname, inpname)
 implicit none
 integer :: i, SYSTEM, RENAME
 character(len=240) :: std_fch, std_inp, std_orb, orbname
 character(len=240), intent(in) :: fchname, inpname
!f2py intent(in) :: fchname, inpname

 call find_specified_suffix(inpname, '.inp', i)
 orbname = inpname(1:i-1)//'.scforb'

 call find_specified_suffix(fchname, '.', i)
 std_fch = fchname(1:i-1)//'_std.fch'
 std_inp = fchname(1:i-1)//'_std_bdf.inp'
 std_orb = fchname(1:i-1)//'_std_bdf.scforb'

 call standardize_fch(fchname)
 i = SYSTEM('fch2bdf '//TRIM(std_fch))
 call delete_file(TRIM(std_fch))
 i = RENAME(TRIM(std_inp), TRIM(inpname))
 i = RENAME(TRIM(std_orb), TRIM(orbname))
end subroutine qchem2bdf

! Q-Chem .fch(k) -> CFOUR (ZMAT, OLDMOS, GENBAS, ECPDATA)
subroutine qchem2cfour(fchname)
 use util_wrapper, only: fch2cfour_wrap
 implicit none
 integer :: i
 character(len=240) :: std_fch
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 call find_specified_suffix(fchname, '.', i)
 std_fch = fchname(1:i-1)//'_std.fch'
 call standardize_fch(fchname)
 call fch2cfour_wrap(std_fch)
 call delete_file(TRIM(std_fch))
end subroutine qchem2cfour

! Q-Chem .fch(k) -> Dalton (.dal, .mol)
subroutine qchem2dalton(fchname, dalname)
 use util_wrapper, only: fch2dal_wrap
 implicit none
 integer :: i
 character(len=240) :: std_fch, molname
 character(len=240), intent(in) :: fchname, dalname
!f2py intent(in) :: fchname, dalname

 call find_specified_suffix(dalname, '.dal', i)
 molname = dalname(1:i-1)//'.mol'

 call find_specified_suffix(fchname, '.', i)
 std_fch = fchname(1:i-1)//'_std.fch'

 call standardize_fch(fchname)
 call fch2dal_wrap(std_fch, dalname)
 call delete_file(TRIM(std_fch))
end subroutine qchem2dalton

! Q-Chem .fch(k) -> GAMESS (.inp)
subroutine qchem2gms(fchname, inpname)
 use util_wrapper, only: fch2inp_wrap
 implicit none
 integer :: i, RENAME
 character(len=240) :: std_fch, std_inp
 character(len=240), intent(in) :: fchname, inpname
!f2py intent(in) :: fchname, inpname

 call find_specified_suffix(fchname, '.', i)
 std_fch = fchname(1:i-1)//'_std.fch'
 std_inp = fchname(1:i-1)//'_std.inp'
 call standardize_fch(fchname)
 call fch2inp_wrap(std_fch, .false., 0, 0, .false.)
 call delete_file(TRIM(std_fch))
 i = RENAME(TRIM(std_inp), TRIM(inpname))
end subroutine qchem2gms

! Q-Chem .fch(k) -> (Open)Molcas (.inporb, .INPORB)
subroutine qchem2molcas(fchname, inpname)
 use util_wrapper, only: fch2inporb_wrap
 implicit none
 integer :: i
 character(len=240) :: std_fch
 character(len=240), intent(in) :: fchname, inpname
!f2py intent(in) :: fchname, inpname

 call find_specified_suffix(inpname, '.input', i)
 call find_specified_suffix(fchname, '.', i)
 std_fch = fchname(1:i-1)//'_std.fch'
 call standardize_fch(fchname)
 call fch2inporb_wrap(std_fch, .false., inpname)
 call delete_file(TRIM(std_fch))
end subroutine qchem2molcas

! Q-Chem .fch(k) -> Molpro (.com, .a, .b)
subroutine qchem2molpro(fchname, inpname)
 use util_wrapper, only: fch2com_wrap
 implicit none
 integer :: i
 character(len=240) :: std_fch
 character(len=240), intent(in) :: fchname, inpname
!f2py intent(in) :: fchname, inpname

 call find_specified_suffix(inpname, '.com',i)
 call find_specified_suffix(fchname, '.', i)
 std_fch = fchname(1:i-1)//'_std.fch'
 call standardize_fch(fchname)
 call fch2com_wrap(std_fch, inpname)
 call delete_file(TRIM(std_fch))
end subroutine qchem2molpro

! Q-Chem .fch(k) -> PSI4 (.inp, .A, .B)
subroutine qchem2psi(fchname, inpname)
 use util_wrapper, only: fch2psi_wrap
 implicit none
 integer :: i
 character(len=240) :: std_fch
 character(len=240), intent(in) :: fchname, inpname
!f2py intent(in) :: fchname, inpname

 call find_specified_suffix(inpname, '.inp', i)
 call find_specified_suffix(fchname, '.', i)
 std_fch = fchname(1:i-1)//'_std.fch'
 call standardize_fch(fchname)
 call fch2psi_wrap(std_fch, inpname)
 call delete_file(TRIM(std_fch))
end subroutine qchem2psi

! Q-Chem .fch(k) -> PySCF (.py, .fch)
subroutine qchem2pyscf(fchname, pyname)
 use util_wrapper, only: bas_fch2py_wrap
 implicit none
 integer :: i
 character(len=240) :: std_fch
 character(len=240), intent(in) :: fchname, pyname
!f2py intent(in) :: fchname, pyname

 call find_specified_suffix(pyname, '.py', i)
 call find_specified_suffix(fchname, '.', i)
 std_fch = fchname(1:i-1)//'_std.fch'
 call standardize_fch(fchname)
 call bas_fch2py_wrap(std_fch, .false., pyname)
end subroutine qchem2pyscf

! Q-Chem .fch(k) -> ORCA (.inp, .mkl, .gbw)
subroutine qchem2orca(fchname, inpname)
 use util_wrapper, only: fch2mkl_wrap, mkl2gbw
 implicit none
 integer :: i, RENAME
 character(len=240) :: std_fch, std_inp, mklname, gbwname
 character(len=240), intent(in) :: fchname, inpname
!f2py intent(in) :: fchname, inpname

 call find_specified_suffix(inpname, '.inp', i)
 mklname = inpname(1:i-1)//'.mkl'
 gbwname = inpname(1:i-1)//'.gbw'
 call find_specified_suffix(fchname, '.', i)
 std_fch = fchname(1:i-1)//'_std.fch'
 std_inp = fchname(1:i-1)//'_std_o.inp'

 call standardize_fch(fchname)
 call fch2mkl_wrap(std_fch, mklname)
 call delete_file(TRIM(std_fch))
 i = RENAME(TRIM(std_inp), TRIM(inpname))
 call mkl2gbw(mklname, gbwname)
end subroutine qchem2orca

! standardize a .fch file, usually for Q-Chem generated .fchk file
subroutine standardize_fch(fchname)
 use fch_content
 implicit none
 integer :: i, ne1, ne2
 character(len=240) :: new_fch
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname
 logical, external :: has_pople_sp

 i = INDEX(fchname, '.', back=.true.)
 if(i == 0) then
  write(6,'(/,A)') "ERROR in subroutine standardize_fch: no '.' symbol found in&
                   & filename "//TRIM(fchname)
  stop
 end if
 new_fch = fchname(1:i-1)//'_std.fch'

 call read_charge_and_mult_from_fch(fchname, charge, mult)
 call read_na_and_nb_from_fch(fchname, na, nb)
 call read_natom_from_fch(fchname, natom)
 allocate(ielem(natom), coor(3,natom))
 call read_nuc_from_fch(fchname, natom, ielem)
 call read_coor_from_fch(fchname, natom, coor)

 call read_ncontr_from_fch(fchname, ncontr)
 allocate(shell_type(ncontr), prim_per_shell(ncontr), shell2atom_map(ncontr))
 call read_shltyp_and_shl2atm_from_fch(fchname, ncontr, shell_type, &
                                       shell2atom_map)
 call read_prim_per_shell(fchname, ncontr, prim_per_shell)

 call read_nprim_from_fch(fchname, nprim)
 allocate(prim_exp(nprim), contr_coeff(nprim))
 call read_prim_and_contr_coeff_from_fch(fchname, nprim, prim_exp, contr_coeff)
 if(has_pople_sp(fchname)) then
  allocate(contr_coeff_sp(nprim))
  call read_contr_coeff_sp_from_fch(fchname, nprim, contr_coeff_sp)
 end if

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(eigen_e_a(nif), alpha_coeff(nbf,nif), tot_dm(nbf,nbf))
 call read_eigenvalues_from_fch(fchname, nif, 'a', eigen_e_a)
 call read_mo_from_fch(fchname, nbf, nif, 'a', alpha_coeff)
 call read_dm_from_fch(fchname, 1, nbf, tot_dm)

 call check_uhf_in_fch(fchname, is_uhf)
 if(is_uhf) then
  allocate(eigen_e_b(nif), beta_coeff(nbf,nif), spin_dm(nbf,nbf))
  call read_eigenvalues_from_fch(fchname, nif, 'b', eigen_e_b)
  call read_mo_from_fch(fchname, nbf, nif, 'b', beta_coeff)
  call read_dm_from_fch(fchname, 2, nbf, spin_dm)
 end if

 call write_fch(new_fch)

 ! check whether ECP/PP is used
 call calc_ne_from_fch(fchname, ne1)
 call read_ne_from_fch(fchname, ne2)
 if(ne1 > ne2) then
  write(6,'(A)') REPEAT('-',79)
  write(6,'(A)') 'Warning from subroutine standardize_fch: ECP/PP is probably u&
                 &sed during calculation.'
  write(6,'(A)') 'The .fch(k) file generated by Q-Chem does not include any ECP&
                 &/PP data. So you'
  write(6,'(A)') 'need to write/add the ECP/PP information in the generated inp&
                 &ut file. For example:'
  write(6,'(/,A)') '#p RHF genecp nosymm int=nobasistransform guess=read geom=a&
                   &llcheck'
  write(6,'(/,A,/,A,/,A)') 'H C 0','cc-pVDZ','****'
  write(6,'(A,/,A,/,A)') 'I 0','SDD','****'
  write(6,'(/,A,/,A,//)') 'I 0', 'SDD'
  write(6,'(A)') REPEAT('-',79)
 end if
end subroutine standardize_fch

subroutine read_prim_per_shell(fchname, ncontr, prim_per_shell)
 implicit none
 integer :: i, fid
 integer, intent(in) :: ncontr
 integer, intent(out) :: prim_per_shell(ncontr)
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname

 prim_per_shell = 0
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:24) == 'Number of primitives per') exit
 end do ! for while

 read(fid,'(6(1X,I11))') (prim_per_shell(i), i=1,ncontr)
 close(fid)
end subroutine read_prim_per_shell

! Check whether there exists 'P(S=P) Contraction coefficients' section in a
! specified .fch file
! If this section exists but contains all zero elements, has_sp is still .T.
function has_pople_sp(fchname) result(has_sp)
 implicit none
 integer :: fid
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
 logical :: has_sp

 has_sp = .false.
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:6) == 'P(S=P)') then
   has_sp = .true.
   exit
  end if
  if(buf(1:19) == 'Coordinates of each') exit
 end do ! for while

 close(fid)
end function has_pople_sp

subroutine read_contr_coeff_sp_from_fch(fchname, nprim, contr_coeff_sp)
 implicit none
 integer :: i, fid
 integer, intent(in) :: nprim
 real(kind=8), intent(out) :: contr_coeff_sp(nprim)
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:6) == 'P(S=P)') exit
 end do ! for while

 if(i /= 0) then
  close(fid)
  write(6,'(/,A)') "ERROR in subroutine read_contr_coeff_sp_from_fch: no 'P(S=P&
                   &)' found in file"
  write(6,'(A)') TRIM(fchname)
  stop
 end if

 read(fid,'(5(1X,ES15.8))') (contr_coeff_sp(i),i=1,nprim)
 close(fid)
end subroutine read_contr_coeff_sp_from_fch

subroutine read_nprim_from_fch(fchname, nprim)
 implicit none
 integer :: fid
 integer, intent(out) :: nprim
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname

 nprim = 0
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:14) == 'Number of prim') exit
 end do ! for while

 BACKSPACE(fid)
 read(fid,'(A49,2X,I10)') buf, nprim
 close(fid)
end subroutine read_nprim_from_fch

subroutine read_prim_and_contr_coeff_from_fch(fchname, nprim, prim_exp, contr_coeff)
 implicit none
 integer :: i, fid
 integer, intent(in) :: nprim
 real(kind=8), intent(out) :: prim_exp(nprim), contr_coeff(nprim)
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:13) == 'Primitive exp') exit
 end do ! for while

 read(fid,'(5(1X,ES15.8))') (prim_exp(i),i=1,nprim)
 read(fid,'(A)') buf
 read(fid,'(5(1X,ES15.8))') (contr_coeff(i),i=1,nprim)
 close(fid)
end subroutine read_prim_and_contr_coeff_from_fch

! Calculate the total number of electrons using elements in a .fch(k) file.
! This value is not affected by ECP/PP.
subroutine calc_ne_from_fch(fchname, ne)
 implicit none
 integer :: charge, mult, natom
 integer, intent(out) :: ne
 integer, allocatable :: nuc(:)
 character(len=240), intent(in) :: fchname

 ne = 0
 call read_charge_and_mult_from_fch(fchname, charge, mult)
 call read_natom_from_fch(fchname, natom)
 allocate(nuc(natom))
 call read_nuc_from_fch(fchname, natom, nuc)
 ne = SUM(nuc) - charge
end subroutine calc_ne_from_fch

! read the type of HF from a .fch file
subroutine read_hf_type_from_fch(fchname, hf_type)
 use fch_content, only: check_ghf_in_fch, check_uhf_in_fch
 implicit none
 integer :: fid, na, nb
 integer, intent(out) :: hf_type
!f2py intent(out) :: hf_type
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname
 logical :: ghf, uhf

 hf_type = 1 ! default: real RHF, IOP(3/116=1)

 call check_ghf_in_fch(fchname, ghf)
 if(ghf) then
  hf_type = 7 ! complex GHF, IOP(3/116=7)
  return
 end if

 call check_uhf_in_fch(fchname, uhf)
 if(uhf) then
  hf_type = 2 ! real UHF, IOP(3/116=2)
  return
 end if

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:15) == 'Number of alpha') exit
 end do ! for while

 BACKSPACE(fid)
 read(fid,'(A49,2X,I10)') buf, na
 read(fid,'(A49,2X,I10)') buf, nb
 close(fid)

 if(na /= nb) hf_type = 101 ! real ROHF, IOP(3/116=101)
end subroutine read_hf_type_from_fch

! transfer Q-Chem AO-based density matrix from file 54.0 to Gaussian .fch(k)
! file
subroutine transfer_dm_qchem2gau(fchname, deleted)
 implicit none
 integer :: nbf, nif
 real(kind=8), allocatable :: dm(:,:)
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname
 logical, intent(in) :: deleted
!f2py intent(in) :: deleted

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(dm(nbf,nbf))
 call read_dm_from_qchem54(nbf, dm, deleted)
 call write_qchem_dm_into_fch(fchname, nbf, dm)
 deallocate(dm)
end subroutine transfer_dm_qchem2gau

! read Q-Chem AO-based density matrix from file 54.0 (assuming it is in the
!  current directory)
subroutine read_dm_from_qchem54(nbf, dm, deleted)
 implicit none
 integer :: i, fid
 integer, intent(in) :: nbf
 real(kind=8), intent(out) :: dm(nbf,nbf)
 logical, intent(in) :: deleted

 open(newunit=fid,file='54.0',status='old',access='stream')
 read(unit=fid,iostat=i) dm

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_dm_from_qchem54: failed to read fi&
                   &le 54.0'
  close(fid)
  stop
 else
  if(deleted) then
   close(fid,status='delete')
  else
   close(fid)
  end if
 end if

 dm = dm*2d0
end subroutine read_dm_from_qchem54

! write Q-Chem AO-based density matrix to a specified Gaussian .fch(k) file
subroutine write_qchem_dm_into_fch(fchname, nbf, dm)
 implicit none
 integer :: i, j, ncontr
 integer, intent(in) :: nbf
!f2py intent(in) :: nbf
 integer, allocatable :: shltyp(:), shl2atm(:)
 integer, allocatable :: idx(:)
 real(kind=8), intent(in) :: dm(nbf,nbf) ! in Q-Chem basis convention
!f2py intent(in) :: dm
!f2py depend(nbf) :: dm
 real(kind=8), allocatable :: dm0(:,:)
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname
 logical :: sph ! spherical harmonic or Cartesian-type basis functions

 call check_sph_in_fch(fchname, sph)
 call read_ncontr_from_fch(fchname, ncontr)
 allocate(shltyp(ncontr), shl2atm(ncontr))
 call read_shltyp_and_shl2atm_from_fch(fchname, ncontr, shltyp, shl2atm)
 deallocate(shl2atm)

 allocate(idx(nbf))
 call get_fch2qchem_permute_idx(sph, ncontr, shltyp, nbf, idx)
 deallocate(shltyp)

 ! the array idx contains Gaussian -> Qchem conversion relationship, here we
 ! need the inverse relationship
 allocate(dm0(nbf,nbf))
 forall(i=1:nbf, j=1:nbf) dm0(idx(j),idx(i)) = dm(j,i)
 deallocate(idx)

 call write_dm_into_fch(fchname, .true., nbf, dm0)
 deallocate(dm0)
end subroutine write_qchem_dm_into_fch

!subroutine freq_dep_diag(npyname, nacto, nacte, init_e)
! implicit none
! integer :: i, j, nocc_a, nocc_b, nvir_a, nvir_b, n_cisd, nlarge, nconf
! integer, parameter :: max_it = 99
! integer, intent(in) :: nacto, nacte
!!f2py intent(in) :: nacto, nacte
! integer(kind=4) :: t0, t1, TIME
! real(kind=8) :: old_e
! real(kind=8), parameter :: thres = 1d-7
! real(kind=8), allocatable :: ham(:,:), R_w(:,:), x(:,:), S0(:,:), S(:,:), &
!  C(:,:), CCT(:,:), CR_w(:,:), w(:)
! real(kind=8), intent(in) :: init_e
!!f2py intent(in) :: init_e
! character(len=240), intent(in) :: npyname
!!f2py intent(in) :: npyname
!
! t0 = TIME()
! nocc_b = nacte/2        ! number of beta occupied spin orbitals
! nocc_a = nacte - nocc_b ! number of alpha occupied spin orbitals
! nvir_a = nacto - nocc_a ! number of alpha unoccupied spin orbitals
! nvir_b = nacto - nocc_b ! number of beta unoccupied spin orbitals
! i = nocc_a*nvir_a; j = nocc_b*nvir_b
! n_cisd = 1 + i*j + (i*(i-nacto+5) + j*(j-nacto+5))/4
!
! call read_sym_mat_size_from_npy(npyname, nconf)
! write(6,'(A,2I7)') 'n_cisd, nconf=', n_cisd, nconf
! if(n_cisd >= nconf) then
!  write(6,'(/,A)') 'ERROR in subroutine freq_dep_diag: n_cisd >= nconf.'
!  write(6,'(A)') 'The excitation level is <= CISD level.'
!  write(6,'(A)') 'npyname='//TRIM(npyname)
!  stop
! end if
!
! allocate(ham(nconf,nconf), w(n_cisd))
! call read_sym_mat_from_npy(npyname, nconf, ham)
! i = n_cisd + 1
! nlarge = nconf - n_cisd
! allocate(S0(n_cisd,n_cisd), source=ham(1:n_cisd,1:n_cisd))
! allocate(C(n_cisd,nlarge), source=ham(1:n_cisd,i:nconf))
! allocate(R_w(nlarge,nlarge), source=ham(i:nconf,i:nconf))
! deallocate(ham)
! allocate(CCT(n_cisd,n_cisd), source=0d0)
! call dgemm('N', 'T', n_cisd, n_cisd, nlarge, 1d0, C, n_cisd, C, n_cisd, 0d0, &
!            CCT, n_cisd)
!
! old_e = init_e
! write(6,'(A,F18.12)') 'i=  0, ini_e = ', init_e
! t1 = TIME()
! write(6,'(I3)') t1-t0
! t0 = t1
!
! do i = 1, max_it, 1
!  ! R-w
!  forall(j = 1:nlarge) R_w(j,j) = R_w(j,j) - old_e
!  ! C(R-w)
!  allocate(CR_w(n_cisd,nlarge), source=0d0)
!  call dgemm('N', 'N', n_cisd, nlarge, nlarge, 1d0, C, n_cisd, R_w, nlarge, &
!             0d0, CR_w, n_cisd)
!  ! (R-w)^(-1) (C^T)
!  allocate(x(nlarge,n_cisd))
!  call solve_multi_lin_eqs(n_cisd, nlarge, CR_w, n_cisd, CCT, x)
!  deallocate(CR_w)
!  t1 = TIME()
!  write(6,'(/,I3)') t1-t0
!  t0 = t1
!  ! S - C (R-w)^(-1) (C^T)
!  allocate(S(n_cisd,n_cisd), source=S0)
!  call dgemm('N', 'N', n_cisd, n_cisd, nlarge, 1d0, C, n_cisd, -x, nlarge, 1d0,&
!             S, n_cisd)
!  deallocate(x)
!  ! diagonalize (S - C1 (R-w)^(-1) C)
!  call diag_get_e_and_vec(n_cisd, S, w)
!  deallocate(S)
!  write(6,'(A,I3,A,F18.12)') 'i=', i, ', new_e = ', w(1)
!  if(DABS(w(1)-old_e) < thres) exit
!  forall(j = 1:nlarge) R_w(j,j) = R_w(j,j) + old_e ! recover R
!  old_e = w(1) ! update old_e
!  t1 = TIME()
!  write(6,'(I3)') t1-t0
!  t0 = t1
! end do ! for i
!
! deallocate(R_w, S0, C, CCT, w)
!end subroutine freq_dep_diag

!subroutine freq_dep_diag2(npyname, nacto, nacte, init_e)
! implicit none
! integer :: i, j, nocc_a, nocc_b, nvir_a, nvir_b, n_cisd, nlarge, nconf
! integer, parameter :: max_it = 99
! integer, intent(in) :: nacto, nacte
!!f2py intent(in) :: nacto, nacte
! integer(kind=4) :: t0, t1, TIME
! real(kind=8) :: old_e
! real(kind=8), parameter :: thres = 1d-7
! real(kind=8), allocatable :: ham(:,:), R_w(:,:), x(:,:), S0(:,:), S(:,:), &
!  C(:,:), w(:)
! real(kind=8), intent(in) :: init_e
!!f2py intent(in) :: init_e
! character(len=240), intent(in) :: npyname
!!f2py intent(in) :: npyname
!
! t0 = TIME()
! nocc_b = nacte/2        ! number of beta occupied spin orbitals
! nocc_a = nacte - nocc_b ! number of alpha occupied spin orbitals
! nvir_a = nacto - nocc_a ! number of alpha unoccupied spin orbitals
! nvir_b = nacto - nocc_b ! number of beta unoccupied spin orbitals
! i = nocc_a*nvir_a; j = nocc_b*nvir_b
! n_cisd = 1 + i*j + (i*(i-nacto+5) + j*(j-nacto+5))/4
!
! call read_sym_mat_size_from_npy(npyname, nconf)
! write(6,'(A,2I7)') 'n_cisd, nconf=', n_cisd, nconf
! if(n_cisd >= nconf) then
!  write(6,'(/,A)') 'ERROR in subroutine freq_dep_diag: n_cisd >= nconf.'
!  write(6,'(A)') 'The excitation level is <= CISD level.'
!  write(6,'(A)') 'npyname='//TRIM(npyname)
!  stop
! end if
!
! allocate(ham(nconf,nconf), w(n_cisd))
! call read_sym_mat_from_npy(npyname, nconf, ham)
! i = n_cisd + 1
! nlarge = nconf - n_cisd
! allocate(S0(n_cisd,n_cisd), source=ham(1:n_cisd,1:n_cisd))
! allocate(C(n_cisd,nlarge), source=ham(1:n_cisd,i:nconf))
! allocate(R_w(nlarge,nlarge), source=ham(i:nconf,i:nconf))
! deallocate(ham)
!
! old_e = init_e
! write(6,'(A,F18.12)') 'i=  0, ini_e = ', init_e
! t1 = TIME()
! write(6,'(I3)') t1-t0
! t0 = t1
!
! do i = 1, max_it, 1
!  ! R-w
!  forall(j = 1:nlarge) R_w(j,j) = R_w(j,j) - old_e
!  ! (R-w)^(-1) (C^T)
!  allocate(x(nlarge,n_cisd))
!  call solve_multi_lin_eqs(nlarge, nlarge, R_w, n_cisd, TRANSPOSE(C), x)
!  t1 = TIME()
!  write(6,'(/,I3)') t1-t0
!  t0 = t1
!  ! S - C (R-w)^(-1) (C^T)
!  allocate(S(n_cisd,n_cisd), source=S0)
!  call dgemm('N', 'N', n_cisd, n_cisd, nlarge, 1d0, C, n_cisd, -x, nlarge, 1d0,&
!             S, n_cisd)
!  deallocate(x)
!  ! diagonalize (S - C1 (R-w)^(-1) C)
!  call diag_get_e_and_vec(n_cisd, S, w)
!  deallocate(S)
!  write(6,'(A,I3,A,F18.12)') 'i=', i, ', new_e = ', w(1)
!  if(DABS(w(1)-old_e) < thres) exit
!  forall(j = 1:nlarge) R_w(j,j) = R_w(j,j) + old_e ! recover R
!  old_e = w(1) ! update old_e
!  t1 = TIME()
!  write(6,'(I3)') t1-t0
!  t0 = t1
! end do ! for i
!
! deallocate(R_w, S0, C, w)
!end subroutine freq_dep_diag2

!subroutine freq_dep_diag3(npyname, nacto, nacte, init_e)
! implicit none
! integer :: i, j, nocc_a, nocc_b, nvir_a, nvir_b, n_cisd, nlarge, nconf
! integer, parameter :: max_it = 99
! integer, intent(in) :: nacto, nacte
!!f2py intent(in) :: nacto, nacte
! integer(kind=4) :: t0, t1, TIME
! real(kind=8) :: old_e
! real(kind=8), parameter :: thres = 1d-7
! real(kind=8), allocatable :: ham(:,:), w(:), R_w(:,:), inv_R_w(:,:), S0(:,:), &
!  S(:,:), C(:,:), inv_R_wCT(:,:)
! real(kind=8), intent(in) :: init_e
!!f2py intent(in) :: init_e
! character(len=240), intent(in) :: npyname
!!f2py intent(in) :: npyname
!
! t0 = TIME()
! nocc_b = nacte/2        ! number of beta occupied spin orbitals
! nocc_a = nacte - nocc_b ! number of alpha occupied spin orbitals
! nvir_a = nacto - nocc_a ! number of alpha unoccupied spin orbitals
! nvir_b = nacto - nocc_b ! number of beta unoccupied spin orbitals
! i = nocc_a*nvir_a; j = nocc_b*nvir_b
! n_cisd = 1 + i*j + (i*(i-nacto+5) + j*(j-nacto+5))/4
!
! call read_sym_mat_size_from_npy(npyname, nconf)
! write(6,'(A,2I7)') 'n_cisd, nconf=', n_cisd, nconf
! if(n_cisd >= nconf) then
!  write(6,'(/,A)') 'ERROR in subroutine freq_dep_diag: n_cisd >= nconf.'
!  write(6,'(A)') 'The excitation level is <= CISD level.'
!  write(6,'(A)') 'npyname='//TRIM(npyname)
!  stop
! end if
!
! allocate(ham(nconf,nconf), w(n_cisd))
! call read_sym_mat_from_npy(npyname, nconf, ham)
! i = n_cisd + 1
! nlarge = nconf - n_cisd
! allocate(S0(n_cisd,n_cisd), source=ham(1:n_cisd,1:n_cisd))
! allocate(C(n_cisd,nlarge), source=ham(1:n_cisd,i:nconf))
! allocate(R_w(nlarge,nlarge), source=ham(i:nconf,i:nconf))
! deallocate(ham)
!
! old_e = init_e
! write(6,'(A,F18.12)') 'i=  0, ini_e = ', init_e
! t1 = TIME()
! write(6,'(I3)') t1-t0
! t0 = t1
!
! do i = 1, max_it, 1
!  ! (R-w)^(-1)
!  forall(j = 1:nlarge) R_w(j,j) = R_w(j,j) - old_e
!  allocate(inv_R_w(nlarge,nlarge))
!  !call inverse(nlarge, R_w, inv_R_w)
!  call newton_inv(nlarge, R_w, inv_R_w)
!  t1 = TIME()
!  write(6,'(/,I3)') t1-t0
!  t0 = t1
!  ! (R-w)^(-1) (C^T)
!  allocate(inv_R_wCT(nlarge,n_cisd), source=0d0)
!  call dsymm('L', 'L', nlarge, n_cisd, 1d0, inv_R_w, nlarge, TRANSPOSE(C), &
!             nlarge, 0d0, inv_R_wCT, nlarge)
!  deallocate(inv_R_w)
!  ! S - C (R-w)^(-1) (C^T)
!  allocate(S(n_cisd,n_cisd), source=S0)
!  call dgemm('N', 'N', n_cisd, n_cisd, nlarge, 1d0, C, n_cisd, -inv_R_wCT, &
!             nlarge, 1d0, S, n_cisd)
!  deallocate(inv_R_wCT)
!  ! diagonalize (S - C1 (R-w)^(-1) C2)
!  call diag_get_e_and_vec(n_cisd, S, w)
!  deallocate(S)
!  write(6,'(A,I3,A,F18.12)') 'i=', i, ', new_e = ', w(1)
!  if(DABS(w(1)-old_e) < thres) exit
!  forall(j = 1:nlarge) R_w(j,j) = R_w(j,j) + old_e
!  old_e = w(1)
!  t1 = TIME()
!  write(6,'(I3)') t1-t0
!  t0 = t1
! end do ! for i
!
! deallocate(S0, C, R_w, w)
!end subroutine freq_dep_diag3

