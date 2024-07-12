! written by jxzou at 20240520
! Note: Amesp currently requires that the atoms which share the same elements
!  use the same basis set.

program main
 use fch_content, only: check_uhf_in_fch
 implicit none
 integer :: i, j, itype ! itype=0/1/2 for R(O)HF/UHF/Some kind of NOs
 character(len=3) :: str3
 character(len=240) :: amoname, fchname
 logical :: alive, uhf

 i = iargc()
 if(i<1 .or. i>3) then
  write(6,'(/,A)') ' ERROR in subroutine amo2fch: wrong command line arguments!'
  write(6,'(A)')   ' Example 1 (for R(O)HF/UHF/CAS): amo2fch a.amo'
  write(6,'(A)')   ' Example 2 (for R(O)HF/UHF/CAS): amo2fch a.amo b.fch'
  write(6,'(A,/)') ' Example 3 (for CAS NOs)       : amo2fch a.amo a.fch -no'
  stop
 end if

 itype = 0; alive = .false.; str3 = ' '; amoname = ' '; fchname = ' '
 call getarg(1, amoname)
 call require_file_exist(amoname)

 if(i == 1) then
  call find_specified_suffix(amoname, '.amo', j)
  fchname = amoname(1:j-1)//'.fch'
 else
  call getarg(2, fchname)
 end if

 inquire(file=TRIM(fchname),exist=alive)
 if(alive) then
  write(6,'(/,A)') 'Remark from program amo2fch: file '//TRIM(fchname)//' exists.'
  write(6,'(A)') 'It will be used directly.'
 else
  write(6,'(/,A)') 'Warning from program amo2fch: file '//TRIM(fchname)//' does&
                   & not exist,'
  write(6,'(A)') 'the program is trying to generate one from scratch...'
  call gen_fch_from_amo(amoname, fchname)
 end if

 call check_uhf_in_fch(fchname, uhf)
 if(uhf) itype = 1

 if(i == 3) then
  call getarg(3, str3)
  if(str3 /= '-no') then
   write(6,'(/,A)') "ERROR in subroutine amo2fch: the 3rd argument can only be &
                    &'-no'."
   stop
  end if
  itype = 2
 end if

 if(uhf .and. itype>1) then
  write(6,'(/,A)') "ERROR in subroutine amo2fch: '-no' argument is imcompatible&
                   & with a UHF-type"
  write(6,'(A)') '.fch file. An R(O)HF-type .fch file is required in such case.'
  stop
 end if

 call amo2fch(amoname, fchname, itype)
 if(.not. alive) write(6,'(/,A,/)') 'Done generation.'
end program main

! write/export/transfer MOs from Amesp .amo to Gaussian .fch
subroutine amo2fch(amoname, fchname, itype)
 implicit none
 integer :: na, nb, nbf, nif, nbf1, nif1, charge, mult
 integer, intent(in) :: itype ! 0/1/2 for R(O)HF/UHF/Some kind of NOs
 real(kind=8), allocatable :: e_a(:), e_b(:), mo_a(:,:), mo_b(:,:), &
  dm_a(:,:), dm_b(:,:), noon(:)
 character(len=240), intent(in) :: amoname, fchname
 logical :: uhf

 call read_charge_and_mult_from_amo(amoname, charge, mult)
 call read_nbf_and_nif_from_amo(amoname, nbf1, nif1)
 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 if(nbf/=nbf1 .or. nif/=nif1) then
  write(6,'(/,A)') 'ERROR in subroutine amo2fch: the array size of MOs is incon&
                   &sistent.'
  write(6,'(A,4I6)') 'nbf1, nif1, nbf, nif=', nbf1, nif1, nbf, nif
  write(6,'(A)') 'amoname='//TRIM(amoname)
  write(6,'(A)') 'fchname='//TRIM(fchname)
  stop
 end if

 call check_uhf_in_amo(amoname, uhf)
 if(uhf .and. itype/=1) then
  write(6,'(/,A)') 'ERROR in subroutine amo2fch: the wave function is UHF-type i&
                   &n .amo, but R(O)HF-'
  write(6,'(A)') 'type in .fch file. Inconsistent.'
  write(6,'(A)') 'amoname='//TRIM(amoname)
  write(6,'(A)') 'fchname='//TRIM(fchname)
  stop
 else if((.not. uhf) .and. itype==1) then
  write(6,'(/,A)') 'ERROR in subroutine amo2fch: the wave function is UHF-type i&
                   &n .fch, but R(O)HF-'
  write(6,'(A)') 'type in .amo file. Inconsistent.'
  write(6,'(A)') 'amoname='//TRIM(amoname)
  write(6,'(A)') 'fchname='//TRIM(fchname)
  stop
 end if

 allocate(e_a(nif), mo_a(nbf,nif))
 call read_ev_from_amo(amoname, nif, 'a', e_a)
 call read_mo_from_amo(amoname, nbf, nif, 'a', mo_a)
 if(itype == 1) then
  allocate(e_b(nif), mo_b(nbf,nif))
  call read_ev_from_amo(amoname, nif, 'b', e_b)
  call read_mo_from_amo(amoname, nbf, nif, 'b', mo_b)
 end if

 call read_na_and_nb_from_amo(amoname, na, nb)
 allocate(noon(nif), source=0d0)
 noon(1:na) = 1d0
 allocate(dm_a(nbf,nbf))
 call calc_dm_using_mo_and_on(nbf, nif, mo_a, noon, dm_a)

 if(itype == 1) then ! UHF
  allocate(dm_b(nbf,nbf))
  if(na > nb) noon(nb+1:na) = 0d0
  call calc_dm_using_mo_and_on(nbf, nif, mo_b, noon, dm_b)
  call write_dm_into_fch(fchname, .true., nbf, dm_a+dm_b)
  call write_dm_into_fch(fchname, .false., nbf, dm_a-dm_b)
  deallocate(dm_b)
 else                ! R(O)HF
  if(mult == 1) then ! RHF
   call write_dm_into_fch(fchname, .true., nbf, 2d0*dm_a)
  else               ! ROHF
   allocate(dm_b(nbf,nbf))
   noon(nb+1:na) = 0d0
   call calc_dm_using_mo_and_on(nbf, nif, mo_a, noon, dm_b)
   call write_dm_into_fch(fchname, .true., nbf, dm_a+dm_b)
   call write_dm_into_fch(fchname, .false., nbf, dm_a-dm_b)
   deallocate(dm_b)
  end if
 end if
 deallocate(noon, dm_a)

 call write_eigenvalues_to_fch(fchname, nif, 'a', e_a, .true.)
 call write_mo_into_fch(fchname, nbf, nif, 'a', mo_a)
 deallocate(e_a, mo_a)
 if(itype == 1) then
  call write_eigenvalues_to_fch(fchname, nif, 'b', e_b, .true.)
  call write_mo_into_fch(fchname, nbf, nif, 'b', mo_b)
  deallocate(e_b, mo_b)
 end if
end subroutine amo2fch

! generate a .fch file from an Amesp .amo file
subroutine gen_fch_from_amo(amoname, fchname)
 use fch_content
 use mkl_content, only: natom1=>natom, ncontr1=>ncontr,shell_type1=>shell_type,&
  shl2atm1=>shl2atm, elem1=>elem, all_pg, read_all_pg_from_amo, &
  un_normalized_all_pg, merge_s_and_p_into_sp
 implicit none
 integer :: i, necpatm
 integer :: ntype ! the number of types of atoms, e.g. 2 for H2O
 integer, allocatable :: ielem0(:)
 character(len=240), intent(in) :: amoname, fchname
 logical :: sfx2c1e, has_sp, ecp

 call read_charge_and_mult_from_amo(amoname, charge, mult)
 call read_na_and_nb_from_amo(amoname, na, nb)
 nopen = na - nb
 call read_nbf_and_nif_from_amo(amoname, nbf, nif)
 call check_sfx2c1e_in_amo(amoname, sfx2c1e)
 if(sfx2c1e) irel = -3
 call check_uhf_in_amo(amoname, is_uhf)
 call check_ecp_in_amo(amoname, ecp)
 call read_natom_from_amo(amoname, natom)
 allocate(elem(natom), ielem(natom), coor(3,natom))
 call read_elem_nuc_coor_from_amo(amoname, natom, elem, ielem, coor)

 call read_ncontr_from_amo(amoname, natom, elem, ncontr, ntype)
 allocate(shell_type(ncontr), prim_per_shell(ncontr), shell2atom_map(ncontr))
 call read_shltyp_prpshl_shl2atm_from_amo(amoname, natom, ntype, ncontr, elem,&
                                    shell_type, prim_per_shell, shell2atom_map)

 natom1 = natom; ncontr1 = ncontr
 allocate(shell_type1(ncontr), source=shell_type)
 allocate(shl2atm1(ncontr), source=shell2atom_map)
 allocate(elem1(natom), source=elem)
 call read_all_pg_from_amo(amoname, ntype, prim_per_shell)
 call un_normalized_all_pg()
 call merge_s_and_p_into_sp()
 if(ncontr1 < ncontr) then
  deallocate(shell_type, shell2atom_map)
  allocate(shell_type(ncontr1), source=shell_type1)
  allocate(shell2atom_map(ncontr1), source=shl2atm1)
  ncontr = ncontr1
 end if

 deallocate(prim_per_shell)
 allocate(prim_per_shell(ncontr))
 call find_nprim_from_all_pg(ncontr, prim_per_shell, nprim, has_sp)
 call all_pg2prim_exp_and_contr_coeff(has_sp)
 deallocate(all_pg)

 if(ecp) then
  call find_LenNCZ_in_amo(amoname, LenNCZ, necpatm)
  allocate(KFirst(natom,10), KLast(natom,10), Lmax(natom), LPSkip(natom), &
   NLP(LenNCZ), RNFroz(natom), CLP(LenNCZ), ZLP(LenNCZ), ielem0(natom))
  forall(i = 1:natom) ielem0(i) = elem2nuc(elem(i))
  forall(i = 1:natom) RNFroz(i) = DBLE(ielem0(i) - ielem(i))
  ielem = ielem0
  deallocate(ielem0)
  call read_ecp_arrays_from_amo(amoname, natom, LenNCZ, necpatm, KFirst, KLast,&
                                Lmax, LPSkip, NLP, CLP, ZLP)
 end if

 allocate(eigen_e_a(nif), source=0d0)
 allocate(alpha_coeff(nbf,nif), source=0d0)
 allocate(tot_dm(nbf,nbf), source=0d0)
 if(is_uhf) then
  allocate(eigen_e_b(nif), source=0d0)
  allocate(beta_coeff(nbf,nif), source=0d0)
 end if

 if(is_uhf .or. ((.not.is_uhf) .and. mult>1)) then
  allocate(spin_dm(nbf,nbf), source=0d0)
 end if
 call write_fch(fchname)
 call free_arrays_in_fch_content()
end subroutine gen_fch_from_amo

! read charge and spin multiplicity from a specified .amo file
subroutine read_charge_and_mult_from_amo(amoname, charge, mult)
 implicit none
 integer :: i, fid
 integer, intent(out) :: charge, mult
 character(len=10) :: buf ! long buf is not needed here
 character(len=240), intent(in) :: amoname

 charge = 0; mult = 1
 open(newunit=fid,file=TRIM(amoname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:8) == '[Charge]') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_charge_and_mult_from_amo: charge n&
                   &ot found in'
  write(6,'(A)') 'file '//TRIM(amoname)
  write(6,'(A)') 'Please make sure that your Amesp version >= May 23, 2024.'
  close(fid)
  stop
 end if

 read(fid,*) charge
 read(fid,'(A)') buf
 if(buf(1:6) /= '[Mult]') then
  write(6,'(/,A)') 'ERROR in subroutine read_charge_and_mult_from_amo: [Mult] i&
                   &s not below [Charge].'
  write(6,'(A)') 'Problematic file: '//TRIM(amoname)
  close(fid)
  stop
 end if

 read(fid,*) mult
 close(fid)
end subroutine read_charge_and_mult_from_amo

! read the number of Alpha/Beta electrons from a specified .amo file
subroutine read_na_and_nb_from_amo(amoname, na, nb)
 implicit none
 integer :: i, fid
 integer, intent(out) :: na, nb
 character(len=10) :: buf
 character(len=240), intent(in) :: amoname

 na = 0; nb = 0
 open(newunit=fid,file=TRIM(amoname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:4) == 'Nocc') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_na_and_nb_from_amo: Nocc not found&
                   & in file '//TRIM(amoname)
  close(fid)
  stop
 end if

 i = INDEX(buf, '=')
 read(buf(i+1:),*) na

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:5) == 'NoccB') exit
 end do ! for while

 close(fid)
 if(i == 0) then
  i = INDEX(buf, '=')
  read(buf(i+1:),*) nb
 else
  nb = na
 end if
end subroutine read_na_and_nb_from_amo

! read nbf and nif from a specified .amo file
subroutine read_nbf_and_nif_from_amo(amoname, nbf, nif)
 implicit none
 integer :: i, nline, ncol, fid
 integer, intent(out) :: nbf, nif
 integer, external :: detect_ncol_in_buf
 character(len=240) :: buf0, buf
 character(len=240), intent(in) :: amoname

 nbf = 0; nif = 0
 open(newunit=fid,file=TRIM(amoname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:3)=='En:' .or. buf(1:4)=='EnA:') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_nbf_and_nif_from_amo: 'En:'/'EnA:'&
                   & not found in"
  write(6,'(A)') 'file '//TRIM(amoname)
  close(fid)
  stop
 end if

 nline = 0
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:4) == 'MoCu') exit
  nline = nline + 1
  buf0 = buf
 end do ! for while

 ncol = detect_ncol_in_buf(buf0)
 nif = (nline-1)*5 + ncol

 nline = 0
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:5) == 'NoccB') exit
  nline = nline + 1
  buf0 = buf
 end do ! for while

 ncol = detect_ncol_in_buf(buf0)
 nbf = ((nline-1)*5 + ncol)/nif
 close(fid)
end subroutine read_nbf_and_nif_from_amo

! check if sfX2C1e is True in a specified .amo file
subroutine check_sfx2c1e_in_amo(amoname, sfx2c1e)
 implicit none
 integer :: i, fid
 character(len=20) :: buf
 character(len=240), intent(in) :: amoname
 logical, intent(out) :: sfx2c1e

 sfx2c1e = .false.
 open(newunit=fid,file=TRIM(amoname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) then
   close(fid)
   return
  end if
  if(buf(1:7) == 'sfx2c1e') exit
 end do ! for while

 close(fid)
 i = INDEX(buf, ':')
 read(buf(i+1:),*) sfx2c1e
end subroutine check_sfx2c1e_in_amo

! check whether the wave function is UHF-type in a specified .amo file
subroutine check_uhf_in_amo(amoname, uhf)
 implicit none
 integer :: i, fid
 character(len=2) :: str2 = '  '
 character(len=3) :: str3 = '   '
 character(len=25) :: buf
 character(len=240), intent(in) :: amoname
 logical, intent(out) :: uhf

 uhf = .false.
 open(newunit=fid,file=TRIM(amoname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:4) == '[MO]') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine check_uhf_in_amo: [MO] not found in fil&
                   &e '//TRIM(amoname)
  stop
 end if

 read(buf(5:),*) str3, str2

 select case(TRIM(str2))
 case('R','RO') ! do nothing
 case('U')
  uhf = .true.
 case default
  write(6,'(/,A)') 'ERROR in subroutine check_uhf_in_amo: wave function type R/&
                   &RO/U cannot be'
  write(6,'(A)') 'recognized in '//TRIM(amoname)
  write(6,'(A)') TRIM(str2)
  stop
 end select
end subroutine check_uhf_in_amo

! check whether ECP/PP is used in a specified .amo file
subroutine check_ecp_in_amo(amoname, ecp)
 implicit none
 integer :: i, fid
 character(len=10) :: buf
 character(len=240), intent(in) :: amoname
 logical, intent(out) :: ecp

 ecp = .false.
 open(newunit=fid,file=TRIM(amoname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:3) == 'ECP') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine check_ecp_in_amo: ECP section is not fo&
                   &und in file '//TRIM(amoname)
  stop
 end if
 read(buf(5:),*) ecp
end subroutine check_ecp_in_amo

! read Alpha/Beta orbital energies from a specified .amo file
subroutine read_ev_from_amo(amoname, nif, ab, ev)
 implicit none
 integer :: i, fid
 integer, intent(in) :: nif
 real(kind=8), intent(out) :: ev(nif)
 character(len=1), intent(in) :: ab
 character(len=10) :: buf ! long buf is not needed here
 character(len=240), intent(in) :: amoname

 ev = 0d0
 open(newunit=fid,file=TRIM(amoname),status='old',position='rewind')

 select case(ab)
 case('a') ! Alpha MOs
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:3)=='En:' .or. buf(1:4)=='EnA:') exit
  end do ! for while
 case('b') ! Beta MOs
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:4) == 'EnB:') exit
  end do ! for while
 case default
  write(6,'(/,A)') 'ERROR in subroutine read_ev_from_amo: MO type cannot be rec&
                   &ognized. ab='//ab
  close(fid)
  stop
 end select

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_ev_from_amo: specified orbital eig&
                   &envalue section'
  write(6,'(A)') 'not found in file '//TRIM(amoname)
  close(fid)
  stop
 end if

 read(fid,*,iostat=i) ev
 close(fid)

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_ev_from_amo: failed to read orbita&
                   &l eigenvalues.'
  write(6,'(A)') 'This file seems problematic: '//TRIM(amoname)
  stop
 end if
end subroutine read_ev_from_amo

! read Alpha/Beta MOs from a specified .amo file
subroutine read_mo_from_amo(amoname, nbf, nif, ab, mo)
 implicit none
 integer :: i, fid
 integer, intent(in) :: nbf, nif
 real(kind=8), intent(out) :: mo(nbf,nif)
 character(len=1), intent(in) :: ab
 character(len=10) :: buf ! long buf is not needed here
 character(len=240), intent(in) :: amoname

 mo = 0d0
 open(newunit=fid,file=TRIM(amoname),status='old',position='rewind')

 select case(ab)
 case('a') ! Alpha MOs
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:5)=='MoCu:' .or. buf(1:6)=='MoCuA:') exit
  end do ! for while
 case('b') ! Beta MOs
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:6) == 'MoCuB:') exit
  end do ! for while
 case default
  write(6,'(/,A)') 'ERROR in subroutine read_mo_from_amo: MO type cannot be rec&
                   &ognized. ab='//ab
  close(fid)
  stop
 end select

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_mo_from_amo: specified MO section &
                   &not found in'
  write(6,'(A)') 'file '//TRIM(amoname)
  close(fid)
  stop
 end if

 read(fid,*,iostat=i) mo
 close(fid)

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_mo_from_amo: failed to read MOs.'
  write(6,'(A)') 'This file seems problematic: '//TRIM(amoname)
  stop
 end if
end subroutine read_mo_from_amo

! read elements, atomic order and Cartesian coordinates from an Amesp .amo file
subroutine read_elem_nuc_coor_from_amo(amoname, natom, elem, nuc, coor)
 use phys_cons, only: Bohr_const
 implicit none
 integer :: i, fid
 integer, intent(in) :: natom
 integer, intent(out) :: nuc(natom)
 real(kind=8), intent(out) :: coor(3,natom)
 character(len=2), intent(out) :: elem(natom)
 character(len=240) :: buf
 character(len=240), intent(in) :: amoname

 elem = ' '; nuc = 0; coor = 0d0
 open(newunit=fid,file=TRIM(amoname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:7) == '[Atoms]') exit
 end do ! for while

 read(fid,'(A)') buf
 do i = 1, natom, 1
  read(fid,*) elem(i), nuc(i), coor(:,i)
 end do ! for i

 close(fid)
 coor = coor*Bohr_const
end subroutine read_elem_nuc_coor_from_amo

! Read the Number of contracted shells from a specified .amo file. Here we
! assume that the elements are known and provided.
subroutine read_ncontr_from_amo(amoname, natom, elem, ncontr, ntype)
 implicit none
 integer :: i, j, k, nline, fid
 integer, intent(in) :: natom
 integer, intent(out) :: ncontr, ntype
 character(len=240) :: buf
 character(len=2) :: str2
 character(len=2), intent(in) :: elem(natom)
 character(len=240), intent(in) :: amoname

 ncontr = 0; ntype = 0
 open(newunit=fid,file=TRIM(amoname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:10) == 'NAtom_Type') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_ncontr_from_amo: no NAtom_Type fou&
                   &nd in file '//TRIM(amoname)
  close(fid)
  stop
 end if

 i = INDEX(buf, ':')
 read(buf(i+1:),*) ntype

 do i = 1, 3
  read(fid,'(A)') buf
 end do

 do i = 1, ntype, 1
  read(fid,*) str2 ! e.g. 'O     8'
  str2 = ADJUSTL(str2)
  read(fid,'(A)') buf
  j = INDEX(buf, ':')
  read(buf(j+1:),*) k
  ncontr = ncontr + k*COUNT(elem==str2)

  do while(.true.)
   read(fid,'(A)') buf
   if(buf(2:6) == 'Bmnl:') exit
  end do ! for while
  read(buf(7:),*) k
  nline = k/5
  if(k - 5*nline > 0) nline = nline + 1
  do j = 1, nline, 1
   read(fid,'(A)') buf
  end do ! for j
 end do ! for i

 close(fid)
end subroutine read_ncontr_from_amo

! read 3 arrays shltyp, prpshl and shl2atm from a specified .amo file
subroutine read_shltyp_prpshl_shl2atm_from_amo(amoname, natom, ntype, ncontr, &
                                               elem, shltyp, prpshl, shl2atm)
 implicit none
 integer :: i, j, k, m, maxnshl, nline, fid
 integer, intent(in) :: natom, ntype, ncontr
 integer, intent(out) :: shltyp(ncontr), prpshl(ncontr), shl2atm(ncontr)
 integer, allocatable :: atmnshl(:), angshl(:,:), bas_num(:,:)
 character(len=2) :: str2
 character(len=3) :: str3
 character(len=2), intent(in) :: elem(natom)
 character(len=2), allocatable :: elem1(:)
 character(len=240) :: buf
 character(len=240), intent(in) :: amoname
 logical :: sph ! spherical harmonic or Cartesian

 str2 = ' '; str3 = ' '; shltyp = 0; prpshl = 0; shl2atm = 0
 allocate(atmnshl(ntype), source=0)
 allocate(elem1(ntype))
 elem1 = ' '

 open(newunit=fid,file=TRIM(amoname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:7) == 'maxNShl') exit
 end do ! for while

 i = INDEX(buf, ':')
 read(buf(i+1:),*) maxnshl
 allocate(angshl(maxnshl,ntype), source=0)
 allocate(bas_num(maxnshl,ntype), source=0)
 read(fid,'(A)') buf
 read(fid,'(A)') buf

 do i = 1, ntype, 1
  read(fid,*) str2 ! e.g. 'O     8'
  str2 = ADJUSTL(str2)
  elem1(i) = str2
  read(fid,'(A)') buf
  j = INDEX(buf, ':')
  read(buf(j+1:),*) k
  atmnshl(i) = k
  read(fid,'(A)') buf
  read(fid,*) angshl(1:k,i)
  read(fid,'(A)') buf
  read(fid,*) bas_num(1:k,i)

  do while(.true.)
   read(fid,'(A)') buf
   if(buf(2:6) == 'Bmnl:') exit
  end do ! for while
  read(buf(7:),*) k
  nline = k/5
  if(k - 5*nline > 0) nline = nline + 1
  do j = 1, nline, 1
   read(fid,'(A)') buf
  end do ! for j
 end do ! for i

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:4) == '[MO]') exit
 end do ! for while
 close(fid)
 read(buf(5:),*) str3

 select case(str3)
 case('sph') ! spherical harmonic 5D 7F
  sph = .true.
 case('car') ! Cartesian 6D 10F
  sph = .false.
 case default
  write(6,'(/,A)') 'ERROR in subroutine read_shltyp_prpshl_shl2atm_from_amo: ba&
                   &sis function type'
  write(6,'(A)') "cannot be recognized: '"//str3//"'"
  stop
 end select
 m = 0

 do i = 1, natom, 1
  str2 = elem(i)
  do j = 1, ntype, 1
   if(elem1(j) == str2) exit
  end do ! for j
  k = atmnshl(j)
  shltyp(m+1:m+k) = angshl(1:k,j)
  prpshl(m+1:m+k) = bas_num(1:k,j)
  shl2atm(m+1:m+k) = i
  m = m + k
 end do ! for i

 deallocate(elem1, atmnshl, angshl, bas_num)
 if(sph) then
  forall(i=1:ncontr, (shltyp(i)>1)) shltyp(i) = -shltyp(i)
 end if
end subroutine read_shltyp_prpshl_shl2atm_from_amo

! find two variables LenNCZ and necpatm in a specified .amo file
subroutine find_LenNCZ_in_amo(amoname, LenNCZ, necpatm)
 implicit none
 integer :: fid
 integer, intent(out) :: LenNCZ, necpatm
 integer, allocatable :: ecpnum(:,:)
 character(len=240) :: buf
 character(len=240), intent(in) :: amoname

 LenNCZ = 0; necpatm = 0
 open(newunit=fid,file=TRIM(amoname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:7) == 'NECPAtm') exit
 end do ! for while

 read(fid,*) necpatm
 allocate(ecpnum(6,necpatm), source=0)

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:6) == 'ECPNum') exit
 end do ! for while

 read(fid,*) ecpnum
 close(fid)
 LenNCZ = SUM(ecpnum)
 deallocate(ecpnum)
end subroutine find_LenNCZ_in_amo

subroutine read_ecp_arrays_from_amo(amoname, natom, LenNCZ, necpatm, KFirst, &
                                    KLast, Lmax, LPSkip, NLP, CLP, ZLP)
 implicit none
 integer :: i, j, k, m, fid
 integer, intent(in) :: natom, LenNCZ, necpatm
 integer, intent(out) :: KFirst(natom,10), KLast(natom,10), Lmax(natom), &
  LPSkip(natom), NLP(LenNCZ)
 integer, allocatable :: ecpidx(:), ecpblock(:), ecpnum(:,:), ecppow(:,:,:)
 real(kind=8), intent(out) :: CLP(LenNCZ), ZLP(LenNCZ)
 real(kind=8), allocatable :: ecpexp(:,:,:), ecpcof(:,:,:)
 character(len=240) :: buf
 character(len=240), intent(in) :: amoname
 logical, allocatable :: lecpatm(:)

 KFirst = 0; KLast = 0; Lmax = 0; LPSkip = 1; NLP = 0; CLP = 0d0; ZLP = 0d0
 allocate(lecpatm(natom),ecpidx(necpatm), ecpblock(necpatm), ecpnum(6,necpatm),&
          ecppow(10,6,necpatm), ecpexp(10,6,necpatm), ecpcof(10,6,necpatm))

 open(newunit=fid,file=TRIM(amoname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:7) == 'lECPAtm') exit
 end do ! for while
 read(fid,*) lecpatm

 read(fid,'(A)') buf ! ECPBlock:
 read(fid,*) ecpblock

 read(fid,'(A)') buf ! ECPNum:
 read(fid,*) ecpnum

 read(fid,'(A)') buf ! ECPPow:
 read(fid,*) ecppow

 read(fid,'(A)') buf ! ECPExp:
 read(fid,*) ecpexp

 read(fid,'(A)') buf ! ECPCof:
 read(fid,*) ecpcof
 close(fid)

 forall(i=1:natom, lecpatm(i)) LPSkip(i) = 0

 j = 0
 do i = 1, natom, 1
  if(lecpatm(i)) then
   j = j + 1
   Lmax(i) = ecpblock(j) - 1
   ecpidx(j) = i
  end if
 end do ! for i
 deallocate(lecpatm, ecpblock)

 m = 0
 do i = 1, necpatm, 1
  do j = 1, 6
   k = ecpnum(j,i)
   if(k == 0) cycle
   NLP(m+1:m+k) = ecppow(1:k,j,i)
   ZLP(m+1:m+k) = ecpexp(1:k,j,i)
   CLP(m+1:m+k) = ecpcof(1:k,j,i)
   KFirst(ecpidx(i),j) = m + 1
   KLast(ecpidx(i),j) = m + k
   m = m + k
  end do ! for j
 end do ! for i

 deallocate(ecppow, ecpexp, ecpcof, ecpnum)
end subroutine read_ecp_arrays_from_amo

