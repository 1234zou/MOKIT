
! Read nbf and nif from a GAMESS .dat file. The number of basis functions in
! GAMESS .dat is always at Cartesian basis, no matter the actual calculation
! is performed at Cartesian or spherical harmonic basis.
subroutine read_cart_nbf_nif_from_dat(datname, only_nbf, nbf, nif)
 implicit none
 integer :: i, j, k, fid
 integer, intent(out) :: nbf, nif
!f2py intent(out) :: nbf, nif
 character(len=240) :: buf
 character(len=240), intent(in) :: datname
!f2py character(len=240), intent(in) :: datname
 logical, intent(in) :: only_nbf
!f2py intent(in) :: only_nbf

 nbf = 0; nif = 0
 open(newunit=fid,file=TRIM(datname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:2) == '$') then
   call upper(buf(3:5))
   if(buf(2:5) == '$VEC') exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cart_nbf_nif_from_dat: no '$VEC' f&
                   &ound in file "//TRIM(datname)
  close(fid)
  stop
 end if

 j = 0
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:2) /= ' 1') exit
  j = j + 1
 end do ! for while

 k = j ! backup

 BACKSPACE(fid)
 BACKSPACE(fid)
 read(fid,'(A)') buf
 nbf = (j-1)*5 + (LEN_TRIM(buf)-5)/15

 if(only_nbf) then ! no need to find nif
  close(fid)
  return
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:2) == '$') then
   call upper(buf(3:5))
   if(buf(2:5) == '$END') exit
  end if
  j = j + 1
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cart_nbf_nif_from_dat: no '$END' c&
                   &orresponds to '$VEC'"
  write(6,'(A)') 'in file '//TRIM(datname)
  stop
 end if

 nif = j/k
end subroutine read_cart_nbf_nif_from_dat

! read nbf and nif from .Orb file of MOLCAS/OpenMOLCAS
subroutine read_nbf_and_nif_from_orb(orbname, nbf, nif)
 implicit none
 integer :: fid
 integer, intent(out) :: nbf, nif
 character(len=240) :: buf
 character(len=240), intent(in) :: orbname

 open(newunit=fid,file=TRIM(orbname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:5) == '#INFO') exit
 end do

 read(fid,'(A)') buf
 read(fid,'(A)') buf
 read(fid,*) nbf
 read(fid,*) nif
 close(fid)
end subroutine read_nbf_and_nif_from_orb

! read the total charge and the spin mltiplicity from a given .fch(k) file
subroutine read_charge_and_mult_from_fch(fchname, charge, mult)
 implicit none
 integer :: i, fid
 integer, intent(out) :: charge, mult
!fp2y intent(out) :: charge, mult
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 charge = 0; mult = 1
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:6) == 'Charge') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_charge_and_mult_from_fch: no 'Char&
                   &ge' found in"
  write(6,'(A)') 'file '//TRIM(fchname)
  close(fid)
  stop
 end if

 BACKSPACE(fid)
 read(fid,'(A49,2X,I10)') buf, charge
 read(fid,'(A49,2X,I10)') buf, mult
 close(fid)
end subroutine read_charge_and_mult_from_fch

! read spin multiplicity from a specified ORCA input file
subroutine read_mult_from_orca_inp(inpname, mult)
 implicit none
 integer :: i, j, itype, fid
 integer, intent(out) :: mult
!f2py intent(out) :: mult
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname
!f2py intent(in) :: inpname

 mult = 1; itype = 1
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:1) == '*') exit       ! * xyz 0 1
  if(buf(1:7) == '%coords') then ! Mult = 1
   itype = 2
   exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_mult_from_orca_inp: spin multiplic&
                   &ity cannot be'
  write(6,'(A)') 'found in file '//TRIM(inpname)
  close(fid)
  stop
 end if

 select case(itype)
 case(1)
  close(fid)
  i = LEN_TRIM(buf)
  j = INDEX(buf(1:i), ' ', back=.true.)
  read(buf(j+1:i),*) mult
 case(2)
  do i = 1, 4
   read(fid,'(A)') buf
   if(INDEX(buf(1:5),'Mult') > 0) exit
  end do ! for i
  close(fid)
  if(i == 5) then
   write(6,'(/,A,I0)') 'ERROR in subroutine read_mult_from_orca_inp: spin multi&
                       &plicity cannot be'
   write(6,'(A)') 'found in file '//TRIM(inpname)
   stop
  end if
  i = INDEX(buf, '=')
  read(buf(i+1:),*) mult
 case default
  write(6,'(/,A,I0)') 'ERROR in subroutine read_mult_from_orca_inp: invalid ity&
                      &pe=', itype
  write(6,'(A)') 'inpname='//TRIM(inpname)
  close(fid)
  stop
 end select
end subroutine read_mult_from_orca_inp

! read spin multiplicity from a specified GAMESS output file (.gms)
subroutine read_mult_from_gms_gms(gmsname, mult)
 implicit none
 integer :: i, fid
 integer, intent(out) :: mult
!f2py intent(out) :: mult
 character(len=240) :: buf
 character(len=240), intent(in) :: gmsname
!f2py intent(in) :: gmsname

 mult = 1
 open(newunit=fid,file=TRIM(gmsname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:10) == 'SPIN MULT') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_mult_from_gms_gms: 'SPIN MULT' not&
                   & found in file "//TRIM(gmsname)
  stop
 end if

 i = INDEX(buf, '=')
 read(buf(i+1:),*) mult
end subroutine read_mult_from_gms_gms

! check whether pure Cartesian functions
subroutine check_cart_compatibility_in_fch(fchname, cart)
 implicit none
 integer :: icart
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname
 logical, intent(in) :: cart
!f2py intent(in) :: cart

 call find_icart_in_fch(fchname, .false., icart)

 if(cart .and. icart==1) then
  write(6,'(/,A)') 'ERROR in subroutine check_cart_compatibility_in_fch: Cartes&
                   &ian functions are'
  write(6,'(A)') 'requested by the user. But you provided a .fch(k) file which &
                 &uses spherical harmonic'
  write(6,'(A)') "functions. Two possible solutions: 1) delete the keyword 'Car&
                 &t' in MOKIT{}"
  write(6,'(A)') '2) provide another .fch file which uses pure Cartesian functi&
                 &ons.'
  write(6,'(A)') 'fchname='//TRIM(fchname)
  stop
 end if

 if((.not.cart) .and. icart==2) then
  write(6,'(/,A)') 'ERROR in subroutine check_cart_compatibility_in_fch: spheri&
                   &cal harmonic functions'
  write(6,'(A)') 'are set as default. But you provided a .fch(k) file which has&
                 & Cartesian functions.'
  write(6,'(A)') "Two possible solutions: 1) add the keyword 'Cart' in MOKIT{};&
                 & 2) provide another .fch"
  write(6,'(A)') 'file which uses pure spherical harmonic functions.'
  write(6,'(A)') 'fchname='//TRIM(fchname)
  stop
 end if
end subroutine check_cart_compatibility_in_fch

! check whether UHF is used in a specified CFOUR output file
subroutine check_uhf_in_cfour_out(outname, uhf)
 implicit none
 integer :: i, k, fid
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
!f2py intent(in) :: outname
 logical, intent(out) :: uhf
!f2py intent(out) :: uhf

 uhf = .false.
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(9:23) == 'SCF reference f') then
   k = INDEX(buf, ':')
   buf = ADJUSTL(buf(k+1:))
   if(buf(1:3) == 'UHF') uhf = .true.
   exit
  end if
 end do ! for while

 close(fid)
end subroutine check_uhf_in_cfour_out

subroutine read_npair_from_uno_out(unofile,nbf,nif,ndb,npair,nopen,lin_dep)
 implicit none
 integer :: i, fid, idx(3), nvir
 integer, intent(out) :: nbf, nif, ndb, npair, nopen
!f2py intent(out) :: nbf, nif, ndb, npair, nopen
 character(len=240), intent(in) :: unofile
!f2py intent(in) :: unofile
 character(len=240) :: buf
 logical, intent(out) :: lin_dep
!f2py intent(out) :: lin_dep

 nbf = 0; nif = 0; ndb = 0; npair = 0; nopen = 0
 lin_dep = .false.; buf = ' '
 open(newunit=fid,file=TRIM(unofile),status='old',position='rewind')

 read(fid,'(A)') buf
 call get_int_after_flag(buf, '=', .true., nbf)

 read(fid,'(A)') buf
 call get_int_after_flag(buf, '=', .true., nif)

 if(nbf > nif) then
  lin_dep = .true.
 else if(nbf < nif) then
  write(6,'(/,A)') 'ERROR in subroutine read_npair_from_uno_out: nbf<nif. This &
                   &is impossible.'
  write(6,'(A)') 'Please check unofile: '//TRIM(unofile)
  close(fid)
  stop
 end if

 close(fid)
 open(newunit=fid,file=TRIM(unofile),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:3) == 'ndb') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_npair_from_uno_out: 'ndb' not found."
  write(6,'(A)') 'Please open file '//TRIM(unofile)//' and check.'
  close(fid)
  stop
 end if
 call get_int_after_flag(buf, '=', .true., ndb)

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:3) == 'idx') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_npair_from_uno_out: 'idx' not found."
  write(6,'(A)') 'Please open file '//TRIM(unofile)//' and check.'
  close(fid)
  stop
 end if

 idx = 0
 i = INDEX(buf,'=')
 read(buf(i+1:),*) idx(1:3)
 close(fid,status='delete')

 npair = (idx(2) - idx(1) - idx(3))/2
 nvir = nif - ndb - 2*npair - idx(3)
 nopen = idx(3)
 write(6,'(6(A,I0))') 'nbf=', nbf, ', nif=', nif, ', doubly_occ=', idx(1)-1, &
                      ', npair=', npair, ', nopen=', idx(3), ', nvir=', nvir
end subroutine read_npair_from_uno_out

! find npair0: the number of active pairs (|C2| > 0.1)
! (assuming that the pair coefficients haven been sorted)
subroutine find_npair0_from_dat(datname, npair, npair0)
 implicit none
 integer :: i, k, datid
 integer, intent(in) :: npair
!f2py intent(in) :: npair
 integer, intent(out) :: npair0
!f2py intent(out) :: npair0
 real(kind=8) :: rtmp
 real(kind=8), allocatable :: pair_coeff(:,:)
 character(len=240) :: buf
 character(len=240), intent(in) :: datname
!f2py intent(in) :: datname

 npair0 = 0 ! initialization
 if(npair == 0) then
  write(6,'(/,A)') 'Warning in subroutine find_npair0_from_dat: npair=npair0=0.'
  return
 end if

 ! find pair coefficients
 open(newunit=datid,file=TRIM(datname),status='old',position='rewind')
 do while(.true.)
  read(datid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(INDEX(buf,'CICOEF(') /= 0) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine find_npair0_from_dat: 'CICOEF(' keyword&
                   & not found in"
  write(6,'(A)') 'file '//TRIM(datname)
  close(datid)
  stop
 end if

 ! read pair coefficients
 BACKSPACE(datid)
 allocate(pair_coeff(2,npair), source=0d0)
 do i = 1, npair, 1
  read(datid,'(A)') buf
  k = INDEX(buf,'=')
  read(buf(k+1:),*) pair_coeff(1,i)
  k = INDEX(buf,',')
  read(buf(k+1:),*) pair_coeff(2,i)
 end do ! for i
 ! pair coefficients read done

 close(datid)
 npair0 = COUNT(pair_coeff(2,:) <= -0.1d0)

 do i = 1, npair, 1
  rtmp = pair_coeff(2,i) + 0.1d0
  if(rtmp>0d0 .and. rtmp<2.6d-3) then
   write(6,'(/,A)') REPEAT('-',79)
   write(6,'(A)') 'Warning from subroutine find_npair0_from_dat: some occupatio&
                  &n number is very'
   write(6,'(A)') 'close to ON_thres. You may consider enlarge the active space&
                  & size slightly in'
   write(6,'(A,I0)') 'later CAS calculation. Check No. pair: ', i
   write(6,'(A)') REPEAT('-',79)
  end if
 end do ! for i

 deallocate(pair_coeff)
end subroutine find_npair0_from_dat

! read ncore, nopen and npair from a GAMESS .gms file
subroutine read_npair_from_gms_gms(gmsname, ncore, nopen, npair)
 implicit none
 integer :: i, fid
 integer, intent(out) :: ncore, nopen, npair
!f2py intent(out) :: ncore, nopen, npair
 character(len=240) :: buf
 character(len=240), intent(in) :: gmsname
!f2py intent(in) :: gmsname

 ncore = 0; nopen = 0; npair = 0; buf = ' '
 open(newunit=fid,file=TRIM(gmsname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)', iostat=i) buf
  if(i /= 0) exit
  if(buf(34:41) == 'NCO    =') exit
 end do

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_npair_from_gms_gms: no 'NCO    =' &
                   &found in file"
  write(6,'(A)') TRIM(gmsname)
  close(fid)
  stop
 end if

 read(buf(42:),*) ncore
 read(fid,'(A)') buf
 read(buf(19:),*) npair
 read(buf(42:),*) nopen
 close(fid)
end subroutine read_npair_from_gms_gms

! find the target CASCI root in a specified PySCF CASCI output file
subroutine read_target_root_from_pyscf_out(outname, target_root, found)
 implicit none
 integer :: i, fid
 integer, intent(out) :: target_root
!f2py intent(out) :: target_root
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
!f2py intent(in) :: outname
 logical, intent(out) :: found
!f2py intent(out) :: found

 found = .false.; target_root = 0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:12) == 'target_root=') then
   read(buf(13:),*) target_root
   found = .true.
   exit
  end if
 end do ! for while

 close(fid)
end subroutine read_target_root_from_pyscf_out

! find pNMR isotropic shieldings of target atoms in ORCA pNMR output, and
! calculate the average value
subroutine average_pnmr_shield_in_orca_pnmr_out(pnmr_out, natom, atom_list, &
                                                ave_val)
 implicit none
 integer :: i, iatom, fid
 integer, intent(in) :: natom
!f2py intent(in) :: natom
 integer, intent(in) :: atom_list(natom)
!f2py intent(in) :: atom_list
!f2py depend(natom) :: atom_list
 real(kind=8), intent(out) :: ave_val ! in ppm
!f2py intent(out) :: ave_val
 real(kind=8) :: iso_shield
 real(kind=8), allocatable :: pnmr_shielding(:)
 character(len=2) :: elem
 character(len=11) :: iatom_elem
 character(len=240), intent(in) :: pnmr_out
!f2py intent(in) :: pnmr_out
 character(len=240) :: buf
 logical, allocatable :: found(:) ! size natom

 ave_val = 0d0
 buf = pnmr_out ! buf has length declared as 240
 call require_file_exist(buf)
 open(newunit=fid,file=TRIM(pnmr_out),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:22) == 'Paramagnetic shielding') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine average_pnmr_shield_in_orca_pnmr_out: '&
                   &Paramagnetic shielding'"
  write(6,'(A)') 'not located in file '//TRIM(pnmr_out)
  close(fid)
  stop
 end if

 do i = 1, 5
  read(fid,'(A)') buf
 end do

 allocate(found(natom), pnmr_shielding(natom))
 found = .false.; pnmr_shielding = 0d0

 do while(.true.)
  read(fid,'(A)') buf
  read(fid,'(A)') buf
  read(buf,*) iatom_elem, iso_shield
  call split_iatom_elem(iatom_elem, iatom, elem)
  iatom = iatom + 1
  do i = 1, natom, 1
   if(atom_list(i) == iatom) exit
  end do ! for i
  if(i < natom+1) then
   found(i) = .true.
   pnmr_shielding(i) = iso_shield
  end if
  if(ALL(found .eqv. .true.)) exit
  read(fid,'(A)') buf
  read(fid,'(A)') buf
  if(buf(1:5) == '-----') exit
 end do ! for while

 close(fid)
 if(ANY(found .eqv. .false.)) then
  do i = 1, natom, 1
   if(.not. found(i)) exit
  end do ! for i
  write(6,'(/,A)') 'ERROR in subroutine average_pnmr_shield_in_orca_pnmr_out: t&
                   &he pNMR isotropic'
  write(6,'(A,I0,A)') 'shielding of atom label ',atom_list(i),' is not found in&
                      & file '//TRIM(pnmr_out)
  deallocate(found, pnmr_shielding)
  stop
 end if

 deallocate(found)
 if(natom == 1) then
  ave_val = pnmr_shielding(1)
 else
  ave_val = SUM(pnmr_shielding)/DBLE(natom)
 end if

 deallocate(pnmr_shielding)
end subroutine average_pnmr_shield_in_orca_pnmr_out

! find NMR isotropic shieldings of target atoms in ORCA NMR output, and
! calculate the average value
subroutine average_nmr_shield_in_orca_out(outname, natom, atom_list, ave_val)
 implicit none
 integer :: i, iatom, fid
 integer, intent(in) :: natom
!f2py intent(in) :: natom
 integer, intent(in) :: atom_list(natom)
!f2py intent(in) :: atom_list
!f2py depend(natom) :: atom_list
 real(kind=8), intent(out) :: ave_val
!f2py intent(out) :: ave_val
 real(kind=8) :: iso_shield
 real(kind=8), allocatable :: nmr_shielding(:)
 character(len=2) :: elem
 character(len=240), intent(in) :: outname
!f2py intent(in) :: outname
 character(len=240) :: buf
 logical, allocatable :: found(:) ! size natom

 ave_val = 0d0
 buf = outname ! buf has length declared as 240
 call require_file_exist(buf)
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:22) == 'CHEMICAL SHIELDING SUM') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine average_nmr_shield_in_orca_out: 'CHEMIC&
                   &AL SHIELDING SUM'"
  write(6,'(A)') 'not located in file '//TRIM(outname)
  close(fid)
  stop
 end if

 allocate(found(natom), nmr_shielding(natom))
 found = .false.; nmr_shielding = 0d0
 do i = 1, 5
  read(fid,'(A)') buf
 end do ! for i

 do while(.true.)
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
  read(buf,*) iatom, elem, iso_shield
  iatom = iatom + 1
  do i = 1, natom, 1
   if(atom_list(i) == iatom) exit
  end do ! for i
  if(i < natom+1) then
   found(i) = .true.
   nmr_shielding(i) = iso_shield
  end if
  if(ALL(found .eqv. .true.)) exit
 end do ! for while

 close(fid)
 if(ANY(found .eqv. .false.)) then
  do i = 1, natom, 1
   if(.not. found(i)) exit
  end do ! for i
  write(6,'(/,A)') 'ERROR in subroutine average_nmr_shield_in_orca_out: the NMR&
                   & isotropic sh-'
  write(6,'(A,I0,A)') 'ielding of atom label ',atom_list(i),' is not found in f&
                      &ile '//TRIM(outname)
  deallocate(found, nmr_shielding)
  stop
 end if

 deallocate(found)
 if(natom == 1) then
  ave_val = nmr_shielding(1)
 else
  ave_val = SUM(nmr_shielding)/DBLE(natom)
 end if

 deallocate(nmr_shielding)
end subroutine average_nmr_shield_in_orca_out

! find NMR isotropic shieldings of target atoms in Gaussian NMR output, and
! calculate the average value
subroutine average_nmr_shield_in_gau_log(logname, natom, atom_list, ave_val)
 implicit none
 integer :: i, j, iatom, fid
 integer, intent(in) :: natom
!f2py intent(in) :: natom
 integer, intent(in) :: atom_list(natom)
!f2py intent(in) :: atom_list
!f2py depend(natom) :: atom_list
 real(kind=8), intent(out) :: ave_val
!f2py intent(out) :: ave_val
 real(kind=8), allocatable :: nmr_shielding(:)
 character(len=240), intent(in) :: logname
!f2py intent(in) :: logname
 character(len=240) :: buf
 logical, allocatable :: found(:) ! size natom

 ave_val = 0d0
 buf = logname ! buf has length declared as 240
 call require_file_exist(buf)
 open(newunit=fid,file=TRIM(logname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:13) == 'SCF GIAO Mag') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine average_nmr_shield_in_gau_log: 'SCF GIA&
                   &O Mag' not"
  write(6,'(A)') 'located in file '//TRIM(logname)
  close(fid)
  stop
 end if

 allocate(found(natom), nmr_shielding(natom))
 found = .false.; nmr_shielding = 0d0

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:11)=='g value of' .or. buf(2:11)=='End of Min') exit
  read(buf,*) iatom
  do i = 1, natom, 1
   if(atom_list(i) == iatom) exit
  end do ! for i
  if(i < natom+1) then
   found(i) = .true.
   j = INDEX(buf, '=')
   read(buf(j+1:),*) nmr_shielding(i)
  end if
  if(ALL(found .eqv. .true.)) exit
  do i = 1, 4
   read(fid,'(A)') buf
  end do ! for i
 end do ! for while

 close(fid)
 if(ANY(found .eqv. .false.)) then
  do i = 1, natom, 1
   if(.not. found(i)) exit
  end do ! for i
  write(6,'(/,A)') 'ERROR in subroutine average_nmr_shield_in_gau_log: the NMR &
                   &isotropic sh-'
  write(6,'(A,I0,A)') 'ielding of atom label ',atom_list(i),' is not found in f&
                      &ile '//TRIM(logname)
  deallocate(found, nmr_shielding)
  stop
 end if

 deallocate(found)
 if(natom == 1) then
  ave_val = nmr_shielding(1)
 else
  ave_val = SUM(nmr_shielding)/DBLE(natom)
 end if

 deallocate(nmr_shielding)
end subroutine average_nmr_shield_in_gau_log

