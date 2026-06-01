
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

! read total charge and spin multiplicity from a given .gjf file
subroutine read_charge_and_mult_from_gjf(gjfname, charge, mult)
 implicit none
 integer :: i, nblank, fid
 integer, intent(out) :: charge, mult
!f2py intent(out) :: charge, mult
 character(len=240) :: buf
 character(len=240), intent(in) :: gjfname
!f2py intent(in) :: gjfname

 charge = 0; mult = 1
 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
 nblank = 0

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) nblank = nblank + 1
  if(nblank == 2) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_charge_and_mult_from_gjf: problema&
                   &tic file '//TRIM(gjfname)
  stop
  close(fid)
 end if

 read(fid,*,iostat=i) charge, mult
 close(fid)

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_charge_and_mult_from_gjf: problema&
                   &tic charge or'
  write(6,'(A)') 'spin multiplicity in file '//TRIM(gjfname)
  stop
 end if
end subroutine read_charge_and_mult_from_gjf

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

! read the total charge and the spin mltiplicity from a given .mkl file
subroutine read_charge_and_mult_from_mkl(mklname, charge, mult)
 implicit none
 integer :: i, fid
 integer, intent(out) :: charge, mult
!f2py intent(out) :: charge, mult
 character(len=240) :: buf
 character(len=240), intent(in) :: mklname
!f2py intent(in) :: mklname

 charge = 0; mult = 1
 call require_file_exist(mklname)
 open(newunit=fid,file=TRIM(mklname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:7) == '$CHAR_M') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_charge_and_mult_from_mkl: no '$CHA&
                   &R_M' found in"
  write(6,'(A)') 'file '//TRIM(mklname)
  close(fid)
  stop
 end if

 read(fid,*) charge, mult
 close(fid)
end subroutine read_charge_and_mult_from_mkl

! read the charge and spin multiplicity from a specified CP2K .inp file
subroutine read_charge_and_mult_from_cp2k_inp(inpname, charge, mult)
 implicit none
 integer :: i, k, fid
 integer, intent(out) :: charge, mult
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname
 logical :: found_charge, found_mult

 charge = 0; mult = 1 ! default value
 found_charge = .false.; found_mult = .false.
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  buf = ADJUSTL(buf)
  k = LEN_TRIM(buf)
  call upper(buf(1:k))
  if(buf(1:6) == 'CHARGE') then
   read(buf(7:),*) charge
   found_charge = .true.
  end if
  if(buf(1:12) == 'MULTIPLICITY') then
   read(buf(i+13:),*) mult
   found_mult = .true.
  end if
  if(found_charge .and. found_mult) exit
 end do ! for while

 close(fid)
end subroutine read_charge_and_mult_from_cp2k_inp

! read spin multiplicity from a specified ORCA input file
subroutine read_mult_from_orca_inp(inpname, mult)
 implicit none
 integer :: i, j, itype, fid
 integer, intent(out) :: mult
!f2py intent(out) :: mult
 character(len=45), parameter :: error_warn = 'ERROR in subroutine read_mult_fr&
                                              &om_orca_inp: '
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
  write(6,'(/,A)') error_warn//'spin multiplicity cannot be'
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
   write(6,'(/,A,I0)') error_warn//'the spin multiplicity cannot be'
   write(6,'(A)') 'found in file '//TRIM(inpname)
   stop
  end if
  i = INDEX(buf, '=')
  read(buf(i+1:),*) mult
 case default
  write(6,'(/,A,I0)') error_warn//'invalid itype=', itype
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
  write(6,'(/,A)') 'ERROR in subroutine read_mult_from_gms_gms: "SPIN MULT" not&
                   & found in file'
  write(6,'(A)') TRIM(gmsname)
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
 character(len=45), parameter :: error_warn = 'ERROR in subroutine read_npair_f&
                                              &rom_uno_out: '
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
  write(6,'(/,A)') error_warn//'nbf<nif. This is impossible.'
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
  write(6,'(/,A)') error_warn//'"ndb" not found.'
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
  write(6,'(/,A)') error_warn//'"idx" not found.'
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
 character(len=58), parameter :: error_warn = 'ERROR in subroutine average_pnmr&
                                              &_shield_in_orca_pnmr_out: '
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
  write(6,'(/,A)') error_warn//'"Paramagnetic shielding"'
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
  write(6,'(/,A)') error_warn//'the pNMR isotropic'
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
 character(len=52), parameter :: error_warn = 'ERROR in subroutine average_nmr_&
                                              &shield_in_orca_out: '
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
  write(6,'(/,A)') error_warn//'"CHEMICAL SHIELDING SUM"'
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
  write(6,'(/,A)') error_warn//'the NMR isotropic sh-'
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
 character(len=51), parameter :: error_warn = 'ERROR in subroutine average_nmr_&
                                              &shield_in_gau_log: '
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
  write(6,'(/,A)') error_warn//'"SCF GIAO Mag" not'
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
  write(6,'(/,A)') error_warn//'the NMR isotropic sh-'
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

! read the number of frames from xyz file
subroutine read_nframe_from_xyz(xyzname, nframe)
 implicit none
 integer :: i, fid, natom
 integer, intent(out) :: nframe
!f2py intent(out) :: nframe
 character(len=240) :: buf
 character(len=240), intent(in) :: xyzname
!f2py intent(in) :: xyzname

 nframe = 0
 open(newunit=fid,file=TRIM(xyzname),status='old',position='rewind')

 do while(.true.)
  read(fid,*,iostat=i) natom
  if(i /= 0) exit
  do i = 1, natom+1, 1
   read(fid,'(A)') buf
  end do ! for i
  nframe = nframe + 1
 end do ! for while

 close(fid)
end subroutine read_nframe_from_xyz

! read the number of frames from pdb file
subroutine read_nframe_from_pdb(pdbname, nframe)
 implicit none
 integer :: i, fid
 integer, intent(out) :: nframe
!f2py intent(out) :: nframe
 character(len=240) :: buf
 character(len=240), intent(in) :: pdbname
!f2py intent(in) :: pdbname

 nframe = 1 ! initialization
 open(newunit=fid,file=TRIM(pdbname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(1:5) == 'MODEL') exit
 end do ! for while

 close(fid)
 if(i /= 0) return ! assume 1 frame

 i = INDEX(buf, ' ')
 read(buf(i+1:),*) nframe
end subroutine read_nframe_from_pdb

! read elements and Cartesian coordinates from a Dalton .mol file
subroutine read_elem_and_coor_from_dalton_mol(molname, natom, elem, coor, nline)
 implicit none
 integer :: i, fid
 integer, intent(in) :: natom
 integer, intent(out) :: nline(natom) ! No. of lines of basis set data per atom
 real(kind=8), intent(out) :: coor(3,natom)
 character(len=2), intent(out) :: elem(natom)
 character(len=240) :: buf
 character(len=240), intent(in) :: molname

 elem = ' '; coor = 0d0; nline = 0
 open(newunit=fid,file=TRIM(molname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:10) == 'AtomTypes=') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_elem_and_coor_from_dalton_mol: no &
                   &'AtomTypes=' found"
  write(6,'(A)') 'in file '//TRIM(molname)
  close(fid)
  stop
 end if

 read(buf(11:),*) i
 if(i /= natom) then
  close(fid)
  write(6,'(/,A)') 'ERROR in subroutine read_elem_and_coor_from_dalton_mol: No.&
                   & atoms in .mol is'
  write(6,'(A)') 'not equal to input natom.'
  stop
 end if
 read(fid,'(A)') buf

 do i = 1, natom, 1
  read(fid,*) elem(i), coor(:,i)
  do while(.true.)
   read(fid,'(A)') buf
   if(i == natom) then
    if(LEN_TRIM(buf) == 0) exit
   end if
   if(buf(1:7) == 'Charge=') exit
   nline(i) = nline(i) + 1
  end do ! for while
 end do ! for i

 close(fid)
end subroutine read_elem_and_coor_from_dalton_mol

subroutine read_elem_and_coor_from_cp2k_inp(inpname, natom, elem, coor)
 implicit none
 integer :: i, k, fid
 integer, intent(in) :: natom
!f2py intent(in) :: natom
 real(kind=8), intent(out) :: coor(3,natom)
!f2py intent(out) :: coor
!f2py depend(natom) :: coor
 character(len=2), intent(out) :: elem(natom)
!f2py intent(out) :: elem
!f2py depend(natom) :: elem
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname
!f2py intent(in) :: inpname

 elem = '  '; coor = 0d0
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  buf = ADJUSTL(buf)
  k = LEN_TRIM(buf)
  call upper(buf(1:k))
  if(buf(1:6) == '&COORD') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_elem_and_coor_from_cp2k_inp: faile&
                   &d to locate `&COORD`'
  write(6,'(A)') 'in file '//TRIM(inpname)
  close(fid)
  stop
 end if

 do i = 1, natom, 1
  read(fid,*) elem(i), coor(:,i)
 end do ! for i

 close(fid)
end subroutine read_elem_and_coor_from_cp2k_inp

! write/create a Gaussian .EOu file
subroutine write_EOu(EOu, e, natom, grad)
 implicit none
 integer :: fid
 integer, intent(in) :: natom
 real(kind=8), intent(in) :: e, grad(3*natom)
 character(len=720), intent(in) :: EOu

 open(newunit=fid,file=TRIM(EOu),status='replace')
 write(fid,'(4D20.12)') e, 0d0,0d0,0d0
 write(fid,'(3D20.12)') grad
 close(fid)
end subroutine write_EOu

! write/create a .gjf file
subroutine write_gjf(gjfname, charge, mult, natom, elem, coor)
 implicit none
 integer :: i, fid
 integer, intent(in) :: charge, mult, natom
!f2py intent(in) :: charge, mult, natom
 real(kind=8), intent(in) :: coor(3,natom)
!f2py intent(in) :: coor
!f2py depend(natom) :: coor
 character(len=240) :: chkname
 character(len=2), intent(in) :: elem(natom)
!f2py intent(in) :: elem
!f2py depend(natom) :: elem
 character(len=240), intent(in) :: gjfname
!f2py intent(in) :: gjfname

 call find_specified_suffix(gjfname, '.gjf', i)
 chkname = gjfname(1:i-1)//'.chk'

 open(newunit=fid,file=TRIM(gjfname),status='replace')
 write(fid,'(A)') '%chk='//TRIM(chkname)
 write(fid,'(A)') '%nprocshared=1'
 write(fid,'(A)') '%mem=2GB'
 write(fid,'(A)') '#p B3LYP/6-31G(d,p) em=GD3BJ nosymm int=nobasistransform'
 write(fid,'(/,A,/)') 'Title'
 write(fid,'(I0,1X,I0)') charge, mult

 do i = 1, natom, 1
  write(fid,'(A2,3(1X,F18.8))') elem(i), coor(:,i)
 end do ! for i

 write(fid,'(/)',advance='no')
 close(fid)
end subroutine write_gjf

! write a frame of molecule into a given .pdb file
subroutine write_frame_into_pdb(pdbname, iframe, natom, cell, elem, resname, &
                                coor, append)
 implicit none
 integer :: i, fid
 integer, intent(in) :: iframe, natom
!f2py intent(in) :: iframe, natom
 character(len=2), intent(in) :: elem(natom)
!f2py intent(in) :: elem
!f2py depend(natom) :: elem
 character(len=3), intent(in) :: resname(natom)
!f2py intent(in) :: resname
!f2py depend(natom) :: resname
 character(len=240), intent(in) :: pdbname
!f2py intent(in) :: pdbname
 real(kind=8), intent(in) :: cell(6), coor(3,natom)
!f2py intent(in) :: cell, coor
!f2py depend(natom) :: coor
 logical, intent(in) :: append
!f2py intent(in) :: append

 if(append) then
  open(newunit=fid,file=TRIM(pdbname),status='old',position='append')
 else
  open(newunit=fid,file=TRIM(pdbname),status='replace')
 end if

 write(fid,'(A)') 'REMARK   1 File created by rwgeom of MOKIT'
 if(ANY(cell > 1d-4)) then
  write(fid,'(A,3(1X,F8.3),3(1X,F6.2),A)') 'CRYST1',cell(1:3),cell(4:6),' P 1           1'
 end if
 if(iframe > 0) write(fid,'(A,1X,I8)') 'MODEL',iframe

 do i = 1, natom, 1
  if(LEN_TRIM(resname(i)) == 0) then
   write(fid,'(A6,I5,2X,A2,10X,I1,4X,3F8.3,A)') 'HETATM', i, elem(i), 0, &
    coor(1:3,i),'  1.00  0.00'
  else
   write(fid,'(A4,I7,2X,A2,2X,A3,5X,I1,4X,3F8.3,A)') 'ATOM', i, elem(i), &
    resname(i), 0, coor(1:3,i), '  1.00  0.00'
  end if
 end do ! for i

 write(fid,'(A)') 'END'
 close(fid)
end subroutine write_frame_into_pdb

