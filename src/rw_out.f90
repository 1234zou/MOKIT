! read output files (.out/.log, etc.) to extract information other than geom, wfn
! also contains 'write' subroutines if some

! whether the Davidson Q CORRECTION can be found in a GAMESS output file
function has_davidson_q(outname) result(alive)
 implicit none
 integer :: i, fid
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical :: alive

 alive = .false.
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:18) == 'CALC. OF DAVIDSON') exit
 end do ! for while

 close(fid)
 if(i /= 0) return

 i = INDEX(buf, '=')
 read(buf(i+1:),*) alive
end function has_davidson_q

! find the target CASCI root in a specified PySCF CASCI output file
subroutine read_target_root_from_pyscf_out(outname, target_root, found)
 implicit none
 integer :: i, fid
 integer, intent(out) :: target_root
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(out) :: found

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