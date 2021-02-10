! written by jxzou at 20200623

! find the number of atoms in Gaussian .gjf file
subroutine read_natom_from_gjf(gjfname, natom)
 implicit none
 integer :: i, fid, nblank
 integer, intent(out) :: natom
 integer, parameter :: iout = 6
 character(len=240) :: buf
 character(len=240), intent(in) :: gjfname

 nblank = 0
 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) nblank = nblank + 1
  if(nblank == 2) exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine read_natom_from_gjf: incomplete file '//TRIM(gjfname)
  stop
 end if

 read(fid,'(A)') buf ! skip charge and mult

 natom = 0
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) exit
  natom = natom + 1
 end do ! for while

 close(fid)
 return
end subroutine read_natom_from_gjf

! read the number of atoms from a given .fch file
subroutine read_natom_from_fch(fchname, natom)
 implicit none
 integer :: i, fid
 integer, intent(out) :: natom
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname

 natom = 0
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(buf(1:14) == 'Atomic numbers') exit
 end do ! for while

 BACKSPACE(fid)
 read(fid,'(A49,2X,I10)') buf, natom

 close(fid)
 return
end subroutine read_natom_from_fch

! read 3 arrays elem, nuc, coor, and the total charge as well as multiplicity
! from a given .gjf file
subroutine read_elem_and_coor_from_gjf(gjfname, natom, elem, nuc, coor, charge, mult)
 use fch_content, only: elem2nuc
 implicit none
 integer :: i, j, k, fid, nblank
 integer, intent(in) :: natom
 integer, intent(out) :: charge, mult, nuc(natom)
 integer, parameter :: iout = 6
 real(kind=8), intent(out) :: coor(3,natom)
 real(kind=8), parameter :: Bohr_const = 0.52917721092d0
 character(len=2), intent(out) :: elem(natom)
 character(len=240) :: buf
 character(len=240), intent(in) :: gjfname
 logical :: bohr

 charge = 0; mult = 1; coor = 0d0; bohr = .false.
 nblank = 0

 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) nblank = nblank + 1
  if(nblank == 2) exit

  if(buf(1:1) =='#') then
   if(index(buf,'units=au') /= 0) bohr = .true.
  end if
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine read_elem_and_coor_from_gjf: incomplete file '//TRIM(gjfname)
  stop
 end if

 read(fid,*,iostat=k) charge, mult
 if(k /= 0) then
  write(iout,'(A)') 'ERROR in subroutine read_elem_and_coor_from_gjf: failed&
                   & to read charge and mult.'
  write(iout,'(A)') 'There exists syntax error in file '//TRIM(gjfname)
  stop
 end if

 do i = 1, natom, 1
  read(fid,'(A)') buf
  j = index(buf,'('); k = index(buf,')') ! in case for the fragment guess wfn
  if(j*k /= 0) buf(j:k) = ' '

  read(buf,*,iostat=k) elem(i), coor(1:3,i)
  if(k /= 0) then
   write(iout,'(A)') 'ERROR in subroutine read_elem_and_coor_from_gjf: only 4-column&
                    & format is supported.'
   stop
  end if
 end do ! for i

 close(fid)

 if(bohr) coor = coor*Bohr_const ! convert Bohr to Angstrom

 forall(i=1:natom)
  elem(i) = ADJUSTL(elem(i))
  nuc(i) = elem2nuc(elem(i))
 end forall
 return
end subroutine read_elem_and_coor_from_gjf

! read charge, spin multiplicities and atom2frag from a given .gjf file
subroutine read_frag_guess_from_gjf(gjfname, natom, atom2frag, nfrag, frag_char_mult)
 implicit none
 integer :: i, j, k, nblank, charge, mult, fid
 integer, parameter :: iout = 6
 integer, intent(in) :: natom, nfrag
 integer, intent(out) :: atom2frag(natom), frag_char_mult(2,nfrag)
 character(len=240) :: buf
 character(len=240), intent(in) :: gjfname

 atom2frag = 0; frag_char_mult = 0
 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')

 nblank = 0
 do while(.true.)
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) nblank = nblank + 1
  if(nblank == 2) exit
 end do ! for while

 read(fid,*,iostat=i) charge, mult, ((frag_char_mult(j,i),j=1,2),i=1,nfrag)
 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine read_frag_guess_from_gjf: failed to read&
                   & charges and spin multiplicities of fragments.'
  write(iout,'(A)') 'Please check syntax in file '//TRIM(gjfname)
  close(fid)
  stop
 end if

 do i = 1, natom, 1
  read(fid,'(A)') buf
  j = index(buf,'='); k = index(buf,')')

  if(j*k == 0) then
   write(iout,'(A)') 'ERROR in subroutine read_frag_guess_from_gjf: failed to read&
                     & atom2frag.'
   write(iout,'(A)') 'Problematic line: '//TRIM(buf)
   write(iout,'(A)') 'Please check syntax in file '//TRIM(gjfname)
   close(fid)
   stop
  end if

  read(buf(j+1:k-1),*) atom2frag(i)
 end do ! for i

 close(fid)
 return
end subroutine read_frag_guess_from_gjf

! read 3 arrays elem, nuc, coor, and the total charge as well as multiplicity
! from a given .fch file
subroutine read_elem_and_coor_from_fch(fchname, natom, elem, nuc, coor, charge, mult)
 use fch_content, only: nuc2elem
 implicit none
 integer :: i, fid
 integer, intent(in) :: natom
 integer, intent(out) :: charge, mult, nuc(natom)
 integer, parameter :: iout = 6
 real(kind=8), intent(out) :: coor(3,natom)
 real(kind=8), allocatable :: coor0(:)
 real(kind=8), parameter :: Bohr_const = 0.52917721092d0
 character(len=2), intent(out) :: elem(natom)
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 ! find charge and mult
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:6) == 'Charge') exit
 end do
 BACKSPACE(fid)
 read(fid,'(A49,2X,I10)') buf, charge
 read(fid,'(A49,2X,I10)') buf, mult

 ! find atomic numbers/nuclear charges
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:14) == 'Atomic numbers') exit
 end do ! for i
 read(fid,'(6(1X,I11))') (nuc(i),i=1,natom)

 ! find and read coordinates
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:12) == 'Current cart') exit
 end do
 allocate(coor0(3*natom), source=0.0d0)
 read(fid,'(5(1X,ES15.8))') (coor0(i),i=1,3*natom)
 coor = RESHAPE(coor0,[3,natom])

 deallocate(coor0)
 close(fid)

 coor = coor*Bohr_const ! convert Bohr to Angstrom
 forall(i=1:natom) elem(i) = nuc2elem(nuc(i))
 return
end subroutine read_elem_and_coor_from_fch

! read Cartesian gradient from a given PySCF output file
subroutine read_grad_from_pyscf_out(outname, natom, grad)
 implicit none
 integer :: i, k, fid
 integer, intent(in) :: natom
 real(kind=8), intent(out) :: grad(3*natom)
 character(len=3) :: elem
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 grad = 0.0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(index(buf,'gradients') /= 0) exit
 end do ! for while

 read(fid,'(A)') buf

 do i = 1, natom, 1
  read(fid,*) k, elem, grad(3*i-2:3*i)
 end do ! for i

 close(fid)
 return
end subroutine read_grad_from_pyscf_out

! read Cartesian gradient from a given .fch file
subroutine read_grad_from_fch(fchname, natom, grad)
 implicit none
 integer :: i, fid
 integer, intent(in) :: natom
 real(kind=8), intent(out) :: grad(3*natom)
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname

 grad = 0.0d0
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:26) == 'Opt point       1 Gradient') exit
 end do ! for while

 read(fid,'(5(1X,ES15.8))') (grad(i),i=1,3*natom)
 close(fid)
 return
end subroutine read_grad_from_fch

! read Cartesian gradient from a given Gaussian .log file
subroutine read_grad_from_gau_log(logname, natom, grad)
 implicit none
 integer :: i, k1, k2, fid
 integer, intent(in) :: natom
 real(kind=8), intent(out) :: grad(3*natom)
 character(len=240) :: buf
 character(len=240), intent(in) :: logname

 grad = 0.0d0
 open(newunit=fid,file=TRIM(logname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(38:52) == 'Forces (Hartree') exit
 end do ! for while

 read(fid,'(A)') buf
 read(fid,'(A)') buf

 do i = 1, natom, 1
  read(fid,*) k1, k2, grad(3*i-2:3*i)
 end do ! for i

 close(fid)

 grad = -grad !!! VIP
 return
end subroutine read_grad_from_gau_log

! read Cartesian gradient from a given GAMESS .gms file
subroutine read_grad_from_gms_gms(outname, natom, grad)
 implicit none
 integer :: i, k, fid
 integer, intent(in) :: natom
 real(kind=8), intent(out) :: grad(3*natom)
 character(len=2) :: elem = ' '
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 grad = 0.0d0

 open(newunit=fid,file=TRIM(outname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(buf(26:47) == 'GRADIENT OF THE ENERGY') exit
 end do ! for while

 do i = 1, 3
  read(fid,'(A)') buf
 end do ! for i

 do i = 1, natom, 1
  read(fid,*) k, elem, grad(3*i-2:3*i)
 end do ! for i

 close(fid)
 return
end subroutine read_grad_from_gms_gms

! read Cartesian gradient from a given GAMESS .dat file
subroutine read_grad_from_gms_dat(datname, natom, grad)
 implicit none
 integer :: i, fid
 integer, intent(in) :: natom
 real(kind=4) :: r
 real(kind=8), intent(out) :: grad(3*natom)
 character(len=2) :: elem
 character(len=240) :: buf
 character(len=240), intent(in) :: datname

 grad = 0.0d0

 open(newunit=fid,file=TRIM(datname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(buf(2:6) == '$GRAD') exit
 end do ! for while

 read(fid,'(A)') buf
 do i = 1, natom, 1
  read(fid,*) elem, r, grad(3*i-2:3*i)
 end do ! for i

 close(fid)
 return
end subroutine read_grad_from_gms_dat

! read Cartesian gradient from a given MOLCAS/OpenMolcas .out file
subroutine read_grad_from_molcas_out(outname, natom, grad)
 implicit none
 integer :: i, fid
 integer, intent(in) :: natom
 real(kind=8), intent(out) :: grad(3*natom)
 character(len=2) :: elem = ' '
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 grad = 0.0d0

 open(newunit=fid,file=TRIM(outname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(buf(17:35) == 'Molecular gradients') exit
 end do ! for while

 do i = 1, 7, 1
  read(fid,'(A)') buf
 end do ! for i

 do i = 1, natom, 1
  read(fid,*) elem, grad(3*i-2:3*i)
 end do ! for i

 close(fid)
 return
end subroutine read_grad_from_molcas_out

! read Cartesian gradient from a given ORCA .out file
subroutine read_grad_from_orca_out(outname, natom, grad)
 implicit none
 integer :: i, k, fid
 integer, intent(in) :: natom
 real(kind=8), intent(out) :: grad(3*natom)
 character(len=1) :: str = ' '
 character(len=2) :: elem = ' '
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 grad = 0.0d0

 open(newunit=fid,file=TRIM(outname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(buf(1:14) == 'CARTESIAN GRAD') exit
 end do ! for while

 read(fid,'(A)') buf
 read(fid,'(A)') buf

 do i = 1, natom, 1
  read(fid,*) k, elem, str, grad(3*i-2:3*i)
 end do ! for i

 close(fid)
 return
end subroutine read_grad_from_orca_out

! read Cartesian gradient from a given Molpro .out file
subroutine read_grad_from_molpro_out(outname, natom, grad)
 implicit none
 integer :: i, k, fid
 integer, intent(in) :: natom
 real(kind=8), intent(out) :: grad(3*natom)
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 grad = 0.0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(buf(2:8) == 'MC GRAD') exit
 end do ! for while

 do i = 1, 3
  read(fid,'(A)') buf
 end do

 do i = 1, natom, 1
  read(fid,*) k, grad(3*i-2:3*i)
 end do ! for i

 close(fid)
 return
end subroutine read_grad_from_molpro_out

! read Cartesian gradient from a given BDF .out file
subroutine read_grad_from_bdf_out(outname, natom, grad)
 implicit none
 integer :: i, k, fid
 integer, intent(in) :: natom
 real(kind=8), intent(out) :: grad(3*natom)
 character(len=10) :: str
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 grad = 0.0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(index(buf,'Molecular gradient - Mol') /= 0) exit
 end do ! for while

 read(fid,'(A)') buf
 read(fid,'(A)') buf

 do i = 1, natom, 1
  read(fid,*) str, grad(3*i-2:3*i)
 end do ! for i

 close(fid)
 return
end subroutine read_grad_from_bdf_out

