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

! find the number of atoms in GAMESS .inp file
subroutine read_natom_from_gms_inp(inpname, natom)
 implicit none
 integer :: i, fid, nline
 integer, intent(out) :: natom
 integer, parameter :: iout = 6
 character(len=1) :: str
 character(len=240):: buf
 character(len=240), intent(in) :: inpname

 natom = 0
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  call upper(buf(2:6))
  if(buf(2:6) == '$DATA') exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine read_natom_from_gms_inp: wrong format&
                   & in file '//TRIM(inpname)
  stop
 end if
 read(fid,'(A)') buf
 read(fid,'(A)') buf

 do while(.true.)
  read(fid,'(A)') buf
  call upper(buf(2:5))
  if(buf(2:5) == '$END') exit

  do while(.true.)
   read(fid,'(A)') buf
   read(buf,*,iostat=i) str, nline
   if(i /= 0) exit

   do i = 1, nline, 1
    read(fid,'(A)') buf
   end do ! for i
  end do ! for while

  natom = natom + 1
 end do ! for while

 close(fid)
 return
end subroutine read_natom_from_gms_inp

! read 3 arrays elem, nuc, coor, and the total charge as well as multiplicity
! from a given .gjf file
subroutine read_elem_and_coor_from_gjf(gjfname, natom, elem, nuc, coor, charge, mult)
 use fch_content, only: elem2nuc
 implicit none
 integer :: i, k, fid, nblank
 integer, intent(in) :: natom
 integer, intent(out) :: charge, mult, nuc(natom)
 integer, parameter :: iout = 6
 real(kind=8), intent(out) :: coor(3,natom)
 real(kind=8), parameter :: Bohr_const = 0.52917721092d0
 character(len=2), intent(out) :: elem(natom)
 character(len=240) :: buf
 character(len=240), intent(in) :: gjfname
 logical :: bohr

 charge = 0; mult = 1; coor = 0.0d0; bohr = .false.
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
  read(fid,*,iostat=k) elem(i), coor(1:3,i)
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

! read elements, nuclear charges and Cartesian coordinates from a GAMESS .inp file
subroutine read_elem_nuc_coor_from_gms_inp(inpname, natom, elem, nuc, coor)
 implicit none
 integer :: i, k, fid, nline
 integer, intent(in) :: natom
 integer, intent(out) :: nuc(natom)
 real(kind=8), parameter :: Bohr_const = 0.52917721092d0
 real(kind=8), allocatable :: nuc1(:)
 real(kind=8), intent(out) :: coor(3,natom)
 character(len=1) :: str
 character(len=2), intent(out) :: elem(natom)
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname
 logical :: bohrs

 allocate(nuc1(natom), source=0.0d0)
 nuc = 0; coor = 0.0d0; elem = ' '

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 ! find in the first 6 lines whether the coordinates are in Angstrom or Bohr
 bohrs = .false.
 do i = 1, 6
  read(fid,'(A)') buf
  call upper(buf)
  if(INDEX(buf,'UNITS=BOHR') /= 0) then
   bohrs = .true.
   exit
  end if
 end do ! for i
 ! Angstrom/Bohr determined

 rewind(fid)
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:6) == '$DATA') exit
 end do ! for while

 read(fid,'(A)') buf
 read(fid,'(A)') buf

 do k = 1, natom, 1
  read(fid,*) elem(k), nuc1(k), coor(1:3,k)

  do while(.true.)
   read(fid,'(A)') buf
   read(buf,*,iostat=i) str, nline
   if(i /= 0) exit

   do i = 1, nline, 1
    read(fid,'(A)') buf
   end do ! for do
  end do ! for while

 end do ! for k

 close(fid)

 forall(i = 1:natom) nuc(i) = DNINT(nuc1(i))
 deallocate(nuc1)
 if(bohrs) coor = coor*Bohr_const
 return
end subroutine read_elem_nuc_coor_from_gms_inp

! generate a RHF/UHF .gjf file
subroutine generate_hf_gjf(gjfname, natom, elem, coor, charge, mult, basis,&
                           uhf, cart, DKH2, mem, nproc)
 implicit none
 integer :: i, fid
 integer, intent(in) :: natom, charge, mult, mem, nproc
 integer, parameter :: iout = 6
 real(kind=8), intent(in) :: coor(3,natom)
 character(len=2), intent(in) :: elem(natom)
 character(len=13), intent(in) :: basis
 character(len=240) :: chkname
 character(len=240), intent(in) :: gjfname
 logical, intent(in) :: uhf, cart, DKH2

 i = index(gjfname, '.gjf', back=.true.)
 chkname = gjfname(1:i-1)//'.chk'

 open(newunit=fid,file=TRIM(gjfname),status='replace')
 write(fid,'(A)') '%chk='//TRIM(chkname)
 write(fid,'(A,I0,A)') '%mem=',mem,'GB'
 write(fid,'(A,I0)') '%nprocshared=', nproc
 write(fid,'(A)',advance='no') '#p scf(xqc,maxcycle=128) nosymm '

 if(DKH2) then
  write(fid,'(A)',advance='no') 'int(nobasistransform,DKH2) '
 else
  write(fid,'(A)',advance='no') 'int=nobasistransform '
 end if

 if(mult == 1) then ! singlet
  if(uhf) then
   write(fid,'(A)',advance='no') 'UHF/'//TRIM(basis)//' guess=mix stable=opt '
  else
   write(fid,'(A)',advance='no') 'RHF/'//TRIM(basis)//' '
  end if
 else               ! not singlet
  if(uhf) then
   write(fid,'(A)',advance='no') 'UHF/'//TRIM(basis)//' stable=opt '
  else
   write(iout,'(A)') 'ERROR in subroutine generate_hf_gjf: this molecule is&
                    & not singlet, but UHF is not specified.'
   stop
  end if
 end if

 if(cart) then
  write(fid,'(A)') '6D 10F'
 else
  write(fid,'(A)') '5D 7F'
 end if

 write(fid,'(/,A,/)') 'HF file generated by AutoMR of MOKIT'
 write(fid,'(I0,1X,I0)') charge, mult

 do i = 1, natom, 1
  write(fid,'(A2,3X,3F15.8)') elem(i), coor(1:3,i)
 end do ! for i

 ! Gaussian default    : Gaussian function distribution
 ! Gaussian iop(3/93=1): point nuclei charge distribution
 ! GAMESS default      : point nuclei charge distribution
 ! I found that if iop(3/93=1) is used initially, SCF sometimes converges slowly,
 ! so I use a --Link1-- to add iop(3/93=1) later
 if(DKH2) then
  write(fid,'(/,A)') '--Link1--'
  write(fid,'(A)') '%chk='//TRIM(chkname)
  write(fid,'(A,I0,A)') '%mem=',mem,'GB'
  write(fid,'(A,I0)') '%nprocshared=', nproc
  write(fid,'(A)',advance='no') '#p scf(xqc,maxcycle=128)'
  if(uhf) then
   write(fid,'(A)',advance='no') ' UHF stable=opt'
  else
   write(fid,'(A)',advance='no') ' RHF'
  end if
  write(fid,'(A)') ' chkbasis nosymm guess=read geom=allcheck iop(3/93=1)&
                   & int(nobasistransform,DKH2)'
 end if

 write(fid,'(/)',advance='no')
 close(fid)
 return
end subroutine generate_hf_gjf

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

! detect the number of columns of data in a string buf
function detect_ncol_in_buf(buf) result(ncol)
 implicit none
 integer :: i, ncol
 character(len=24), allocatable :: sbuf(:)
 character(len=240), intent(in) :: buf

 ncol = 0
 do while(.true.)
  read(buf,*,iostat=i) 
  if(i /= 0) exit
  ncol = ncol + 1
 end do ! for while

 return
end function detect_ncol_in_buf

