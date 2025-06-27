! moved from rwgeom.f90

! read Cartesian gradient from a given file
subroutine read_grad_from_output(prog_name, outname, natom, grad)
 implicit none
 integer :: i
 integer, intent(in) :: natom
!f2py intent(in) :: natom
 real(kind=8) :: e
 real(kind=8), intent(out) :: grad(3*natom)
!f2py intent(out) :: grad
!f2py depend(natom) :: grad
 character(len=10), intent(in) :: prog_name
!f2py intent(in) :: prog_name
 character(len=240), intent(in) :: outname
!f2py intent(in) :: outname

 select case(TRIM(prog_name))
 case('bdf')
  call read_grad_from_bdf_out(outname, natom, grad)
 case('cfour')
  call read_grad_from_cfour_out(outname, natom, grad)
 case('dalton')
  call read_grad_from_dalton_out(outname, natom, grad)
 case('gamess')
  i = LEN_TRIM(outname)
  select case(outname(i-3:i))
  case('.dat')
   call read_grad_from_dat(outname, natom, grad)
  case('.gms')
   call read_grad_from_gms_gms(outname, natom, grad)
  case default
   write(6,'(/,A)') 'ERROR in subroutine read_grad_from_output: filetype not re&
                    &cognized.'
   write(6,'(A)') 'Suffix='//outname(i-3:i)
   stop
  end select
 case('gaussian')
  call read_grad_from_gau_log(outname, natom, grad)
 case('molpro')
  call read_grad_from_molpro_out(outname, natom, grad)
 case('orca')
  call read_grad_from_engrad(outname, natom, e, grad)
 case('psi4')
  call read_grad_from_psi4_out(outname, natom, grad)
 case('pyscf')
  call read_grad_from_pyscf_out(outname, natom, grad)
 case('qchem')
  call read_grad_from_qchem131(natom, grad, .true.)
 case('openmolcas')
  call read_grad_from_molcas_out(outname, natom, grad)
 case default
  write(6,'(/,A)') 'ERROR in subroutine read_grad_from_output: program cannot b&
                   &e identified.'
  write(6,'(A)') 'prog_name='//TRIM(prog_name)
  stop
 end select

 write(6,'(/,A)') 'Cartesian gradients (HARTREE/BOHR):'
 write(6,'(5(1X,ES15.8))') grad
end subroutine read_grad_from_output

! read Cartesian gradient from a given PySCF output file
subroutine read_grad_from_pyscf_out(outname, natom, grad)
 implicit none
 integer :: i, k, fid
 integer, intent(in) :: natom
 real(kind=8), intent(out) :: grad(3*natom)
 character(len=3) :: elem
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 grad = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(INDEX(buf,'gradients') /= 0) exit
 end do ! for while

 read(fid,'(A)') buf

 do i = 1, natom, 1
  read(fid,*) k, elem, grad(3*i-2:3*i)
 end do ! for i

 close(fid)
end subroutine read_grad_from_pyscf_out

! read Cartesian gradient from a given .fch file
subroutine read_grad_from_fch(fchname, natom, grad)
 implicit none
 integer :: i, fid
 integer, intent(in) :: natom
 real(kind=8), intent(out) :: grad(3*natom)
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname

 grad = 0d0
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:26) == 'Opt point       1 Gradient') exit
 end do ! for while

 read(fid,'(5(1X,ES15.8))') (grad(i),i=1,3*natom)
 close(fid)
end subroutine read_grad_from_fch

! read Cartesian gradient from a given Gaussian .log file
subroutine read_grad_from_gau_log(logname, natom, grad)
 implicit none
 integer :: i, k1, k2, fid
 integer, intent(in) :: natom
 real(kind=8), intent(out) :: grad(3*natom)
 character(len=240) :: buf
 character(len=240), intent(in) :: logname

 grad = 0d0
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
end subroutine read_grad_from_gau_log

! read Cartesian gradient from a given GAMESS .gms file
subroutine read_grad_from_gms_gms(outname, natom, grad)
 implicit none
 integer :: i, k, icase, fid
 integer, intent(in) :: natom
 real(kind=8) :: r
 real(kind=8), intent(out) :: grad(3*natom)
 character(len=2) :: elem = ' '
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 grad = 0d0; icase = 0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(buf(26:47) == 'GRADIENT OF THE ENERGY') then
   icase = 1
   exit
  end if
  if(buf(34:47) == 'GRADIENT (HART') then
   icase = 2
   exit
  end if
 end do ! for while

 do i = 1, 3
  read(fid,'(A)') buf
 end do ! for i

 select case(icase)
 case(1)
  do i = 1, natom, 1
   read(fid,*) k, elem, grad(3*i-2:3*i)
  end do ! for i
 case(2)
  do i = 1, natom, 1
   read(fid,*) k, elem, r, grad(3*i-2:3*i)
  end do ! for i
 case default
  write(6,'(/,A)') 'ERROR in subroutine read_grad_from_gms_gms: icase out of ra&
                   &nge!'
  write(6,'(A)') 'outname='//TRIM(outname)
  close(fid)
  stop
 end select

 close(fid)
end subroutine read_grad_from_gms_gms

! read Cartesian gradient from a given GAMESS .dat file
subroutine read_grad_from_dat(datname, natom, grad)
 implicit none
 integer :: i, fid
 integer, intent(in) :: natom
 real(kind=4) :: r
 real(kind=8), intent(out) :: grad(3*natom)
 character(len=2) :: elem
 character(len=240) :: buf
 character(len=240), intent(in) :: datname

 grad = 0d0
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
end subroutine read_grad_from_dat

! read Cartesian gradient from a given (Open)Molcas .out file
subroutine read_grad_from_molcas_out(outname, natom, grad)
 implicit none
 integer :: i, fid
 integer, intent(in) :: natom
 real(kind=8), intent(out) :: grad(3*natom)
 character(len=8) :: str = ' '
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 grad = 0d0
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
  read(fid,*) str, grad(3*i-2:3*i)
 end do ! for i

 close(fid)
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

 grad = 0d0
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
end subroutine read_grad_from_orca_out

! read Cartesian gradients from a given ORCA .engrad file
subroutine read_grad_from_engrad(engrad, natom, e, grad)
 implicit none
 integer :: i, fid
 integer, intent(in) :: natom
 real(kind=8), intent(out) :: e, grad(3*natom)
 character(len=240) :: buf
 character(len=240), intent(in) :: engrad

 e = 0d0; grad = 0d0
 open(newunit=fid,file=TRIM(engrad),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:17) == '# The current tot') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_grad_from_engrad: no '# The curren&
                   &t tot' found in"
  write(6,'(A)') 'file '//TRIM(engrad)
  close(fid)
  stop
 end if

 read(fid,'(A)') buf
 read(fid,*) e

 do i = 1, 3
  read(fid,'(A)') buf
 end do
 read(fid,*) grad

 close(fid)
end subroutine read_grad_from_engrad

! read CASSCF or CASPT2 Cartesian gradients from a given Molpro output file
subroutine read_grad_from_molpro_out(outname, natom, grad)
 implicit none
 integer :: i, k, fid
 integer, intent(in) :: natom
 real(kind=8), intent(out) :: grad(3*natom)
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 grad = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(INDEX(buf,'GRADIENT FOR S') > 0) exit
 end do ! for while

 do i = 1, 3
  read(fid,'(A)') buf
 end do ! for i

 do i = 1, natom, 1
  read(fid,*) k, grad(3*i-2:3*i)
 end do ! for i

 close(fid)
end subroutine read_grad_from_molpro_out

! read Cartesian gradient from a given BDF .out file
subroutine read_grad_from_bdf_out(outname, natom, grad)
 implicit none
 integer :: i, fid
 integer, intent(in) :: natom
 real(kind=8), intent(out) :: grad(3*natom)
 character(len=10) :: str
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 grad = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(INDEX(buf,'Molecular gradient - Mol') /= 0) exit
 end do ! for while

 read(fid,'(A)') buf
 read(fid,'(A)') buf

 do i = 1, natom, 1
  read(fid,*) str, grad(3*i-2:3*i)
 end do ! for i

 close(fid)
end subroutine read_grad_from_bdf_out

subroutine read_grad_from_psi4_out(outname, natom, grad)
 implicit none
 integer :: i, fid
 integer, intent(in) :: natom
 character(len=10) :: str
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8), intent(out) :: grad(3*natom)

 grad = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(4:13)=='Total Grad' .or. buf(4:13)=='Total grad') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_grad_from_psi4_out: no 'Total Grad&
                   &' found in"
  write(6,'(A)') 'file '//TRIM(outname)
  close(fid)
  stop
 end if

 read(fid,'(A)') buf
 read(fid,'(A)') buf
 do i = 1, natom, 1
  read(fid,*) str, grad(3*i-2:3*i)
 end do ! for i
 close(fid)
end subroutine read_grad_from_psi4_out

subroutine read_grad_from_qchem131(natom, grad, deleted)
 implicit none
 integer :: i, fid
 integer, intent(in) :: natom
 real(kind=8), intent(out) :: grad(3*natom)
 logical, intent(in) :: deleted

 grad = 0d0
 open(newunit=fid,file='131.0',status='old',access='stream')
 read(fid,iostat=i) grad

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_grad_from_qchem_out: failed to rea&
                   &d nuclear gradients'
  write(6,'(A)') 'from file 131.0'
  close(fid)
  stop
 else ! read successfully
  if(deleted) then
   close(fid,status='delete')
  else
   close(fid)
  end if
 end if
end subroutine read_grad_from_qchem131

subroutine read_grad_from_cfour_out(outname, natom, grad)
 implicit none
 integer :: i, fid
 integer, intent(in) :: natom
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8), intent(out) :: grad(3*natom)

 grad = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(3:10) == 'gradient') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_grad_from_cfour_out: keyword 'grad&
                   &ient' not found"
  write(6,'(A)') 'in file '//TRIM(outname)
  close(fid)
  stop
 end if

 read(fid,*) grad
 close(fid)
end subroutine read_grad_from_cfour_out

subroutine read_grad_from_dalton_out(outname, natom, grad)
 implicit none
 integer :: i, fid
 integer, intent(in) :: natom
 character(len=3) :: str
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8), intent(out) :: grad(3*natom)

 grad = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(30:43)=='Molecular grad' .or. buf(33:46)=='Molecular grad') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_grad_from_dalton_out: keyword 'Mol&
                   &ecular grad' not"
  write(6,'(A)') 'found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 read(fid,'(A)') buf
 read(fid,'(A)') buf

 do i = 1, natom, 1
  read(fid,*) str, grad(3*i-2:3*i)
 end do ! for i

 close(fid)
end subroutine read_grad_from_dalton_out

! read electronic energy, atomic forces and stress tensor from a CP2K output file
! Note: forces are negative gradients.
subroutine read_efs_from_cp2k_out(outname, natom, e, force, stress)
 implicit none
 integer :: i, j, k, fid
 integer, intent(in) :: natom
 real(kind=8), intent(out) :: e, force(3,natom), stress(3,3) ! a.u.
 !                      Hartree, Hartree/Bohr  , GPa
 character(len=1) :: str1
 character(len=2) :: elem
 character(len=7) :: str7
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical :: old_cp2k

 e = 0d0; force = 0d0; stress = 0d0
 call require_file_exist(outname)
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:13) == 'ENERGY| Tota') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_force_from_cp2k_out: failed to loc&
                   &ate 'ENERGY| Tota'"
  write(6,'(A)') 'in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, ":")
 read(buf(i+1:),*) e

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:13) == 'ATOMIC FORCE') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_force_from_cp2k_out: failed to loc&
                   &ate 'ATOMIC FORCE'"
  write(6,'(A)') 'in file '//TRIM(outname)
  close(fid)
  stop
 end if

 read(fid,'(A)') buf
 read(fid,'(A)') buf
 do i = 1, natom, 1
  read(fid,*) j, k, elem, force(:,i)
 end do ! for i

 old_cp2k = .false.
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:13) == 'STRESS| Anal') exit
  if(buf(2:13) == 'STRESS TENSO') then
   read(fid,'(A)') buf
   old_cp2k = .true.
   exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_force_from_cp2k_out: failed to loc&
                   &ate stress tensor'
  write(6,'(A)') 'in file '//TRIM(outname)
  close(fid)
  stop
 end if

 read(fid,'(A)') buf
 if(old_cp2k) then
  do i = 1, 3
   read(fid,*) str1, stress(:,i)
  end do ! for i
 else
  do i = 1, 3
   read(fid,*) str7, str1, stress(:,i)
  end do ! for i
 end if

 close(fid)
end subroutine read_efs_from_cp2k_out

! read Cartesian Force Constants from a .fch(k) file
subroutine read_cart_force_const_from_fch(fchname, natom, cfc)
 implicit none
 integer :: i, j, k, fid
 integer, intent(in) :: natom
!f2py intent(in) :: natom
 real(kind=8), intent(out) :: cfc(3*natom,3*natom)
!f2py intent(out) :: cfc
!f2py depend(natom) :: cfc
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 cfc = 0d0
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:21) == 'Cartesian Force Const') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cart_force_const_from_fch: failed &
                   &to locate 'Cartesian"
  write(6,'(A)') " Force Const' in file "//TRIM(fchname)
  close(fid)
  stop
 end if

 k = 3*natom
 read(fid,*) ((cfc(j,i),j=1,i),i=1,k)
 close(fid)
 call symmetrize_dmat(k, cfc)
end subroutine read_cart_force_const_from_fch

! Write/create an ORCA .hess file which only contains Cartesian Force Constants
! and Cartesian coordinates. Such a .hess file can be used as the initial Hessian
! in an ORCA geometry optimization job.
! A full .hess file includes also spin multiplicity, normal modes, etc, which
! will not be printed in this subroutine.
! Note: the input coor must be in unit Angstrom.
subroutine write_orca_hess(natom, elem, ram, coor, cfc, hessname)
 use phys_cons, only: Bohr_const
 implicit none
 integer :: i, j, k, m, n, nb, fid
 integer, intent(in) :: natom
!f2py intent(in) :: natom
 real(kind=8), allocatable :: coor_b(:,:)
 real(kind=8), intent(in) :: ram(natom), coor(3,natom), cfc(3*natom,3*natom)
!f2py intent(in) :: ram, coor, cfc
!f2py depend(natom) :: ram, coor, cfc
 character(len=2), intent(in) :: elem(natom)
!f2py intent(in) :: elem
!f2py depend(natom) :: elem
 character(len=240), intent(in) :: hessname
!f2py intent(in) :: hessname

 k = 3*natom
 open(newunit=fid,file=TRIM(hessname),status='replace')
 write(fid,'(/,A)') "$orca_hessian_file"
 write(fid,'(/,A)') "$hessian"
 write(fid,'(I0)') k
 nb = k/5

 do i = 1, nb, 1
  write(fid,'(2X,5(9X,I10))') (5*i-6+j, j=1,5)
  m = 5*i - 4
  do j = 1, k, 1
   write(fid,'(I5,3X,5(1X,ES18.10))') j-1, cfc(j,m:m+4)
  end do ! for j
 end do ! for i

 n = k - nb*5
 if(n > 0) then
  write(fid,'(2X,5(9X,I10))') (5*nb-1+j, j=1,n)
  m = 5*nb + 1
  do j = 1, k, 1
   write(fid,'(I5,3X,5(1X,ES18.10))') j-1, cfc(j,m:m+n-1)
  end do ! for j
 end if

 write(fid,'(/,A)') '#'
 write(fid,'(A)') '# The atoms: label  mass x y z (in bohrs)'
 write(fid,'(A)') '#'
 write(fid,'(A)') '$atoms'
 write(fid,'(I0)') natom

 allocate(coor_b(3,natom))
 coor_b = coor/Bohr_const
 do i = 1, natom, 1
  write(fid,'(A2,3X,F10.5,3(1X,F20.10))') elem(i), ram(i), coor_b(:,i)
 end do ! for i

 close(fid)
 deallocate(coor_b)
end subroutine write_orca_hess

! Find whether there is a $HESS section in a specified GAMESS .inp/.dat file, if
! yes, delete the $HESS section.
subroutine del_hess_in_gms_dat(datname)
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=240) :: buf, datname1
 character(len=240), intent(in) :: datname
!f2py intent(in) :: datname
 logical :: find_hess

 find_hess = .false.
 i = INDEX(datname, '.', back=.true.)
 datname1 = datname(1:i-1)//'.t'

 open(newunit=fid,file=TRIM(datname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(buf(2:6) == '$HESS') then
   find_hess = .true.
   exit
  end if
  if(i /= 0) exit
 end do ! for while

 if(find_hess) then
  rewind(fid)
  open(newunit=fid1,file=TRIM(datname1),status='replace')
  do while(.true.)
   read(fid,'(A)') buf
   if(buf(2:6) == '$HESS') exit
   write(fid1,'(A)') TRIM(buf)
  end do ! for while
  do while(.true.)
   read(fid,'(A)') buf
   if(buf(2:5) == '$END') exit
  end do ! for while
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   write(fid1,'(A)') TRIM(buf)
  end do ! for while
  close(fid,status='delete')
  close(fid1)
  i = RENAME(TRIM(datname1), TRIM(datname))
 else
  close(fid)
 end if
end subroutine del_hess_in_gms_dat

! Write/print Cartesian Force Constants into a specified GAMESS .inp/.dat file
subroutine write_hess2gms_inp(n, cfc, inpname)
 implicit none
 integer :: i, j, k, nline, fid
 integer, intent(in) :: n ! = 3*natom
!f2py intent(in) :: n
 real(kind=8), intent(in) :: cfc(n,n)
!f2py intent(in) :: cfc
!f2py depend(n) :: cfc
 character(len=240), intent(in) :: inpname
!f2py intent(in) :: inpname

 call del_hess_in_gms_dat(inpname)
 nline = n/5
 if(n > 5*nline) nline = nline + 1
 open(newunit=fid,file=TRIM(inpname),status='old',position='append')
 write(fid,'(A)') ' $HESS'
 write(fid,'(2(A,F19.10))') 'ENERGY IS ',0d0,' E(NUC) IS ',0d0

 do i = 1, n, 1
  k = MOD(i, 100)
  do j = 1, nline-1, 1
   write(fid,'(I2,I3,5ES15.8)') k, MOD(j,1000), cfc(5*j-4:5*j,i)
  end do ! for j
  write(fid,'(I2,I3,5ES15.8)') k, MOD(j,1000), cfc(5*j-4:n,i)
 end do ! for i

 write(fid,'(A)') ' $END'
 close(fid)
end subroutine write_hess2gms_inp

