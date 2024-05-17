! moved from rwgeom.f90

! read Cartesian gradient from a given file
subroutine read_grad_from_output(prog_name, outname, natom, grad)
 implicit none
 integer, intent(in) :: natom
 real(kind=8) :: e
 real(kind=8), intent(out) :: grad(3*natom)
 character(len=10), intent(in) :: prog_name
 character(len=240), intent(in) :: outname

 select case(TRIM(prog_name))
 case('bdf')
  call read_grad_from_bdf_out(outname, natom, grad)
 case('cfour')
  call read_grad_from_cfour_out(outname, natom, grad)
 case('dalton')
  call read_grad_from_dalton_out(outname, natom, grad)
 case('gamess')
  call read_grad_from_dat(outname, natom, grad)
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
 integer :: i, k, fid
 integer, intent(in) :: natom
 real(kind=8), intent(out) :: grad(3*natom)
 character(len=2) :: elem = ' '
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 grad = 0d0

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

