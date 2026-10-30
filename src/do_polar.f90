! written by jxzou at 20260619: calculate the the molecular polarizability

! Currently only the CASSCF polarizability is supported. We will add more
! interfaces gradually. CASCI polarizability is unsupported.

! In the future, we might rename this to subroutine do_prop. But currently
! we only have the property `polarizability`, so we use do_polar.

!TODO: read CASSCF energy for Polar_prog=Dalton/ORCA, and compare with the
!      CASSCF energy from previous do_cas
!TODO: signal errors when CASCI is required in do_cas
!TODO: support MP2 and CCSD polar

subroutine do_polar()
 use mr_keyword, only: mem, nproc, iroot, polar, molcas_omp, dalton_mpi, bgchg,&
  chgname, casnofch, casscf_polar, casscf_prog, polar_prog, gau_path, orca_path,&
  molcas_path
 use mol, only: nacto, nacte, casscf_e, polarizability
 use util_wrapper, only: add_bgcharge2inp_wrap, unfchk, fch2mkl_wrap, &
  mkl2gbw, fch2inporb_wrap
 implicit none
 integer :: i
 real(kind=8) :: alpha_iso, alpha_aniso(3), e(2)
 real(kind=8), parameter :: e_diff_thres = 1d-4
 character(len=24) :: data_string
 character(len=30), parameter :: error_warn='ERROR in subroutine do_polar: '
 character(len=240) :: buf, proname, inpname, outname, mklname, cas_out
 logical :: copy_output

 if(.not. polar) return
 write(6,'(//,A)') 'Enter subroutine do_polar...'

 if(casscf_polar .and. TRIM(casscf_prog)=='orca' .and. TRIM(polar_prog)=='orca') then
  copy_output = .true.
 else
  copy_output = .false.
 end if

 i = LEN_TRIM(casnofch)
 if(casnofch(i-6:i) == '_NO.fch') then
  proname = casnofch(1:i-6)
 else
  write(6,'(/,A)') error_warn//'unrecognized casnofch.'
  write(6,'(A)') 'casnofch='//TRIM(casnofch)
  stop
 end if

 write(6,'(A,2(I0,A))') 'CASSCF(',nacte,'e,',nacto,'o) polarizability using pro&
                        &gram '//TRIM(polar_prog)

 select case(TRIM(polar_prog))
 case('gaussian')
  call check_exe_exist(gau_path)
  inpname = TRIM(proname)//'pol.gjf'
#ifdef _WIN32
  outname = TRIM(proname)//'pol.out'
#else
  outname = TRIM(proname)//'pol.log'
#endif
  mklname = TRIM(proname)//'pol.chk'
  call prt_cas_gjf(inpname, nacto, nacte, .true., .false., .true.)
  if(bgchg) call add_bgcharge2inp_wrap(chgname, inpname)
  call unfchk(casnofch, mklname)
  call submit_gau_job(gau_path, inpname, .true.)
  call delete_file(TRIM(mklname))

 case('orca')
  outname = TRIM(proname)//'pol.out'
  if(copy_output) then
   write(6,'(A)') 'CASSCF_prog=Polar_prog=ORCA detected. Directly copy and use &
                  &CASSCF output file.'
   cas_out = proname(1:i-7)//'.out'
   call sys_copy_file(TRIM(cas_out), TRIM(outname), .false.)
  else
   call check_exe_exist(orca_path)
   inpname = TRIM(proname)//'pol.inp'
   mklname = TRIM(proname)//'pol.mkl'
   call fch2mkl_wrap(casnofch, mklname, REPEAT(' ',30), .true.)
   call prt_orca_cas_inp(inpname, .true., .true.)
   if(bgchg) call add_bgcharge2inp_wrap(chgname, inpname)
   ! if bgchg = .True., .inp and .mkl file will be updated
   call mkl2gbw(mklname)
   call delete_file(TRIM(mklname))
   call submit_orca_job(orca_path, inpname, .true., .false., .false.)
  end if

 case('openmolcas')
  call check_exe_exist(molcas_path)
  inpname = TRIM(proname)//'pol.input'
  outname = TRIM(proname)//'pol.out'
  call fch2inporb_wrap(casnofch, .false., inpname)
  call prt_cas_molcas_inp(inpname, .true., .true.)
  if(bgchg) call add_bgcharge2inp_wrap(chgname, inpname)
  if(polar) then
   buf = 'echo -e "&LOPROP\n" >> '//TRIM(inpname)
   call run_command(TRIM(buf), .false., .false.)
  end if
  call submit_molcas_job(inpname, mem, nproc, molcas_omp)

 case('dalton')
  i = LEN_TRIM(proname)
  buf = proname(1:i)//'pol'
  inpname = proname(1:i)//'pol.dal'
  outname = proname(1:i)//'pol.out'
  call prt_cas_dalton_prop_inp(casnofch, .true., .false., .true., iroot, i)
  if(bgchg) call add_bgcharge2inp_wrap(chgname, inpname)
  call submit_dalton_job(buf, mem, nproc, dalton_mpi, .false., .false., .false.)
  call read_cas_energy_from_dalton_out(outname, e, .true.)
  if(DABS(e(2)-casscf_e) > e_diff_thres) then
   write(6,'(/,A)') error_warn//'it seems that Dalton CASSCF is not converged'
   write(6,'(A)') 'to the same CASSCF solution from previous step.'
   write(6,'(A,F18.6)') 'Energy difference threshold :', e_diff_thres
   write(6,'(A,F18.6)') 'E(CASSCF) from previous step:', casscf_e
   write(6,'(A,F18.6)') 'E(CASSCF) from Dalton polar :', e(2)
   stop
  end if

 case default
  write(6,'(/,A)') error_warn//'unrecognized Polar_prog='//TRIM(polar_prog)
  stop
 end select

 call read_polarizability_from_output(outname, polar_prog, polarizability)
 write(6,'(/,A)') 'The molecular polarizability (atomic units):'
 do i = 1, 3
  write(6,'(3(1X,ES15.8))') polarizability(:,i)
 end do ! for i

 call calc_alpha_iso_aniso_from_polar(polarizability, alpha_iso, alpha_aniso)
 write(6,'(/,A,F12.3)') 'Isotropic polarizability alpha_iso: ', alpha_iso
 write(6,'(A)') 'Anisotropic polarizability:'
 write(6,'(A)') 'Below `a` means alpha and `e` means eigenvalues of 3*3 polar t&
                &ensor.'
 write(6,'(A)') 'Definition 1): sqrt([(axx-ayy)^2+(axx-azz)^2+(ayy-azz)^2+6(axy&
                &^2+axz^2+ayz^2)]/2)'
 write(6,'(A)') '               Chem. Phys. 2013, 410, 90, DOI: 10.1016/j.chemp&
                &hys.2012.11.005'
 write(6,'(A)') 'Definition 2): sqrt([(axx-ayy)^2+(axx-azz)^2+(ayy-azz)^2)]/2)'
 write(6,'(A)') '               J. Chem. Phys. 1993, 98, 3022, DOI: 10.1063/1.4&
                &64129'
 write(6,'(A)') 'Definition 3): e3 - (e1+e2)/2'
 write(6,'(A,F10.3)') 'alpha_aniso using 1): ', alpha_aniso(1)
 write(6,'(A,F10.3)') 'alpha_aniso using 2): ', alpha_aniso(2)
 write(6,'(A,F10.3)') 'alpha_aniso using 3): ', alpha_aniso(3)

 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine do_polar at '//TRIM(data_string)
end subroutine do_polar

subroutine read_polarizability_from_output(outname, polar_prog, polarizability)
 implicit none
 real(kind=8), intent(out) :: polarizability(3,3)
 character(len=10), intent(in) :: polar_prog
 character(len=240), intent(in) :: outname

 select case(TRIM(polar_prog))
 case('orca')
  call read_polarizability_from_orca_out(outname, polarizability)
 case('dalton')
  call read_polarizability_from_dalton_out(outname, polarizability)
 case('openmolcas')
  call read_polarizability_from_molcas_out(outname, polarizability)
 case('gaussian')
  call read_polarizability_from_gau_log(outname, polarizability)
 case default
  write(6,'(/,A)') 'ERROR in subroutine read_polarizability_from_output: unreco&
                   &gnized Polar_prog='//TRIM(polar_prog)
  write(6,'(A)') 'outname='//TRIM(outname)
  stop
 end select
end subroutine read_polarizability_from_output

subroutine read_polarizability_from_orca_out(outname, polarizability)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: polarizability(3,3)
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 polarizability = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:11) == 'The raw car') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_polarizability_from_orca_out: "The&
                   & raw car"'
  write(6,'(A)') 'not located in file '//TRIM(outname)
  close(fid)
  stop
 end if

 read(fid,*) polarizability(:,1)
 read(fid,*) polarizability(:,2)
 read(fid,*) polarizability(:,3)
 close(fid)
end subroutine read_polarizability_from_orca_out

subroutine read_polarizability_from_dalton_out(outname, polarizability)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: polarizability(3,3)
 character(len=2) :: str2
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 polarizability = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(28:54) == 'Static polarizabilities (au') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_polarizability_from_dalton_out: "S&
                   &tatic polarizabilities (au"'
  write(6,'(A)') 'not located in file '//TRIM(outname)
  close(fid)
  stop
 end if

 do i = 1, 4
  read(fid,'(A)') buf
 end do ! for i

 read(fid,*) str2, polarizability(:,1)
 read(fid,*) str2, polarizability(:,2)
 read(fid,*) str2, polarizability(:,3)
 close(fid)
end subroutine read_polarizability_from_dalton_out

subroutine read_polarizability_from_molcas_out(outname, polarizability)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: polarizability(3,3)
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 polarizability = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(3:17) == 'Molecular Polar') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_polarizability_from_molcas_out: "M&
                   &olecular Polar"'
  write(6,'(A)') 'not located in file '//TRIM(outname)
  close(fid)
  stop
 end if

 read(fid,'(A)') buf
 read(fid,'(A)') buf
 read(fid,*) polarizability(1,1)
 read(fid,*) polarizability(1:2,2)
 read(fid,*) polarizability(1:3,3)
 close(fid)
 polarizability(2,1) = polarizability(1,2)
 polarizability(3,1) = polarizability(1,3)
 polarizability(3,2) = polarizability(2,3)
end subroutine read_polarizability_from_molcas_out

subroutine read_polarizability_from_gau_log(logname, polarizability)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: polarizability(3,3)
 character(len=240) :: buf
 character(len=240), intent(in) :: logname

 polarizability = 0d0
 open(newunit=fid,file=TRIM(logname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:16) == 'Isotropic polar') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_polarizability_from_gau_log: "Isot&
                   &ropic polar" not'
  write(6,'(A)') 'located in file '//TRIM(logname)
  close(fid)
  stop
 end if

 read(fid,'(A)') buf
 read(fid,*) i, polarizability(1,1)
 read(fid,*) i, polarizability(1:2,2)
 read(fid,*) i, polarizability(1:3,3)
 close(fid)
 polarizability(2,1) = polarizability(1,2)
 polarizability(3,1) = polarizability(1,3)
 polarizability(3,2) = polarizability(2,3)
end subroutine read_polarizability_from_gau_log

