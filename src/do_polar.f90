! written by jxzou at 20260619: calculate the the molecular polarizability

! Currently only the CASSCF polarizability is supported. We will add more
! interfaces gradually. CASCI polarizability is unsupported.

! In the future, we might rename this to subroutine do_prop. But currently
! we only have the property `polarizability`, so we use do_polar.

!TODO: automatically set RIC auxbasis
!TODO: support cc-pVnZ-PP and their corresponding RIC auxbasis
!TODO: signal errors when CASCI is required in do_cas
!TODO: support CCSD polar

subroutine do_polar()
 use mr_keyword, only: mem, nproc, polar, molcas_omp, bgchg, chgname, gau_path,&
  orca_path, molcas_path, casnofch, polar_prog
 use mol, only: nacto, nacte, polarizability
 use util_wrapper, only: add_bgcharge2inp_wrap, unfchk, fch2mkl_wrap, &
  mkl2gbw, fch2inporb_wrap
 implicit none
 integer :: i
 character(len=24) :: data_string
 character(len=30), parameter :: error_warn='ERROR in subroutine do_polar: '
 character(len=240) :: buf, proname, inpname, outname, mklname

 if(.not. polar) return
 write(6,'(//,A)') 'Enter subroutine do_polar...'

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
  call check_exe_exist(orca_path)
  inpname = TRIM(proname)//'pol.inp'
  mklname = TRIM(proname)//'pol.mkl'
  outname = TRIM(proname)//'pol.out'
  call fch2mkl_wrap(casnofch, mklname, REPEAT(' ',30), .true.)
  call prt_cas_orca_inp(inpname, .true., .true.)
  if(bgchg) call add_bgcharge2inp_wrap(chgname, inpname)
  ! if bgchg = .True., .inp and .mkl file will be updated
  call mkl2gbw(mklname)
  call delete_file(TRIM(mklname))
  call submit_orca_job(orca_path, inpname, .true., .false., .false.)
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
 case default
  write(6,'(/,A)') error_warn//'unrecognized Polar_prog='//TRIM(polar_prog)
  stop
 end select

 call read_polarizability_from_output(outname, polar_prog, polarizability)
 write(6,'(/,A)') 'The molecular polarizability (atomic units):'
 do i = 1, 3
  write(6,'(3(1X,ES15.8))') polarizability(:,i)
 end do ! for i

 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine do_polar at '//TRIM(data_string)
end subroutine do_polar

subroutine read_polarizability_from_output(outname, polar_prog, polarizability)
 implicit none
 real(kind=8), intent(out) :: polarizability(3,3)
 character(len=10), intent(in) :: polar_prog
 character(len=240), intent(in) :: outname

 select case(TRIM(polar_prog))
 case('gaussian')
  call read_polarizability_from_gau_log(outname, polarizability)
 case('orca')
  call read_polarizability_from_orca_out(outname, polarizability)
 case('openmolcas')
  call read_polarizability_from_molcas_out(outname, polarizability)
 case default
  write(6,'(/,A)') 'ERROR in subroutine read_polarizability_from_output: unreco&
                   &gnized Polar_prog='//TRIM(polar_prog)
  write(6,'(A)') 'outname='//TRIM(outname)
  stop
 end select
end subroutine read_polarizability_from_output

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
  write(6,'(/,A)') 'ERROR in subroutine read_polarizability_from_molcas_out: "T&
                   &he raw car" not'
  write(6,'(A)') 'located in file '//TRIM(outname)
  close(fid)
  stop
 end if

 read(fid,*) polarizability(:,1)
 read(fid,*) polarizability(:,2)
 read(fid,*) polarizability(:,3)
 close(fid)
end subroutine read_polarizability_from_orca_out

