! read some kind of energy from a specified file

! read GVB electronic energy from a given GAMESS .gms file
subroutine read_gvb_energy_from_gms_gms(gmsname, gvb_e)
 implicit none
 integer :: i, j, fid
 real(kind=8), intent(out) :: gvb_e
 character(len=240) :: buf
 character(len=240), intent(in) :: gmsname

 gvb_e = 0d0
 open(newunit=fid,file=TRIM(gmsname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:10) == 'FINAL GVB') exit
 end do

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_gvb_energy_from_gms_gms: no GVB en&
                   &ergy found in'
  write(6,'(A)') 'file '//TRIM(gmsname)
  write(6,'(/,A)') 'You can open this file and check whether the SCF oscillates.'
  write(6,'(A)') 'If yes, reducing the number of processors and re-run may do&
                 & dome help.'
  write(6,'(A)') "If no, check if there is any error message like 'gamess.01.x&
                 & could not be found'."
  write(6,'(A)') 'In the latter case, you should read Section 4.4.10 in MOKIT&
                 & manual.'
  stop
 end if

 i = INDEX(buf, 'IS')
 j = INDEX(buf, 'AFTER')
 read(buf(i+2:j-1),*) gvb_e

 if(DABS(gvb_e) < 1d-5) then
  write(6,'(/,A)') 'ERROR in subroutine read_gvb_energy_from_gms_gms: it seems &
                   &that GVB'
  write(6,'(A)') 'computation is not converged. gmsname='//TRIM(gmsname)
  stop
 end if
end subroutine read_gvb_energy_from_gms_gms

! read GVB energy from a Gaussian output file
subroutine read_gvb_energy_from_gau_log(logname, gvb_e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: gvb_e
 character(len=240) :: buf
 character(len=240), intent(in) :: logname

 gvb_e = 0d0
 open(newunit=fid,file=TRIM(logname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(8:19) == 'TOTAL ENERGY') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_gvb_energy_from_gau_log: no 'TOTAL&
                   & ENERGY' found"
  write(6,'(A)') 'in file '//TRIM(logname)
  stop
 end if

 read(buf(34:),*) gvb_e
end subroutine read_gvb_energy_from_gau_log

! read GVB energy from a Q-Chem output file
subroutine read_gvb_energy_from_qchem_out(outname, gvb_e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: gvb_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 gvb_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:13) == 'GVB Converge') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_gvb_energy_from_qchem_out: no 'GVB&
                   & Converge' found"
  write(6,'(A)') 'in file '//TRIM(outname)
  stop
 end if

 read(buf(16:),*) i, gvb_e
end subroutine read_gvb_energy_from_qchem_out

! read GVB energy from a GAMESS/Gaussian/Q-Chem output file
subroutine read_gvb_energy_from_output(gvb_prog, outname, gvb_e)
 implicit none
 real(kind=8), intent(out) :: gvb_e
!f2py intent(out) :: gvb_e
 character(len=10), intent(in) :: gvb_prog
!f2py intent(in) :: gvb_prog
 character(len=240), intent(in) :: outname
!f2py intent(in) :: outname

 select case(TRIM(gvb_prog))
 case('gamess')
  call read_gvb_energy_from_gms_gms(outname, gvb_e)
 case('gaussian')
  call read_gvb_energy_from_gau_log(outname, gvb_e)
 case('qchem')
  call read_gvb_energy_from_qchem_out(outname, gvb_e)
 case default
  write(6,'(/,A)') 'ERROR in subroutine read_gvb_energy_from_output: gvb_prog c&
                   &annot be recognized.'
  write(6,'(A)') 'GVB_prog='//TRIM(gvb_prog)
  stop
 end select
end subroutine read_gvb_energy_from_output

! read CASCI/CASSCF energy from a Gaussian .log file
subroutine read_cas_energy_from_gau_log(outname, e, scf)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: e(2)
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: scf

 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(13:22)=='EIGENVALUE' .or. buf(20:29)=='EIGENVALUE' .or. &
     buf(23:32)=='Eigenvalue') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cas_energy_from_gau_log: no 'EIGEN&
                   &VALUE' or 'Eigenvalue'"
  write(6,'(A)') 'found in file '//TRIM(outname)
  stop
 end if

 i = INDEX(buf, 'lue')
 if(i == 0) i = INDEX(buf, 'LUE')

 if(scf) then
  read(buf(i+3:),*) e(2) ! CASSCF
 else
  read(buf(i+3:),*) e(1) ! CASCI
 end if

 if(scf) then ! read CASCI energy
  open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
  do while(.true.)
   read(fid,'(A)') buf
   if(buf(2:8) == 'ITN=  1') exit
  end do ! for while
  close(fid)
  i = INDEX(buf, 'E=')
  read(buf(i+2:),*) e(1)
 end if
end subroutine read_cas_energy_from_gau_log

! read CASCI/CASSCF energy from a PySCF output file
subroutine read_cas_energy_from_pyout(outname, e, scf, spin, dmrg)
 implicit none
 integer :: i, j, k, fid
 integer, intent(in) :: spin ! na - nb
 real(kind=8) :: s_square, expect
 real(kind=8), intent(out) :: e(2)
 real(kind=8), parameter :: max_diff = 1d-3
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 character(len=48), parameter :: err_str = 'ERROR in subroutine read_cas_energy&
                                           &_from_pyout: '
 logical, intent(in) :: scf, dmrg
 logical :: state_specific

 e = 0d0; i = 0; j = 0; k = 0; state_specific = .false.
 s_square = 0d0; expect = 0d0
 expect = 0.5d0*DBLE(spin)
 expect = expect*(expect + 1d0)

 if(scf) then ! (DMRG-)CASSCF
  open(newunit=fid,file=TRIM(outname),status='old',position='append')
  do while(.true.)
   BACKSPACE(fid,iostat=k)
   if(k /= 0) exit
   BACKSPACE(fid,iostat=k)
   if(k /= 0) exit
   read(fid,'(A)',iostat=k) buf
   if(k /= 0) exit
   if(buf(1:23) == '1-step CASSCF converged') then
    i = 1; exit
   end if
   if(buf(1:27) == '1-step CASSCF not converged') then
    j = 1; exit
   end if
   if(buf(1:3) == 'SSS') state_specific = .true.
  end do ! for while

 else ! (DMRG-)CASCI

  if(dmrg) then
   open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
  else
   open(newunit=fid,file=TRIM(outname),status='old',position='append')
   do while(.true.)
    BACKSPACE(fid,iostat=k)
    if(k /= 0) exit
    BACKSPACE(fid,iostat=k)
    if(k /= 0) exit
    read(fid,'(A)',iostat=k) buf
    if(k /= 0) exit
    if(buf(1:15) == 'CASCI converged') then
     i = 1; exit
    end if
    if(buf(1:19) == 'CASCI not converged') then
     j = 1; exit
    end if
   end do ! for while
  end if
 end if

 if(k /= 0) then
  write(6,'(/,A)') TRIM(err_str)//'the file'
  write(6,'(A)') TRIM(outname)//' seems problematic.'
  close(fid)
  stop
 else ! k = 0
  if(j/=0 .and. (.not.state_specific)) then
   write(6,'(/,A)') TRIM(err_str)//'CASCI or CASSCF not converged.'
   close(fid)
   stop
  end if
 end if

 ! read CASCI/CASSCF energy in CASCI/CASSCF job
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:7) == 'CASCI E') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') TRIM(err_str)//"'CASCI E' keyword not found in"
  write(6,'(A)') 'file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, '=')
 if(scf) then
  read(buf(i+1:),*) e(2)
 else
  read(buf(i+1:),*) e(1)
 end if

 i = INDEX(buf, 'S^2 =', back=.true.)
 if(i == 0) then
  write(6,'(/,A)') TRIM(err_str)//"'S^2 =' not found in the desired line."
  write(6,'(A)') "buf='"//TRIM(buf)//"'"
  write(6,'(A)') "Probably your PySCF version is so old that 'S^2 =' is not pri&
                 &nted."
  close(fid)
  stop
 end if
 read(buf(i+5:),*) s_square

 if(DABS(expect - s_square) > max_diff) then
  write(6,'(/,A)') REPEAT('-',79)
  write(6,'(A)') 'Warning from subroutine read_cas_energy_from_pyout: <S**2> de&
                 &viates too much'
  write(6,'(2(A,F11.5))') 'from the expectation value. Expectation=', expect, &
                          ', S_square=', s_square
  write(6,'(A)') 'If this is a ground state calculation, it is probably because&
                 & this CASSCF is'
  write(6,'(A)') 'converged to a wrong spin state. You may try to add the keywo&
                 &rd CrazyWFN in'
  write(6,'(A)') 'mokit{} in .gjf file.'
  write(6,'(A)') 'If this is an excited state calculation where the spin of the&
                 & target excited'
  write(6,'(A)') 'state is different from that of the ground state, you can ign&
                 &ore this warning.'
  write(6,'(A)') REPEAT('-',79)
 end if

 ! Note: in a CASSCF job, there is also a CASCI energy, read it.
 if(scf) then
  rewind(fid)

  do while(.true.)
   read(fid,'(A)') buf
   if(buf(1:9) == 'CASCI E =') exit
  end do ! for while

  close(fid)
  read(buf(10:),*) e(1)
  i = INDEX(buf, '=', back=.true.)
  read(buf(i+1:),*) s_square

  if( DABS(expect - s_square) > max_diff) then
   write(6,'(/,A)') REPEAT('-',79)
   write(6,'(A)') 'Warning in subroutine read_cas_energy_from_pyout: the 0-th s&
                  &tep in this CASSCF'
   write(6,'(A)') 'job, i.e. the CASCI <S**2> deviates too much from the expect&
                  &ation value.'
   write(6,'(2(A,F11.5))') 'Expectation=', expect, ', S_square=', s_square
   write(6,'(A)') 'If this is a ground state calculation, it is probably becaus&
                  &e this CASCI is'
   write(6,'(A)') 'unconverged, or converged to a wrong spin state. If this CAS&
                  &CI energy is'
   write(6,'(A)') 'useless to you, or if the following CASSCF happens to be con&
                  &verged to the'
   write(6,'(A)') 'desired spin, you can ignore this warning. Otherwise, you ma&
                  &y try to add'
   write(6,'(A)') 'the keyword CrazyWFN in mokit{} in .gjf file.'
   write(6,'(A)') 'If this is an excited state calculation where the spin of th&
                  &e target excited'
   write(6,'(A)') 'state is different from that of the ground state, you can ig&
                  &nore this warning.'
   write(6,'(A)') REPEAT('-',79)
  end if
 else
  close(fid)
 end if
end subroutine read_cas_energy_from_pyout

! read CASCI/CASSCF energy from the GAMESS output file
subroutine read_cas_energy_from_gms_gms(outname, e, scf, spin)
 implicit none
 integer :: i, fid
 integer, intent(in) :: spin
 real(kind=8) :: s_square, expect
 real(kind=8), intent(out) :: e(2)
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: scf

 expect = DBLE(spin)/2d0
 expect = expect*(expect + 1d0)
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 if(scf) then  ! CASSCF job
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(2:24) == 'THE DENSITIES ARE STATE') exit
  end do ! for while
 
  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_cas_energy_from_gms_gms: no 'THE &
                    &DENSITIES ARE STATE'"
   write(6,'(A)') 'found in file '//TRIM(outname)
   close(fid)
   stop
  end if

  read(fid,'(A)') buf
  i = INDEX(buf,'ENERGY=')
  read(buf(i+7:),*) e(1)   ! CASCI energy in the CASSCF job
  i = INDEX(buf,'=',back=.true.)
  read(buf(i+1:),*) s_square
  s_square = s_square*(s_square+1d0)
  if( DABS(expect - s_square) > 1D-2) then
   write(6,'(/,A)') 'ERROR in subroutine read_cas_energy_from_gms_gms: in this &
                    &CASSCF job, the 0-th'
   write(6,'(A)') 'step, i.e., the CASCI <S**2> deviates too much from the expe&
                  &ctation value.'
   write(6,'(2(A,F10.6))') 'expectation = ', expect, ', s_square=', s_square
   stop
  end if

  do while(.true.)
   read(fid,'(A)') buf
   if(buf(2:13) == 'STATE   1  E') exit
  end do ! for while
  i = INDEX(buf,'ENERGY=')
  read(buf(i+7:),*) e(2)   ! CASSCF energy
  i = INDEX(buf,'S=')
  read(buf(i+2:),*) s_square
  s_square = s_square*(s_square+1d0)
  if( DABS(expect - s_square) > 1D-2) then
   write(6,'(/,A)') 'ERROR in subroutine read_cas_energy_from_gms_gms: CASSCF <&
                    &S**2> deviates too'
   write(6,'(2(A,F10.6))') 'much from the expectation value. expectation = ', &
                           expect, ', s_square=', s_square
   stop
  end if

 else          ! CASCI job
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(2:20) == 'DENSITY MATRIX WILL') exit
  end do ! for while
 
  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_cas_energy_from_gms_gms: no 'DENS&
                    &ITY MATRIX' found"
   write(6,'(A)') 'in file '//TRIM(outname)
   stop
  end if
 
  i = INDEX(buf,'=', back=.true.)
  read(buf(i+1:),*) s_square
  s_square = s_square*(s_square+1d0)
  if( DABS(expect - s_square) > 1D-2) then
   write(6,'(/,A)') 'ERROR in subroutine read_cas_energy_from_gms_gms: CASCI <S&
                    &**2> deviates too'
   write(6,'(2(A,F10.6))') 'much from the expectation value. expectation = ', &
                           expect, ', s_square=', s_square
   stop
  end if

  read(fid,'(A)') buf
  read(fid,'(A)') buf
  i = INDEX(buf,'=', back=.true.)
  read(buf(i+1:),*) e(1)
 end if

 close(fid)
end subroutine read_cas_energy_from_gms_gms

! read CASCI/CASSCF energy from a given OpenMolcas/Molcas output file
subroutine read_cas_energy_from_molcas_out(outname, e, scf)
 implicit none
 integer :: i, j, fid
 real(kind=8) :: add
 real(kind=8), intent(out) :: e(2)
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 character(len=53), parameter :: error_warn = 'ERROR in subroutine read_cas_en&
                                              &ergy_from_molcas_out: '
 logical, intent(in) :: scf

 e = 0d0; add = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:21) == 'Nr of preliminary CI') exit
 end do ! for while

 read(fid,'(A)') buf
 if(INDEX(buf,'No convergence') > 0) then
  if(scf) then
   write(6,'(/,A)') 'Warning in subroutine read_cas_energy_from_molcas_out:'
   write(6,'(A)') 'The CASCI iterative diagonalization fails to converge. This &
                  &is a defect'
   write(6,'(A)') 'of OpenMolcas when doing CASSCF. If you want a correct CASCI&
                  & energy, please'
   write(6,'(A)') 'run a single CASCI job. This may or may not affect the final&
                  & CASSCF result,'
   write(6,'(A)') 'so the program will continue.'
  else ! CASCI
   write(6,'(/,A)') error_warn
   write(6,'(A)') 'The CASCI iterative diagonalization fails to converge.'
   close(fid)
   stop
  end if
  read(fid,'(A)') buf
 end if

 if(INDEX(buf,'Total energies') > 0) then
  i = INDEX(buf,'Add'); j = INDEX(buf,'au')
  read(buf(i+3:j-1),*) add
  read(fid,'(A)') buf
 end if

 read(buf,*) i, j, i, j, e(1) ! CASCI energy
 e(1) = e(1) + add
 if(.not. scf) then ! CASCI
  close(fid)
  return
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(INDEX(buf,'Convergence after') > 0) exit
 end do ! for while

 if(i /= 0) then
  close(fid)
  write(6,'(/,A)') error_warn
  write(6,'(A)') "'Convergence after' is not found in file "//TRIM(outname)
  stop
 end if

 read(fid,*) i, j, i, j, e(2) ! CASSCF energy
 e(2) = e(2) + add
 close(fid)
end subroutine read_cas_energy_from_molcas_out

! read CASCI/CASSCF energy from a specified ORCA output file
subroutine read_cas_energy_from_orca_out(outname, e, scf)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: e(2)
 character(len=51), parameter :: error_warn = 'ERROR in subroutine read_cas_ene&
                                              &rgy_from_orca_out: '
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: scf

 e = 0d0

 if(scf) then ! CASSCF
  open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:9) == '   E(CAS)') exit
  end do ! for while

  if(i /= 0) then
   write(6,'(/,A)') error_warn//"'   E(CAS)' not found in"
   write(6,'(A)') 'file '//TRIM(outname)
   close(fid)
   stop
  end if
  i = INDEX(buf, '=')
  read(buf(i+1:),*) e(1)

  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:19) == 'Final CASSCF energy') exit
  end do ! for while
 
  if(i /= 0) then
   write(6,'(/,A)') error_warn//"'Final CASSCF energy' not found in"
   write(6,'(A)') 'file '//TRIM(outname)
   close(fid)
   stop
  end if
  i = INDEX(buf, ':')
  read(buf(i+1:),*) e(2)

 else ! CASCI
  open(newunit=fid,file=TRIM(outname),status='old',position='append')

  do while(.true.)
   BACKSPACE(fid)
   BACKSPACE(fid)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:5) == 'STATE') exit
  end do ! for while
 
  if(i /= 0) then
   write(6,'(/,A)') error_warn//"'^STATE' not found in file "//TRIM(outname)
   close(fid)
   stop
  end if

  i = INDEX(buf, '=')
  read(buf(i+1:),*) e(1)
 end if

 close(fid)
end subroutine read_cas_energy_from_orca_out

! read CASCI/CASSCF energy from a given Molpro output file
subroutine read_cas_energy_from_molpro_out(outname, e, scf)
 implicit none
 integer :: i, j, k, fid
 real(kind=8), intent(out) :: e(2)
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: scf

 e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  select case(buf(2:8))
  case('ITER. M','ITER  M','ITE  MI')
   exit
  end select
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cas_energy_from_molpro_out: 'ITER.&
                   & M' not found"
  write(6,'(A)') 'in file '//TRIM(outname)
  write(6,'(A)') 'Error termination of the Molpro CASCI job.'
  close(fid)
  stop
 end if

 do i = 1, 4
  read(fid,fmt=*,iostat=j) k,k,k,k, e(1) ! CASCI energy
  if(j == 0) exit
 end do ! for i
 !read(fid,'(A)') buf
 !read(fid,'(A)') buf
 !if(LEN_TRIM(buf) == 0) then ! Multipassing in transformation
 ! read(fid,'(A)') buf
 ! read(fid,'(A)') buf
 !end if
 !read(buf,*) k,k,k,k, e(1) ! CASCI energy

 if(scf) then
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(2:9) == '!MCSCF S') exit
  end do ! for while
  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_cas_energy_from_molpro_out: '!MCS&
                    &CF S' not found"
   write(6,'(A)') 'in file '//TRIM(outname)
   write(6,'(A)') 'Error termination of the Molpro CASSCF job.'
   close(fid)
   stop
  end if
  i = INDEX(buf, 'Energy')
  read(buf(i+6:),*) e(2) ! CASSCF energy
 end if

 close(fid)
end subroutine read_cas_energy_from_molpro_out

! read CASCI/CASSCF energy from a given BDF output file
! Note: BDF changes output format frequently, one must frequently update this subroutine
subroutine read_cas_energy_from_bdf_out(outname, e, scf)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: e(2)
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: scf

 e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 if(scf) then
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:11) == 'mcscf_eneci') exit
  end do ! for while

  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_cas_energy_from_bdf_out: 'mcscf_e&
                    &neci' not found"
   write(6,'(A)') 'in file '//TRIM(outname)
   write(6,'(A)') 'Error termination of the BDF CASSCF job.'
   close(fid)
   stop
  end if
  read(buf(15:),*) e(1) ! CASCI energy in CASSCF job
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(3:19) == 'CHECKDATA:MCSCF:M') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cas_energy_from_bdf_out: 'CHECKDAT&
                   &A:MCSCF:M' not"
  write(6,'(A)') 'found in file '//TRIM(outname)
  write(6,'(A)') 'Error termination of the BDF CASCI/CASSCF job.'
  close(fid)
  stop
 end if

 i = INDEX(buf,':', back=.true.)
 if(scf) then
  read(buf(i+1:),*) e(2) ! CASSCF energy in CASSCF job
 else
  read(buf(i+1:),*) e(1) ! CASCI energy in CASCI job
 end if

 close(fid)
end subroutine read_cas_energy_from_bdf_out

! read CASCI/CASSCF energy from a given PSI4 output file
subroutine read_cas_energy_from_psi4_out(outname, e, scf)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: e(2)
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: scf

 e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(scf) then
   if(buf(1:12)=='        Iter' .or. buf(1:15)=='           Iter') exit
  else
   if(buf(5:19) == 'Total CI energy') exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cas_energy_from_psi4_out: no 'Iter&
                   &' found in file"
  write(6,'(A)') TRIM(outname)
  close(fid)
  stop
 end if

 if(.not. scf) then ! CASCI
  close(fid)
  i = INDEX(buf,'=')
  read(buf(i+1:),*) e(1)
  return
 end if

 read(fid,'(A)') buf
 ! sometimes there will be extra output in PSI4, e.g.
 ! '(sem_iter): H0block_->H0b_diag'...
 if(INDEX(buf,'MCSCF  1:') == 0) then
  do while(.true.)
   read(fid,'(A)') buf
   if(INDEX(buf,'MCSCF  1:') > 0) exit
  end do ! for while
 end if

 i = INDEX(buf,':')
 read(buf(i+1:),*) e(1)

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(INDEX(buf,'MCSCF Final E') > 0) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cas_energy_from_psi4_out: no 'MCSC&
                   &F Final E' found"
  write(6,'(A)') 'in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf,':')
 read(buf(i+1:),*) e(2)
 close(fid)
end subroutine read_cas_energy_from_psi4_out

! read CASCI/CASSCF energy from a given Dalton output file
subroutine read_cas_energy_from_dalton_out(outname, e, scf)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: e(2)
 character(len=1) :: str
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: scf

 e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:19) == '@ Final CI energies') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cas_energy_from_dalton_out: no '@ &
                   &Final CI energies' found"
  write(6,'(A)') 'in file '//TRIM(outname)
  close(fid)
  stop
 end if

 read(fid,*) str, i, e(1) ! CASCI energy
 if(.not. scf) then
  close(fid)
  return
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:19) == '@    Final MCSCF en') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cas_energy_from_dalton_out: no &
                   &'@    Final MCSCF en' found"
  write(6,'(A)') 'in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf,':')
 read(buf(i+1:),*) e(2)
 close(fid)
end subroutine read_cas_energy_from_dalton_out

! read CASCI/CASSCF energy from a Gaussian/PySCF/GAMESS/OpenMolcas/ORCA output file
subroutine read_cas_energy_from_output(cas_prog, outname, spin, scf, dmrg, &
                                       ptchg_e, nuc_pt_e, e)
 implicit none
 integer, intent(in) :: spin
!f2py intent(in) :: spin
 real(kind=8), intent(in) :: ptchg_e, nuc_pt_e
!f2py intent(in) :: ptchg_e, nuc_pt_e
 real(kind=8), intent(out) :: e(2)
!f2py intent(out) :: e
 character(len=10), intent(in) :: cas_prog
!f2py intent(in) :: cas_prog
 character(len=240), intent(in) :: outname
!f2py intent(in) :: outname
 logical, intent(in) :: scf, dmrg
!f2py intent(in) :: scf, dmrg

 select case(TRIM(cas_prog))
 case('gaussian')
  call read_cas_energy_from_gau_log(outname, e, scf)
 case('gamess')
  call read_cas_energy_from_gms_gms(outname, e, scf, spin)
  e(1) = e(1) + ptchg_e + nuc_pt_e
 case('pyscf')
  call read_cas_energy_from_pyout(outname, e, scf, spin, dmrg)
  e = e + ptchg_e
 case('openmolcas')
  call read_cas_energy_from_molcas_out(outname, e, scf)
  e = e + ptchg_e
 case('orca')
  call read_cas_energy_from_orca_out(outname, e, scf)
 case('molpro')
  call read_cas_energy_from_molpro_out(outname, e, scf)
  e = e + ptchg_e
 case('bdf')
  call read_cas_energy_from_bdf_out(outname, e, scf)
  e = e + ptchg_e + nuc_pt_e
 case('psi4')
  call read_cas_energy_from_psi4_out(outname, e, scf)
  e = e + ptchg_e + nuc_pt_e
 case('dalton')
  call read_cas_energy_from_dalton_out(outname, e, scf)
  e = e + ptchg_e
 case default
  write(6,'(/,A)') 'ERROR in subroutine read_cas_energy_from_output: CAS_prog c&
                   &annot be identified.'
  write(6,'(A)') 'CAS_prog='//TRIM(cas_prog)
  stop
 end select
end subroutine read_cas_energy_from_output

! read CASSCF OVB-MP2 energy from a Gaussian output file
subroutine read_mrpt_energy_from_gau_log(outname, ref_e, corr_e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: ref_e, corr_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 ref_e = 0d0; corr_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(13:22) == 'EIGENVALUE') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_mrpt_energy_from_gau_log: no 'EIGE&
                   &NVALUE' found in"
  write(6,'(A)') 'file '//TRIM(outname)
  close(fid)
 end if
 read(buf(23:),*) ref_e

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(28:34) == 'EUMP2 =') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_mrpt_energy_from_gau_log: no 'EUMP&
                   &2 =' found in"
  write(6,'(A)') 'file '//TRIM(outname)
 end if

 read(buf(35:),*) corr_e
 corr_e = corr_e - ref_e
end subroutine read_mrpt_energy_from_gau_log

! read NEVPT2 energy from PySCF output file
subroutine read_mrpt_energy_from_pyscf_out(outname, troot, ref_e, corr_e)
 implicit none
 integer :: i, k, fid
 integer, intent(in) :: troot ! 0 for the ground state, >0 for excited state
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8), intent(out) :: ref_e, corr_e

 ref_e = 0d0; corr_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:8) == 'Nevpt2 E') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_mrpt_energy_from_pyscf_out: no 'Ne&
                   &vpt2 E' found in"
  write(6,'(A)') 'file '//TRIM(outname)
  close(fid)
  stop
 end if
 call get_dpv_after_flag(buf, '=', .true., corr_e)

 if(troot == 0) then
  do while(.true.)
   BACKSPACE(fid,iostat=i)
   if(i /= 0) exit
   BACKSPACE(fid,iostat=i)
   if(i /= 0) exit
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:7) == 'CASCI E') exit
  end do ! for while
 else
  do while(.true.)
   BACKSPACE(fid,iostat=i)
   if(i /= 0) exit
   BACKSPACE(fid,iostat=i)
   if(i /= 0) exit
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:11) == 'CASCI state') then
    read(buf(12:),*) k
    if(k == troot) exit
   end if
  end do ! for while
 end if

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_mrpt_energy_from_pyscf_out: no CAS&
                   &CI energy found'
  write(6,'(A)') 'in file '//TRIM(outname)
  write(6,'(A,I0)') 'troot=', troot
  stop
 end if

 call get_dpv_after_flag(buf, '=', .true., ref_e)
end subroutine read_mrpt_energy_from_pyscf_out

! read CASTP2 energy from OpenMolcas output file
subroutine read_mrpt_energy_from_molcas_out(outname, itype, ref_e, corr_e)
 implicit none
 integer :: i, j, fid
 integer, intent(in) :: itype ! 1/2/3 for SC-NEVPT2/FIC-NEVPT2/CASPT2
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8) :: rtmp
 real(kind=8), intent(out) :: ref_e, corr_e

 ref_e = 0d0; corr_e = 0d0
 if(itype<1 .or. itype>3) then
  write(6,'(/,A,I0)') 'ERROR in subroutine read_mrpt_energy_from_molcas_out: in&
                      &valid itype=', itype
  write(6,'(A)') 'Allowed values are 1/2/3 for SC-NEVPT2/FIC-NEVPT2/CASPT2.'
  stop
 end if

 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
 select case(itype)
 case(1,2) ! NEVPT2
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(2:30) == 'Energies of zeroth-order DMRG') exit
  end do ! for while

  if(i /= 0) then
   write(6,'(/,A)') 'ERROR in subroutine read_mrpt_energy_from_molcas_out: No r&
                    &eference energy'
   write(6,'(A)') 'found in file '//TRIM(outname)
   close(fid)
   stop
  end if

  read(fid,'(A)') buf
  read(fid,'(A)') buf
  i = INDEX(buf, '=')
  read(buf(i+1:),*) ref_e

  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(2:18) == 'state number:   1') exit
  end do ! for while
 
  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_mrpt_energy_from_molcas_out: No &
                    &'state number:   1'"
   write(6,'(A)') 'found in file '//TRIM(outname)
   close(fid)
   stop
  end if

  do j = 1, 15
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(2:6) == 'Total') exit
  end do ! for j

  if(i/=0 .or. j==16) then
   write(6,'(/,A)') 'ERROR in subroutine read_mrpt_energy_from_molcas_out: No N&
                    &EVPT2 energy found'
   write(6,'(A)') 'in file '//TRIM(outname)
   close(fid)
   stop
  end if

  if(itype == 1) then ! SC-NEVPT2
   read(buf(7:),*) rtmp, corr_e
  else                ! FIC-NEVPT2
   read(buf(7:),*) rtmp, rtmp, rtmp, corr_e
  end if

 case(3) ! CASPT2
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(7:23) == 'Reference energy:') exit
  end do ! for while
 
  if(i /= 0) then
   write(6,'(/,A)') 'ERROR in subroutine read_mrpt_energy_from_molcas_out: No r&
                    &eference energy'
   write(6,'(A)') 'found in file '//TRIM(outname)
   close(fid)
   stop
  end if

  i = INDEX(buf,':')
  read(buf(i+1:),*) ref_e

  do i = 1, 3
   read(fid,'(A)') buf
  end do ! for i 

  i = INDEX(buf, ':')
  read(buf(i+1:),*) corr_e
  corr_e = corr_e - ref_e

 case default
  write(6,'(/,A,I0)') 'ERROR in subroutine read_mrpt_energy_from_molcas_out: in&
                      &valid itype=', itype
  stop
 end select

 close(fid)
end subroutine read_mrpt_energy_from_molcas_out

! read NEVPT2/CASPT2 energy from a Molpro output file
subroutine read_mrpt_energy_from_molpro_out(outname, itype, ref_e, corr_e)
 implicit none
 integer :: i, fid
 integer, intent(in) :: itype ! 1/2/3 for SC-NEVPT2/FIC-NEVPT2/CASPT2
 real(kind=8), intent(out) :: ref_e, corr_e
 character(len=8), parameter :: key(3) = ['Strongly','!NEVPT2 ','!RSPT2 S']
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 ref_e = 0d0; corr_e = 0d0
 if(itype<1 .or. itype>3) then
  write(6,'(/,A)') 'ERROR in subroutine read_mrpt_energy_from_molpro_out: itype&
                   & out of range.'
  write(6,'(A,I0,A)') 'itype=', itype, ', outname='//TRIM(outname)
  stop
 end if

 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:20)=='!MCSCF STATE  1.1 E' .or. buf(2:19)=='!MCSCF STATE 1.1 E') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_mrpt_energy_from_molpro_out: CASSC&
                   &F energy cannot'
  write(6,'(A)') "be found in file "//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf,'ergy')
 read(buf(i+4:),*) ref_e

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:9) == key(itype)) exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_mrpt_energy_from_molpro_out:'
  write(6,'(A)') "'"//key(itype)//"' not found in file "//TRIM(outname)
  stop
 end if

 i = INDEX(buf,'ergy')
 read(buf(i+4:),*) corr_e
 corr_e = corr_e - ref_e
end subroutine read_mrpt_energy_from_molpro_out

! read NEVPT2/CASPT2 energy from a ORCA .out file
subroutine read_mrpt_energy_from_orca_out(outname, ref_e, corr_e)
 implicit none
 integer :: i, fid
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8), intent(out) :: ref_e, corr_e

 ref_e = 0d0; corr_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(4:20) == 'Total Energy Corr') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_mrpt_energy_from_orca_out: no 'Tot&
                   &al Energy Corr'"
  write(6,'(A)') 'found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 call get_dpv_after_flag(buf, '=', .true., corr_e)
 read(fid,'(A)') buf
 read(fid,'(A)') buf
 close(fid)
 call get_dpv_after_flag(buf, '=', .true., ref_e)
end subroutine read_mrpt_energy_from_orca_out

! read MRMP2 energy from GAMESS output file (.gms)
subroutine read_mrpt_energy_from_gms_out(outname, ref_e, corr_e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: ref_e, corr_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 ref_e = 0d0; corr_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:16) == 'TOTAL   (MCSCF)') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_mrpt_energy_from_gms_out: no refer&
                   &ence energy found'
  write(6,'(A)') 'in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, '=')
 read(buf(i+1:),*) ref_e

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(22:42) == '2ND ORDER ENERGY CORR') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_mrpt_energy_from_gms_out: no MRMP2&
                   & energy found in'
  write(6,'(A)') 'file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, '=')
 read(buf(i+1:),*) corr_e
 close(fid)
end subroutine read_mrpt_energy_from_gms_out

! read NEVPT2 energy from a BDF .out file
subroutine read_mrpt_energy_from_bdf_out(outname, itype, ref_e, corr_e, dav_e)
 implicit none
 integer :: i, fid
 integer, intent(in) :: itype   ! 1/2 for SDSPT2/FIC-NEVPT2
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8) :: rtmp(6)
 real(kind=8), intent(out) :: ref_e, corr_e, dav_e
 ! dav_e: Davidson correction energy for size-inconsistency error

 ref_e = 0d0; corr_e = 0d0; dav_e = 0d0
 if(.not. (itype==1 .or. itype==2)) then
  write(6,'(/,A)') 'ERROR in subroutine read_mrpt_energy_from_bdf_out: currentl&
                   &y only reading SDSPT2/'
  write(6,'(A)') 'NEVPT2 energy is supported.'
  stop
 end if

 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:11) == 'Print final') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_mrpt_energy_from_bdf_out: 'Print f&
                   &inal' not found in"
  write(6,'(A)') 'file '//TRIM(outname)
  write(6,'(A)') 'Error termination of the BDF CASSCF in SDSPT2/NEVPT2 job.'
  close(fid)
  stop
 end if

 do i = 1, 4
  read(fid,'(A)') buf
 end do ! for i
 i = INDEX(buf, '=')
 read(buf(i+1:),*) ref_e ! CASCI/CASSCF energy

 if(itype == 1) then
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(2:11) == 'MRPT2 calc') exit
  end do ! for while
  read(fid,'(A)') buf
  read(fid,'(A)') buf
  read(fid,*) i, rtmp(1:6), dav_e
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:13) == 'TARGET_XIANCI') exit
 end do ! for while
 close(fid)

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_mrpt_energy_from_bdf_out: no 'TARG&
                   &ET_XIANCI' found in"
  write(6,'(A)') 'file '//TRIM(outname)
  stop
 end if

 read(buf(26:),*) corr_e
 if(itype == 1) dav_e = dav_e - corr_e
 corr_e = corr_e - ref_e
end subroutine read_mrpt_energy_from_bdf_out

! read CASPT3 electronic energy from Molpro .out file
subroutine read_caspt3_energy_from_molpro_out(outname, ref_e, corr2_e, corr3_e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: ref_e, corr2_e, corr3_e
 character(len=19), parameter :: key(6) = ['!MCSCF STATE  1.1 E', &
 '!MCSCF STATE 1.1 En', '!RSPT2 STATE  1.1 E', '!RSPT2 STATE 1.1 En', &
 '!RSPT3 STATE  1.1 E', '!RSPT3 STATE 1.1 En']
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 ref_e = 0d0; corr2_e = 0d0; corr3_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:20)==key(1) .or. buf(2:20)==key(2)) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_caspt3_energy_from_molpro_out:'
  write(6,'(A)') 'CASSCF energy not found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf,'ergy')
 read(buf(i+4:),*) ref_e

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:20)==key(3) .or. buf(2:20)==key(4)) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_caspt3_energy_from_molpro_out:'
  write(6,'(A)') 'CASPT2 energy not found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf,'ergy')
 read(buf(i+4:),*) corr2_e
 corr2_e = corr2_e - ref_e

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:20)==key(5) .or. buf(2:20)==key(6)) exit
 end do ! for while
 close(fid)

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_caspt3_energy_from_molpro_out:'
  write(6,'(A)') 'CASPT3 energy not found in file '//TRIM(outname)
  stop
 end if

 i = INDEX(buf,'ergy')
 read(buf(i+4:),*) corr3_e
 corr3_e = corr3_e - ref_e
end subroutine read_caspt3_energy_from_molpro_out

! read NEVPT3 electronic energy from BDF .out file
subroutine read_nevpt3_energy_from_bdf_out(outname, ref_e, corr2_e, corr3_e)
 implicit none
 integer :: i, fid
 real(kind=8) :: rtmp
 real(kind=8), intent(out) :: ref_e, corr2_e, corr3_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 ref_e = 0d0; corr2_e = 0d0; corr3_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:6) == 'NROOT') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_mrpt_nevpt3_from_bdf_out: no 'NROO&
                   &T' found in"
  write(6,'(A)') 'file '//TRIM(outname)
  close(fid)
  stop
 end if

 read(fid,*) i, ref_e, corr2_e, rtmp, corr3_e
 close(fid)
 corr2_e = corr2_e - ref_e
 corr3_e = corr3_e - ref_e
end subroutine read_nevpt3_energy_from_bdf_out

! read FIC-NEVPT3/FIC-NEVPT4(SD) electronic energy from ORCA output file
subroutine read_nevpt34_e_from_orca_out(outname, pt4, ref_e, corr_e)
 implicit none
 integer :: i, fid
 real(kind=8) :: rtmp
 real(kind=8), intent(out) :: ref_e, corr_e(3)
 real(kind=8), parameter :: thres = 1d-5
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: pt4

 ref_e = 0d0; corr_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:5) == 'E(0)') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_nevpt34_e_from_orca_out: no 'E(0)'&
                   & found in"
  write(6,'(A)') 'file '//TRIM(outname)
  close(fid)
  stop
 end if
 call get_dpv_after_flag(buf, '=', .true., ref_e) ! E(0)

 read(fid,'(A)') buf
 call get_dpv_after_flag(buf, '=', .true., rtmp) ! E(1)
 if(DABS(rtmp) > thres) then
  write(6,'(/,A)') 'ERROR in subroutine read_nevpt34_e_from_orca_out: E(1) is n&
                   &ot zero.'
  write(6,'(A)') 'Unexpected case.'
  close(fid)
  stop
 end if

 read(fid,'(A)') buf
 call get_dpv_after_flag(buf, '=', .true., corr_e(1)) ! E(2)

 read(fid,'(A)') buf
 close(fid)
 call get_dpv_after_flag(buf, '=', .true., corr_e(2)) ! E(3)
 corr_e(2) = corr_e(2) + corr_e(1)

 if(pt4) then
  call get_dpv_after_flag(buf, '=', .true., corr_e(3)) ! E(3)
  corr_e(3) = corr_e(3) + corr_e(2)
 end if
end subroutine read_nevpt34_e_from_orca_out

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

! read Davidson correction and MRCISD energy from OpenMolcas, ORCA, Gaussian or
! Molpro output file
subroutine read_mrci_energy_from_output(CtrType, mrcisd_prog, outname, ptchg_e,&
                                        nuc_pt_e, davidson_e, e)
 implicit none
 integer :: i, fid
 integer, intent(in) :: CtrType
 real(kind=8) :: casci_e, ref_weight
 real(kind=8), intent(in) :: ptchg_e, nuc_pt_e
 real(kind=8), intent(out) :: davidson_e, e
 character(len=10), intent(in) :: mrcisd_prog
 character(len=240), intent(in) :: outname
 character(len=240) :: buf
 logical, external :: has_davidson_q

 davidson_e = 0d0; e = 0d0; ref_weight = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 select case(TRIM(mrcisd_prog))
 case('pyscf')
  if(CtrType == 3) then ! DMRG-FIC-MRCISD
   do while(.true.)
    BACKSPACE(fid)
    BACKSPACE(fid)
    read(fid,'(A)') buf
    if(buf(1:7) == 'E(MRCI)') exit
   end do ! for while
   i = INDEX(buf, 'DC =', back=.true.)
   read(buf(i+4:),*) davidson_e
   read(fid,'(A)') buf
   i = INDEX(buf, '=')
   read(buf(i+1:),*) e
  else
   write(6,'(/,A)') 'ERROR in subroutine read_mrci_energy_from_output: invalid &
                    &CtrType.'
   write(6,'(A,I0)') 'Current CtrType=', CtrType
   close(fid)
   stop
  end if
  e = e + ptchg_e
 case('openmolcas')
  if(CtrType == 1) then ! uncontracted MRCISD
   do while(.true.)
    BACKSPACE(fid)
    BACKSPACE(fid)
    read(fid,'(A)') buf
    if(buf(1:29) == '::    RASSCF root number  1 T') exit
   end do ! for while
   i = INDEX(buf,':',back=.true.)
   read(buf(i+1:),*) e
  else if(CtrType == 2) then ! ic-MRCISD
   do while(.true.)
    BACKSPACE(fid)
    BACKSPACE(fid)
    read(fid,'(A)') buf
    if(buf(23:31) == 'CI ENERGY') exit
   end do ! for while
   i = INDEX(buf,':',back=.true.)
   read(buf(i+1:),*) e
   read(fid,'(A)') buf
   i = INDEX(buf,':',back=.true.)
   read(buf(i+1:),*)  davidson_e
  end if
  e = e + ptchg_e

 case('orca')
  if(CtrType == 1) then ! uncontracted MRCISD
   do while(.true.)
    BACKSPACE(fid)
    BACKSPACE(fid)
    read(fid,'(A)') buf
    if(buf(20:27) == 'E(MR-CI)') exit
   end do ! for while
   read(fid,'(A)') buf
   read(fid,'(A)') buf
   read(buf(13:),*) e
   do while(.true.)
    BACKSPACE(fid)
    BACKSPACE(fid)
    read(fid,'(A)') buf
    if(buf(14:26) == 'residual conv') exit
   end do ! for while
   read(fid,'(A)') buf
   read(buf(67:),*) ref_weight
   do while(.true.)
    BACKSPACE(fid)
    BACKSPACE(fid)
    read(fid,'(A)') buf
    if(buf(1:17) == 'Computing the ref') exit
   end do ! for while
   BACKSPACE(fid)
   BACKSPACE(fid)
   BACKSPACE(fid)
   read(fid,'(A)') buf
   read(buf(22:),*) casci_e
   davidson_e = (1d0 - ref_weight)*(e - casci_e)
  else if(CtrType == 3) then ! FIC-MRCISD
   do while(.true.)
    BACKSPACE(fid)
    BACKSPACE(fid)
    read(fid,'(A)') buf
    if(buf(1:5) == 'ETOT ') exit
   end do ! for while
   i = INDEX(buf,'...',back=.true.)
   read(buf(i+3:),*) e
   read(fid,'(A)') buf
   read(fid,'(A)') buf
   i = INDEX(buf,'...',back=.true.)
   read(buf(i+3:),*) davidson_e
  end if

 case('gaussian') ! uncontracted MRCISD
  do while(.true.)
   BACKSPACE(fid)
   BACKSPACE(fid)
   read(fid,'(A)') buf
   !               Davidson                           Lanczos
   if(buf(17:32)=='Final Eigenvalue' .or. buf(20:29)=='EIGENVALUE') exit
  end do ! for while
  i = INDEX(buf,'lue',back=.true.)
  if(i == 0) i = INDEX(buf,'LUE',back=.true.)
  read(buf(i+3:),*) e

 case('molpro')
  do while(.true.)
   BACKSPACE(fid)
   BACKSPACE(fid)
   read(fid,'(A)') buf
   if(INDEX(buf(2:10),'!Total e') > 0) exit
  end do ! for while

  read(buf(17:),*) e
  read(fid,'(A)') buf
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) read(fid,'(A)') buf
  i = INDEX(buf,':')
  read(buf(i+1:),*) davidson_e
  davidson_e = davidson_e - e
  e = e + ptchg_e

 case('psi4')
  do while(.true.)
   BACKSPACE(fid)
   BACKSPACE(fid)
   read(fid,'(A)') buf
   if(buf(5:19) == 'Total CI energy') exit

   if(buf(11:23) == 'Psi4: An Open') then
    write(6,'(/,A)') 'ERROR in subroutine read_mrci_energy_from_output:'
    write(6,'(A)') "No 'Total CI energy' found in file "//TRIM(outname)
    close(fid)
    stop
   end if
  end do ! for while

  i = INDEX(buf,'=')
  read(buf(i+1:),*) e
  e = e + ptchg_e + nuc_pt_e

 case('dalton')
  do while(.true.)
   BACKSPACE(fid)
   BACKSPACE(fid)
   read(fid,'(A)') buf
   if(buf(1:19) == '@ Final CI energies') exit

   if(buf(22:32) == 'Dalton - An') then
    write(6,'(/,A)') 'ERROR in subroutine read_mrci_energy_from_output:'
    write(6,'(A)') "No '@ Final CI energies' found in file "//TRIM(outname)
    close(fid)
    stop
   end if
  end do ! for while

  read(fid,'(A)') buf
  read(buf(7:),*) e
  e = e + ptchg_e

 case('gamess')
  close(fid)

  if(has_davidson_q(outname)) then
   open(newunit=fid,file=TRIM(outname),status='old',position='append')
   do while(.true.)
    BACKSPACE(fid)
    BACKSPACE(fid)
    read(fid,'(A)') buf
    if(buf(1:13) == ' E(MR-CISD) =') exit
   end do ! for while

   read(buf(14:),*) e
   read(fid,'(A)') buf
   read(fid,'(A)') buf
   read(buf(14:),*) ref_weight

   do while(.true.)
    read(fid,'(A)') buf
    if(buf(2:7) == 'E(REF)') exit
   end do ! for while
   i = INDEX(buf, '=')
   read(buf(i+1:),*) casci_e

   ! GAMESS uses renormalized Davidson correction. Here we compute the Davidson
   ! correction
   davidson_e = (1d0 - ref_weight)*(e - casci_e)

  else ! no 'CALC. OF DAVIDSON', i.e. no Davidson Q

   open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
   do while(.true.)
    read(fid,'(A)')  buf
    if(buf(2:14) == 'CI EIGENSTATE') exit
   end do ! for while

   close(fid)
   i = INDEX(buf, '=')
   read(buf(i+1:),*) e
  end if

 case default
  write(6,'(/,A)') 'ERROR in subroutine read_mrci_energy_from_output: invalid m&
                   &rcisd_prog='//TRIM(mrcisd_prog)
  close(fid)
  stop
 end select

 close(fid)
end subroutine read_mrci_energy_from_output

subroutine read_mrcc_energy_from_output(mrcc_prog, mrcc_type, outname, ref_e, &
                                        corr_e)
 implicit none
 integer :: i, fid
 integer, intent(in) :: mrcc_type ! 1/2/3 for FIC-/Mk-/BW-MRCCSD
 ! 4/5 for Mk-/BW-MRCCSD(T), 6~8 for GVB-BCCC2b/3b/4b
 real(kind=8) :: rtmp
 real(kind=8), intent(out) :: ref_e, corr_e(2)
 character(len=240) :: buf
 character(len=10), intent(in) :: mrcc_prog
 character(len=240), intent(in) :: outname

 ref_e = 0d0; corr_e = 0d0
 ! corr_e(1): FIC-NEVPT2 energy
 ! corr_e(2): FIC-/Mk-/BW- MRCCSD energy

 select case(TRIM(mrcc_prog))
 case('orca')
  open(newunit=fid,file=TRIM(outname),status='old',position='append')
  do while(.true.)
   BACKSPACE(fid,iostat=i)
   if(i /= 0) exit
   BACKSPACE(fid,iostat=i)
   if(i /= 0) exit
   read(fid,'(A)') buf
   if(buf(1:17) == ' Summary of multi') exit
  end do ! for while

  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_mrcc_energy_from_output: ' Summary&
                   & of multi' not found"
   write(6,'(A)') 'in file '//TRIM(outname)
   close(fid)
   stop
  end if

  read(fid,'(A)') buf
  read(fid,'(A)') buf
  read(fid,*) i,i,rtmp,rtmp,ref_e,corr_e(2)

  do while(.true.)
   BACKSPACE(fid,iostat=i)
   if(i /= 0) exit
   BACKSPACE(fid,iostat=i)
   if(i /= 0) exit
   read(fid,'(A)') buf
   if(buf(1:5) == 'MULT=') exit
  end do ! for while
  close(fid)

  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_mrcc_energy_from_output: 'MULT='&
                   & not found in file "//TRIM(outname)
   stop
  end if

  i = INDEX(buf,'EC=')
  read(buf(i+3:),*) corr_e(1)

 case('gvb_bcci2b','gvb_bccc2b','gvb_bccc3b')
  open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:5) == 'GVB e') exit
  end do ! for while

  close(fid)
  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_mrcc_energy_from_output: no &
                   &'GVB e' found in"
   write(6,'(A)') 'file '//TRIM(outname)
   stop
  end if

  i = INDEX(buf, '=')
  read(buf(i+1:),*) ref_e

  ! find if any 'ERROR: amplitude iterations fail'
  open(newunit=fid,file=TRIM(outname),status='old',position='append')
  do while(.true.)
   BACKSPACE(fid)
   BACKSPACE(fid)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:32) == 'ERROR: amplitude iterations fail') then
    close(fid)
    write(6,'(A)') 'ERROR in subroutine read_mrcc_energy_from_output: amplitud&
                   &e iterations not converge.'
    stop
   end if
   if(buf(1:3) == 'ITN') exit
  end do ! for while

  close(fid)
  open(newunit=fid,file=TRIM(outname),status='old',position='append')

  do while(.true.)
   BACKSPACE(fid)
   BACKSPACE(fid)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(index(buf,'=') > 0) exit
  end do ! for while

  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  i = INDEX(buf, '=')
  read(buf(i+1:),*) corr_e(1) ! GVB-BCCCnb correlation energy

  read(fid,'(A)') buf
  i = INDEX(buf, '=')
  read(buf(i+1:),*) corr_e(2) ! GVB-BCCCnb total energy
  close(fid)

 !case('nwchem')
 case default
  write(6,'(/,A)') 'ERROR in subroutine read_mrcc_energy_from_output: invalid m&
                   &rcc_prog='//TRIM(mrcc_prog)
  stop
 end select
end subroutine read_mrcc_energy_from_output

