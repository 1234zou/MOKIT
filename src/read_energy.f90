! read some kind of energy from a specified file

subroutine prt_pyscf_forge_recommend()
 implicit none
 write(6,'(/,A)') 'If the pyscf-forge package is not installed, it is recommend&
                  &ed to install it.'
 write(6,'(A)') 'If the pyscf-forge package is already installed and used by MO&
                &KIT, it is'
 write(6,'(A)') 'recommended to specify `HardWFN` or `CrazyWFN` in mokit{} and &
                &re-try.'
end subroutine prt_pyscf_forge_recommend

! read HF electronic energy from a Gaussian .log/.out file
subroutine read_hf_e_and_ss_from_gau_log(logname, first, e, ss)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: e, ss ! HF energy and spin square
 character(len=240) :: buf
 character(len=240), intent(in) :: logname
 logical, intent(in) :: first ! read from the first/last `SCF Done`

 e = 0d0; ss = 0d0

 if(first) then
  open(newunit=fid,file=TRIM(logname),status='old',position='rewind')
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(2:10) == 'SCF Done:') exit
  end do ! for while
 else
  open(newunit=fid,file=TRIM(logname),status='old',position='append')
  do while(.true.)
   BACKSPACE(fid)
   BACKSPACE(fid)
   read(fid,'(A)') buf
   if(buf(2:18) == 'Entering Gaussian') then
    i = -1; exit
   end if
   if(buf(2:10) == 'SCF Done:') then
    i = 0; exit
   end if
  end do ! for while
 end if

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_hf_e_and_ss_from_gau_log: "SCF Don&
                   &e" not found'
  write(6,'(A)') 'in file '//TRIM(logname)
  close(fid)
  stop
 end if
 call get_dpv_after_flag(buf, '=', .true., e)

 ! We do not read <S**2> below 'SCF Done' because when the spin is very high
 !  the format here would become <S**2>=******* due to Fortran features.
 ! Instead, we search the 'S**2 before annihilation' below 'SCF Done'
 do i = 1, 7
  read(fid,'(A)') buf
  if(buf(2:14) == 'S**2 before a') then
   read(buf(26:),*) ss
   exit
  end if
 end do ! for i

 close(fid)
end subroutine read_hf_e_and_ss_from_gau_log

! read spin multiplicity from a PSI4 output file
subroutine read_mult_from_psi4_out(outname, mult)
 implicit none
 integer :: i, fid
 integer, intent(out) :: mult
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(3:14) == 'Multiplicity') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_mult_from_psi4_out: `Multiplicity`&
                   & not found in file'
  write(6,'(A)') TRIM(outname)
  stop
 end if

 call get_int_after_flag(buf, '=', .false., mult)
end subroutine read_mult_from_psi4_out

! read HF/DFT electronic energy from a PSI4 .out file
subroutine read_scf_e_and_ss_from_psi4_out(outname, hf_type, e, ssquare)
 implicit none
 integer :: mult, fid
 integer, intent(in) :: hf_type
 real(kind=8), intent(out) :: e, ssquare
 character(len=53), parameter :: error_warn = 'ERROR in subroutine read_scf_e_a&
                                              &nd_ss_from_psi4_out: '
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 e = 0d0; ssquare = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 ! There are two HF energies in this file, we should begin from the end of file
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(buf(11:23) == 'Psi4: An Open') then
   close(fid)
   write(6,'(/,A)') error_warn//'SCF energy not found in file'
   write(6,'(A)') TRIM(outname)
   stop
  end if
  if(buf(5:18) == 'Total Energy =') exit
 end do ! for while

 call get_dpv_after_flag(buf, ':', .true., e)

 select case(hf_type)
 case(1) ! RHF
  close(fid)
  ssquare = 0d0
 case(2) ! ROHF
  close(fid)
  call read_mult_from_psi4_out(outname, mult)
  ssquare = 0.25d0*DBLE(mult*mult-1)
 case(3) ! UHF
  do while(.true.)
   BACKSPACE(fid)
   BACKSPACE(fid)
   read(fid,'(A)') buf
   if(buf(11:23) == 'Psi4: An Open') then
    close(fid)
    write(6,'(/,A)') error_warn//'`S^2 Observed` not found in file'
    write(6,'(A)') TRIM(outname)
    stop
   end if
   if(buf(5:16) == 'S^2 Observed') exit
  end do ! for while
  close(fid)
  call get_dpv_after_flag(buf, ':', .false., ssquare)
 case default
  write(6,'(A,I0)') error_warn//'invalid hf_type = ', hf_type
  stop
 end select
end subroutine read_scf_e_and_ss_from_psi4_out

! read SCF electronic energy from a CP2K output file
subroutine read_scf_e_from_cp2k_out(outname, scf_e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: scf_e ! in Hartree
!f2py intent(out) :: scf_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
!f2py intent(in) :: outname

 scf_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:13) == 'ENERGY| Tota') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_scf_e_from_cp2k_out: failed to loc&
                   &ate 'ENERGY| Tota'"
  write(6,'(A)') 'in file '//TRIM(outname)
  stop
 end if

 call get_dpv_after_flag(buf, ']', .false., scf_e)
end subroutine read_scf_e_from_cp2k_out

! read IP-/EA-ADC energies from a specified PySCF output file
subroutine read_adc_e_from_pyscf_out(outname, nstate, e)
 implicit none
 integer :: i, fid
 integer, intent(in) :: nstate
 real(kind=8) :: scf_e
 real(kind=8), intent(out) :: e(0:nstate)
 character(len=6) :: str6
 character(len=17), parameter :: adc_str = 'ADC calculation s'
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:13) == 'converged SCF') exit
 end do ! for while
 call get_dpv_after_flag(buf, '=', .true., scf_e) ! HF energy

 do while(.true.)
  read(fid,'(A)') buf
  str6 = buf(1:6)
  if(str6=='E_corr' .or. str6=='MP2 co' .or. str6=='MP3 co') exit
 end do ! for while
 call get_dpv_after_flag(buf, '=', .true., e(0)) ! MP2/MP3 correlation energy
 e(0) = e(0) + scf_e

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(9:25)==adc_str .or. buf(13:29)==adc_str .or. buf(20:36)==adc_str) exit
 end do ! for while
 read(fid,'(A)') buf ! skip 1 line

 do i = 1, nstate, 1
  read(fid,'(A)') buf
  call get_dpv_after_flag(buf, '=', .true., e(i))
  e(i) = e(i) + e(0)
 end do ! for i

 close(fid)
end subroutine read_adc_e_from_pyscf_out

! Read oscillator strengths from a specified ORCA output file. This subroutine
! can be used for ORCA 5/6 ADC(2), EOM-CCSD and SA-CASSCF calculations.
subroutine read_fosc_from_orca_out(outname, nstate, fosc)
 implicit none
 integer :: i, k, fid
 integer, intent(in) :: nstate
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8), intent(out) :: fosc(nstate)
 logical :: orca6

 orca6 = .false.; fosc = 0d0; k = 0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(36:48) == 'O   R   C   A') then
   i = -1
   exit
  end if
  if(buf(10:28) == 'ABSORPTION SPECTRUM') then
   k = 26; exit
  end if
  if(buf(33:51) == 'ABSORPTION SPECTRUM') then
   k = 39; exit
  end if
  if(buf(22:40) == 'ABSORPTION SPECTRUM') then
   orca6 = .true.; k = 49; exit
  end if
 end do ! for while

 if(i /= 0) then
  close(fid)
  write(6,'(/,A)') 'ERROR in subroutine read_fosc_from_orca_out: failed to read&
                   & fosc in file'
  write(6,'(A)') TRIM(outname)
  stop
 end if

 do i = 1, 4 ! skip 4 lines
  read(fid,'(A)') buf
 end do

 do i = 1, nstate, 1
  read(fid,'(A)') buf
  read(buf(k:),*) fosc(i)
 end do ! for i

 close(fid)
end subroutine read_fosc_from_orca_out

! read ADC (excitation) energies from a specified ORCA output file
subroutine read_adc_e_from_orca_out(outname, ip_or_ea, nstate, e, fosc)
 implicit none
 integer :: i, j, fid
 integer, intent(in) :: nstate
 real(kind=8), intent(out) :: e(0:nstate), fosc(nstate)
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: ip_or_ea

 buf = ' '; e = 0d0; fosc = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:12) == 'Second-order') exit
 end do ! for while

 if(i /= 0) then
  close(fid)
  write(6,'(/,A)') "ERROR in subroutine read_adc_e_from_orca_out: 'Second-order&
                   &' not found"
  write(6,'(A)') 'in file '//TRIM(outname)
  stop
 end if

 i = INDEX(buf, '...')
 read(buf(i+3:),*) e(0)

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:10) == 'ADC(2) RES') exit
 end do ! for while

 if(i /= 0) then
  close(fid)
  write(6,'(/,A)') "ERROR in subroutine read_adc_e_from_orca_out: 'ADC(2) RES' &
                   &not found in"
  write(6,'(A)') 'file '//TRIM(outname)
  stop
 end if

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:6) == 'IROOT=') then
   i = INDEX(buf, ':')
   read(buf(7:i-1),*) j
   read(buf(i+1:),*) e(j)
   if(j == nstate) exit
  end if
 end do ! for while

 close(fid)
 forall(i = 1:nstate) e(i) = e(i) + e(0)
 if(ip_or_ea) return
 call read_fosc_from_orca_out(outname, nstate, fosc)
end subroutine read_adc_e_from_orca_out

subroutine read_adc_e_from_psi4_out(outname, adc3, nstate, e, fosc)
 implicit none
 integer :: i, j, fid
 integer, intent(in) :: nstate
 real(kind=8) :: scf_e, corr_e2, corr_e3
 real(kind=8), intent(out) :: e(0:nstate), fosc(nstate)
 character(len=1) :: str1
 character(len=46), parameter :: error_warn = 'ERROR in subroutine read_adc_e_f&
                                              &rom_psi4_out: '
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: adc3

 scf_e = 0d0; corr_e2 = 0d0; corr_e3 = 0d0
 e = 0d0; fosc = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 ! read SCF energy, only the 2nd 'Total Energy' is useful
 do i = 1, 2
  do while(.true.)
   read(fid,'(A)',iostat=j) buf
   if(j /= 0) exit
   if(buf(5:16) == 'Total Energy') exit
  end do ! for while

  if(j /= 0) then
   write(6,'(/,A)') error_warn//'no SCF energy found in file'
   write(6,'(A)') TRIM(outname)
   close(fid)
   stop
  end if

  call get_dpv_after_flag(buf, '=', .true., scf_e)
 end do ! for i

 ! read MPn correlation energy
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(5:25) == 'Energy correlation MP') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'no correlation energy found in'
  write(6,'(A)') 'file '//TRIM(outname)
  close(fid)
  stop
 end if

 read(buf(27:),*) corr_e2
 if(adc3) then
  read(fid,'(A)') buf
  read(buf(27:),*) corr_e3
 end if
 e(0) = scf_e + corr_e2 + corr_e3

 ! read ADC excitation energies
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(7:22) == 'Excited states s') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'no "Excited states s" found in'
  write(6,'(A)') 'file '//TRIM(outname)
  close(fid)
  stop
 end if

 do i = 1, 6
  read(fid,'(A)') buf
 end do ! for i

 do i = 1, nstate, 1
  read(fid,*) str1, j, e(i)
 end do ! for i
 forall(i = 1:nstate) e(i) = e(i) + e(0)

 ! read oscillator strengths
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:18) == 'Excited state prop') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'no "Excited state prop" found'
  write(6,'(A)') 'in file '//TRIM(outname)
  close(fid)
  stop
 end if

 do i = 1, nstate, 1
  read(fid,'(A)') buf
  read(fid,'(A)') buf
  j = INDEX(buf, ')')
  read(buf(j+1:),*) fosc(i)
  read(fid,'(A)') buf
 end do ! for i

 close(fid)
end subroutine read_adc_e_from_psi4_out

subroutine read_adc_e_from_mrcc_out(outname, nstate, e, fosc)
 implicit none
 integer :: i, j, k, fid
 integer, intent(in) :: nstate
 real(kind=8) :: rtmp(4)
 real(kind=8), intent(out) :: e(0:nstate), fosc(nstate)
 character(len=46), parameter :: error_warn = 'ERROR in subroutine read_adc_e_f&
                                              &rom_mrcc_out: '
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical :: x_triplet

 e = 0d0; fosc = 0d0; x_triplet = .false.
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:12) == 'Total MP2 e') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"Total MP2 e" not found in file'
  write(6,'(A)') TRIM(outname)
  close(fid)
  stop
 end if

 call get_dpv_after_flag(buf, ':', .true., e(0))

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:27) == 'Final result of the ADC(2)') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"Final result of the ADC(2)" not found'
  write(6,'(A)') 'in file '//TRIM(outname)
  close(fid)
  stop
 end if

 read(fid,'(A)') buf
 select case(buf(2:8))
 case('Singlet')
 case('Triplet')
  x_triplet = .true.
 case default
  write(6,'(/,A)') error_warn//' unrecognized spin: '//buf(2:8)
  stop
 end select

 do i = 1, 4 ! skip 4 lines
  read(fid,'(A)') buf
 end do ! for i

 do i = 1, nstate, 1
  read(fid,*) j, e(i)
  e(i) = e(i) + e(0)
 end do ! for i

 ! S0->Tn transition dipoles are 0, MRCC will not calculate
 if(x_triplet) then
  close(fid)
  return
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(24:46) == 'dipole moment in length') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"dipole moment in length" not found'
  write(6,'(A)') 'in file '//TRIM(outname)
  close(fid)
  stop
 end if

 read(fid,'(A)') buf
 read(fid,'(A)') buf
 do i = 1, nstate, 1
  read(fid,*) j, k, rtmp(1), rtmp(2), rtmp(3), rtmp(4), fosc(i)
 end do ! for i

 close(fid)
end subroutine read_adc_e_from_mrcc_out

! read ADC(2) excited state energies from a specified Turbomole output file
subroutine read_adc_e_from_tm_out(outname, nstate, e, fosc)
 implicit none
 integer :: i, ntimes, fid
 integer, intent(in) :: nstate
 real(kind=8), intent(out) :: e(0:nstate), fosc(nstate)
 character(len=43), parameter :: error_warn = 'ERROR in subroutine read_adc_e_f&
                                              &rom_m_out: '
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical :: x_triplet

 e = 0d0; fosc = 0d0; x_triplet = .false.
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(10:20) == 'Final MP2 e') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"Final MP2 e" not found in file'
  write(6,'(A)') TRIM(outname)
  close(fid)
  stop
 end if

 call get_dpv_after_flag(buf, ':', .true., e(0))

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(35:53)=='ADC(2) excitation e' .and. buf(68:70)=='%t1') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"ADC(2) excitation e" and "%t1" are'
  write(6,'(A)') 'not found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 do i = 1, 3 ! skip 3 lines
  read(fid,'(A)') buf
 end do ! for i

 do i = 1, nstate, 1
  read(fid,'(A)') buf
  read(buf(25:),*) e(i)
  e(i) = e(i) + e(0)
 end do ! for i

 ntimes = 0
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(8:31) == 'oscillator strength (len') then
   ntimes = ntimes + 1
   call get_dpv_after_flag(buf, ':', .true., fosc(ntimes))
   if(ntimes == nstate) exit
  end if
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"oscillator strength (len" are'
  write(6,'(A)') 'not found in file '//TRIM(outname)
  stop
 end if
end subroutine read_adc_e_from_tm_out

subroutine read_adc_e_from_qchem_out(outname, nstate, e, fosc)
 implicit none
 integer :: i, j, k, fid
 integer, intent(in) :: nstate
 real(kind=8), intent(out) :: e(0:nstate), fosc(nstate)
 character(len=47), parameter :: error_warn = 'ERROR in subroutine read_adc_e_f&
                                              &rom_qchem_out: '
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical :: determined, x_triplet

 e = 0d0; fosc = 0d0; determined = .false.; x_triplet = .false.
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(3:15) == 'Total energy:') then
   call get_dpv_after_flag(buf, ':', .true., e(0)) ! MP2/MP3
   exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"Total energy:" not found in file'
  write(6,'(A)') TRIM(outname)
  close(fid)
  stop
 end if
 k = 1

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(3:15) == 'Excited state') then
   if(.not. determined) then
    j = INDEX(buf, '(')
    select case(buf(j+1:j+7))
    case('singlet')
     determined = .true.
    case('triplet')
     determined = .true.
     x_triplet = .true.
    case default
     write(6,'(/,A)') error_warn//'state spin cannot be recognized.'
     write(6,'(A)') 'buf='//TRIM(buf)
     stop
    end select
   end if
   do j = 1, 4
    read(fid,'(A)') buf
   end do ! for j
   call get_dpv_after_flag(buf, ':', .true., e(k))
   if(.not. x_triplet) then
    do j = 1, 3
     read(fid,'(A)') buf
    end do ! for j
    call get_dpv_after_flag(buf, ':', .true., fosc(k))
   end if
   k = k + 1
   if(k == nstate+1) exit
  end if
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') error_warn//'not enough states are found in file'
  write(6,'(A)') TRIM(outname)
  stop
 end if
end subroutine read_adc_e_from_qchem_out

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

! Read CASCI/CASSCF energy from a Gaussian .log file. In an geometry optimization
! or numerical differentiation job, there might be multiple frames.
! first=.True. : try to locate CASCI/CASSCF energies of the 1st frame
!      =.False.: try to locate CASCI/CASSCF energies of the last frame
subroutine read_cas_energy_from_gau_log(logname, scf, first, e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: e(2)
 character(len=50), parameter :: error_warn = 'ERROR in subroutine read_cas_ene&
                                              &rgy_from_gau_log: '
 character(len=240) :: buf
 character(len=240), intent(in) :: logname
 logical, intent(in) :: scf, first

 if(first) then
  open(newunit=fid,file=TRIM(logname),status='old',position='rewind')
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(13:22)=='EIGENVALUE' .or. buf(20:29)=='EIGENVALUE' .or. &
      buf(23:32)=='Eigenvalue') exit
  end do ! for while
 else
  open(newunit=fid,file=TRIM(logname),status='old',position='append')
  do while(.true.)
   BACKSPACE(fid)
   BACKSPACE(fid)
   read(fid,'(A)') buf
   if(buf(2:18) == 'Entering Gaussian') then
    i = -1; exit
   end if
   if(buf(13:22)=='EIGENVALUE' .or. buf(20:29)=='EIGENVALUE' .or. &
      buf(23:32)=='Eigenvalue') then
    i = 0; exit
   end if
  end do ! for while
 end if

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"EIGENVALUE"/"Eigenvalue" not'
  write(6,'(A)') 'located in file '//TRIM(logname)
  close(fid)
  stop
 end if

 i = INDEX(buf, 'lue')
 if(i == 0) i = INDEX(buf, 'LUE')

 if(scf) then
  read(buf(i+3:),*) e(2) ! CASSCF energy
  do while(.true.)
   BACKSPACE(fid)
   BACKSPACE(fid)
   read(fid,'(A)') buf
   if(buf(2:18) == 'Entering Gaussian') then
    i = -1; exit
   end if
   if(buf(2:8) == 'ITN=  1') then
    i = 0; exit
   end if
  end do ! for while
  close(fid)
  if(i /= 0) then
   write(6,'(/,A)') error_warn//'"ITN=  1" not located in file'
   write(6,'(A)') TRIM(logname)
   stop
  end if
  i = INDEX(buf, 'E=')
  read(buf(i+2:),*) e(1) ! CASCI energy in CASSCF job
 else
  read(buf(i+3:),*) e(1) ! CASCI energy in CASCI job
 end if
end subroutine read_cas_energy_from_gau_log

! read (DMRG-)CASCI energy from a PySCF output file
! DO NOT use this subroutine to read (DMRG-)CASSCF energy.
subroutine read_casci_e_from_pyscf_out(outname, target_mult, dmrg, e)
 implicit none
 integer :: i, fid
 integer, intent(in) :: target_mult ! target spin multiplicity
 real(kind=8) :: ssquare, expected
 real(kind=8), intent(out) :: e(2)
 real(kind=8), parameter :: ss_diff_thres = 1d-3
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 character(len=49), parameter :: err_str = 'ERROR in subroutine read_casci_e_fr&
                                           &om_pyscf_out: '
 logical, intent(in) :: dmrg
 logical :: state_specific

 e = 0d0; ssquare = 0d0; state_specific = .false.
 expected = 0.25d0*DBLE(target_mult*target_mult - 1)

 if(dmrg) then ! DMRG-CASCI
  open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
 else          ! CASCI
  open(newunit=fid,file=TRIM(outname),status='old',position='append')
  do while(.true.)
   BACKSPACE(fid)
   BACKSPACE(fid)
   read(fid,'(A)') buf
   if(buf(1:15) == 'CASCI converged') exit
   if(buf(1:19) == 'CASCI not converged') then
    write(6,'(/,A)') err_str//'it seems that CASCI is not'
    write(6,'(A)') 'converged. Please check file '//TRIM(outname)
    call prt_pyscf_forge_recommend()
    close(fid)
    stop
   end if
  end do ! for while
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:7) == 'CASCI E') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') err_str//'`CASCI E` keyword not located'
  write(6,'(A)') 'after `CASCI converged`.'
  stop
 end if
 call get_dpv_after_flag(buf, '=', .true., e(1))

 i = INDEX(buf, 'S^2 =', back=.true.)
 if(i == 0) then
  write(6,'(/,A)') err_str//'`S^2 =` not located in the desired'
  write(6,'(A)') 'line. buf=`'//TRIM(buf)//'`'
  write(6,'(A)') 'Probably your PySCF version is so old that `S^2 =` is not pri&
                 &nted.'
  stop
 end if
 read(buf(i+5:),*) ssquare

 if(DABS(expected - ssquare) > ss_diff_thres) then
  write(6,'(/,A)') err_str//'it seems that CASCI is converged'
  write(6,'(A)') 'to a wrong spin state.'
  call prt_pyscf_forge_recommend()
  stop
 end if
end subroutine read_casci_e_from_pyscf_out

! read (DMRG-)CASSCF energy from a PySCF output file
! DO NOT use this subroutine to read (DMRG-)CASCI energy.
subroutine read_casscf_e_from_pyscf_out(outname, target_mult, dmrg, e)
 implicit none
 integer :: i, j, k, fid
 integer, intent(in) :: target_mult ! target spin multiplicity
 real(kind=8) :: ssquare, expected
 real(kind=8), intent(out) :: e(2)
 real(kind=8), parameter :: ss_diff_thres = 1d-3
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 character(len=50), parameter :: err_str = 'ERROR in subroutine read_casscf_e_f&
                                           &rom_pyscf_out: '
 logical, intent(in) :: dmrg
 logical :: state_specific

 e = 0d0; i = 0; j = 0; k = 0; ssquare = 0d0; state_specific = .false.
 expected = 0.25d0*DBLE(target_mult*target_mult - 1)
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

 ! read CASSCF energy in a CASSCF job
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:7) == 'CASCI E') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') err_str//"'CASCI E' keyword not found in"
  write(6,'(A)') 'file '//TRIM(outname)
  close(fid)
  stop
 end if
 call get_dpv_after_flag(buf, '=', .true., e(2))

 i = INDEX(buf, 'S^2 =', back=.true.)
 if(i == 0) then
  write(6,'(/,A)') err_str//'`S^2 =` not located in the desired line.'
  write(6,'(A)') 'buf=`'//TRIM(buf)//'`'
  write(6,'(A)') 'Probably your PySCF version is so old that `S^2 =` is not pri&
                 &nted.'
  close(fid)
  stop
 end if
 read(buf(i+5:),*) ssquare

 if(DABS(expected - ssquare) > ss_diff_thres) then
  write(6,'(/,A)') REPEAT('-',79)
  write(6,'(A)') 'Warning from subroutine read_cas_energy_from_pyout: <S**2> de&
                 &viates too much'
  write(6,'(A,F8.3,A,F11.5)') 'from the expectation value. Expected <S^2>=', &
                              expected, ', Calc. <S^2>=', ssquare
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
 rewind(fid)
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:9) == 'CASCI E =') exit
 end do ! for while

 close(fid)
 read(buf(10:),*) e(1)
 call get_dpv_after_flag(buf, '=', .false., ssquare)

 if(DABS(expected - ssquare) > ss_diff_thres) then
  write(6,'(/,A)') REPEAT('-',79)
  write(6,'(A)') 'Warning in subroutine read_cas_energy_from_pyout: the 0-th s&
                 &tep in this CASSCF'
  write(6,'(A)') 'job, i.e. the CASCI <S**2> deviates too much from the expect&
                 &ation value.'
  write(6,'(A,F8.3,A,F11.5)') 'Expected <S^2>=',expected,', Calc. <S^2>=',ssquare
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
end subroutine read_casscf_e_from_pyscf_out

! read CASCI/CASSCF energy from the GAMESS output file
subroutine read_cas_energy_from_gms_gms(outname, target_mult, scf, e)
 implicit none
 integer :: i, fid
 integer, intent(in) :: target_mult
 real(kind=8), parameter :: ss_diff_thres = 1d-3
 real(kind=8) :: ssquare, expected
 real(kind=8), intent(out) :: e(2)
 character(len=41), parameter :: error_warn = 'subroutine read_cas_energy_from_&
                                              &gms_gms: '
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: scf

 expected = 0.25d0*DBLE(target_mult*target_mult - 1)
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 if(scf) then  ! CASSCF job
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(2:24) == 'THE DENSITIES ARE STATE') exit
  end do ! for while
 
  if(i /= 0) then
   write(6,'(/,A)') error_warn//'no `THE DENSITIES ARE STATE` found'
   write(6,'(A)') 'in file '//TRIM(outname)
   close(fid)
   stop
  end if

  read(fid,'(A)') buf
  i = INDEX(buf,'ENERGY=')
  read(buf(i+7:),*) e(1)   ! CASCI energy in the CASSCF job
  i = INDEX(buf,'=',back=.true.)
  read(buf(i+1:),*) ssquare
  ssquare = ssquare*(ssquare+1d0)
  if(DABS(expected - ssquare) > ss_diff_thres) then
   write(6,'(/,A)') error_warn//'in this CASSCF job, the 0-th'
   write(6,'(A)') 'step, i.e., the CASCI <S**2> deviates too much from the expe&
                  &ctation value.'
   write(6,'(A,F8.3,A,F11.5)') 'Expected <S^2>=', expected, ', Calc. <S^2>=', &
                               ssquare
   stop
  end if

  do while(.true.)
   read(fid,'(A)') buf
   if(buf(2:13) == 'STATE   1  E') exit
  end do ! for while
  i = INDEX(buf,'ENERGY=')
  read(buf(i+7:),*) e(2)   ! CASSCF energy
  i = INDEX(buf, 'S=')
  read(buf(i+2:),*) ssquare
  ssquare = ssquare*(ssquare+1d0)
  if(DABS(expected - ssquare) > ss_diff_thres) then
   write(6,'(/,A)') error_warn//'CASSCF <S**2> deviates too'
   write(6,'(A,F8.3,A,F11.5)') 'much from the expectation value. Expected <S^2>&
                               &= ', expected, ', Calc. <S^2>=', ssquare
   stop
  end if

 else          ! CASCI job
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(2:20) == 'DENSITY MATRIX WILL') exit
  end do ! for while
 
  if(i /= 0) then
   write(6,'(/,A)') error_warn//'no `DENSITY MATRIX WILL` found'
   write(6,'(A)') 'in file '//TRIM(outname)
   stop
  end if
 
  i = INDEX(buf,'=', back=.true.)
  read(buf(i+1:),*) ssquare
  ssquare = ssquare*(ssquare+1d0)
  if(DABS(expected - ssquare) > ss_diff_thres) then
   write(6,'(/,A)') error_warn//'CASCI <S^2> deviates too'
   write(6,'(A,F8.3,A,F11.5)') 'much from the expectation value. Expected <S^2>&
                               &= ', expected, ', Calc. <S^2>=', ssquare
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
  call get_dpv_after_flag(buf, '=', .true., e(1))

  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:19) == 'Final CASSCF energy') exit
  end do ! for while
 
  if(i /= 0) then
   write(6,'(/,A)') error_warn//'"Final CASSCF energy" not found in'
   write(6,'(A)') 'file '//TRIM(outname)
   close(fid)
   stop
  end if
  call get_dpv_after_flag(buf, ':', .true., e(2))

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
   write(6,'(/,A)') error_warn//'"^STATE" not found in file '//TRIM(outname)
   close(fid)
   stop
  end if

  call get_dpv_after_flag(buf, '=', .true., e(1))
 end if

 close(fid)
end subroutine read_cas_energy_from_orca_out

subroutine find_mcscf_e_in_molpro_out(outname, casscf_e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: casscf_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical :: alive

 casscf_e = 0d0
 inquire(file=TRIM(outname),opened=alive)
 if(alive) then
  inquire(file=TRIM(outname),number=fid)
 else
  open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:9) == '!MCSCF S') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine find_mcscf_e_in_molpro_out: "!MCSCF S" &
                   &not found in'
  write(6,'(A)') 'file '//TRIM(outname)
  stop
 end if

 i = INDEX(buf, 'Energy')
 read(buf(i+6:),*) casscf_e
end subroutine find_mcscf_e_in_molpro_out

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
  case(' Iter  ') ! new output format in Molpro CASCI
   call find_mcscf_e_in_molpro_out(outname, e(1))
   return
  end select
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_cas_energy_from_molpro_out: CASCI &
                   &energy string not'
  write(6,'(A)') 'located in file '//TRIM(outname)
  close(fid)
  stop
 end if

 do i = 1, 4
  read(fid,fmt=*,iostat=j) k,k,k,k, e(1) ! CASCI energy
  if(j == 0) exit
 end do ! for i

 if(scf) call find_mcscf_e_in_molpro_out(outname, e(2))
 ! The file `outname` will be closed in subroutine find_mcscf_e_in_molpro_out,
 ! so there is no need to close it here.
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
 character(len=53), parameter :: error_warn = 'ERROR in subroutine read_cas_ene&
                                              &rgy_from_dalton_out: '
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
  write(6,'(/,A)') error_warn//'"@ Final CI energies" not'
  write(6,'(A)') 'found in file '//TRIM(outname)
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
  write(6,'(/,A)') error_warn//'"@    Final MCSCF en" not'
  write(6,'(A)') 'found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf,':')
 read(buf(i+1:),*) e(2)
 close(fid)
end subroutine read_cas_energy_from_dalton_out

! read CASCI/CASSCF energy from a given MRCC output file
subroutine read_cas_energy_from_mrcc_out(outname, scf, e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: e(2)
 character(len=51), parameter :: error_warn = 'ERROR in subroutine read_cas_ene&
                                              &rgy_from_mrcc_out: '
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: scf

 e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(6:29) == 'MCSCF ENERGY IN STEP   1') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'`MCSCF ENERGY IN STEP   1`'
  write(6,'(A)') 'not found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 read(buf(33:),*) e(1) ! CASCI energy
 if(.not. scf) then
  close(fid)
  return
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(5:23) == 'FINAL MCSCF ENERGY:') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') error_warn//'`FINAL MCSCF ENERGY:` not'
  write(6,'(A)') 'found in file '//TRIM(outname)
  stop
 end if

 call get_dpv_after_flag(buf, ':', .false., e(2))
end subroutine read_cas_energy_from_mrcc_out

! read CASCI/CASSCF energy from a Gaussian/PySCF/GAMESS/OpenMolcas/ORCA output file
subroutine read_cas_energy_from_output(cas_prog, outname, target_mult, scf, &
                                       dmrg, ptchg_e, nuc_pt_e, e)
 implicit none
 integer, intent(in) :: target_mult
 real(kind=8), intent(in) :: ptchg_e, nuc_pt_e
 real(kind=8), intent(out) :: e(2)
 character(len=10), intent(in) :: cas_prog
 character(len=240), intent(in) :: outname
 logical, intent(in) :: scf, dmrg

 select case(TRIM(cas_prog))
 case('gaussian')
  call read_cas_energy_from_gau_log(outname, scf, .false., e)
 case('gamess')
  call read_cas_energy_from_gms_gms(outname, target_mult, scf, e)
  e(1) = e(1) + ptchg_e + nuc_pt_e
 case('pyscf')
  if(scf) then
   call read_casscf_e_from_pyscf_out(outname, target_mult, dmrg, e)
  else
   call read_casci_e_from_pyscf_out(outname, target_mult, dmrg, e)
  end if
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
 case('mrcc')
  call read_cas_energy_from_mrcc_out(outname, scf, e)
  e = e + ptchg_e
 case default
  write(6,'(/,A)') 'ERROR in subroutine read_cas_energy_from_output: CAS_prog c&
                   &annot be'
  write(6,'(A)') 'identified. CAS_prog='//TRIM(cas_prog)
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
subroutine read_mrpt_energy_from_pyscf_out(outname, target_mult, target_root, &
                                           ref_e, corr_e)
 implicit none
 integer :: i, k, fid
 integer, intent(in) :: target_mult, target_root
 ! target_mult: target spin multiplicity
 ! target_root: 0 for the ground state, >0 for excited state
 real(kind=8) :: ssquare0, ssquare
 real(kind=8), parameter :: ss_diff_thres = 1d-3
 real(kind=8), intent(out) :: ref_e, corr_e
 character(len=53), parameter :: error_warn = 'ERROR in subroutine read_mrpt_en&
                                              &ergy_from_pyscf_out: '
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 ssquare0 = 0.25d0*DBLE(target_mult*target_mult - 1)
 ssquare = 0d0; ref_e = 0d0; corr_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  if(buf(1:7) == 'SCF ene') then
   i = -1; exit
  end if
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:8) == 'Nevpt2 E') then
   call get_dpv_after_flag(buf, '=', .true., corr_e)
   exit
  end if
  if(buf(1:17) == 'E(WickICNEVPT2) =') then
   call get_dpv_after_flag(buf, '=', .false., corr_e)
   exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'neither `Nevpt2 E` nor `E(WickICNEVPT2) =`'
  write(6,'(A)') 'is found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 if(target_root == 0) then
  do while(.true.)
   BACKSPACE(fid)
   BACKSPACE(fid)
   if(buf(1:7) == 'SCF ene') then
    i = -1; exit
   end if
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:7) == 'CASCI E') exit
  end do ! for while
 else
  do while(.true.)
   BACKSPACE(fid)
   BACKSPACE(fid)
   if(buf(1:7) == 'SCF ene') then
    i = -1; exit
   end if
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:11) == 'CASCI state') then
    read(buf(12:),*) k
    if(k == target_root) exit
   end if
  end do ! for while
 end if

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'no CASCI energy found in file'
  write(6,'(A)') TRIM(outname)
  write(6,'(A,I0)') 'target_root=', target_root
  close(fid)
  stop
 end if

 call get_dpv_after_flag(buf, '=', .true., ref_e)
 call get_dpv_after_flag(buf, '=', .false., ssquare)

 if(DABS(ssquare0 - ssquare) > ss_diff_thres) then
  write(6,'(/,A,F8.3,A,F10.3)') 'Expected <S^2>=',ssquare0,', Calc. <S^2>=', &
                                ssquare
  write(6,'(A)') error_warn//'it seems that CASCI is'
  write(6,'(A)') 'converged to a different spin state. Please check file'
  write(6,'(A)') TRIM(outname)
  call prt_pyscf_forge_recommend()
  close(fid)
  stop
 end if

 BACKSPACE(fid)
 BACKSPACE(fid)
 read(fid,'(A)') buf
 close(fid)

 if(buf(1:14) == 'CASCI not conv') then
  write(6,'(/,A)') error_warn//'CASCI Davidson iteration is'
  write(6,'(A)') 'not converged in file '//TRIM(outname)
  call prt_pyscf_forge_recommend()
  stop
 end if
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
  call get_dpv_after_flag(buf, '=', .true., ref_e)

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

  call get_dpv_after_flag(buf, ':', .true., corr_e)
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

 call get_dpv_after_flag(buf, '=', .true., ref_e)

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(22:42) == '2ND ORDER ENERGY CORR') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_mrpt_energy_from_gms_out: no MRMP2&
                   & energy found in'
  write(6,'(A)') 'file '//TRIM(outname)
  stop
 end if

 call get_dpv_after_flag(buf, '=', .true., corr_e)
end subroutine read_mrpt_energy_from_gms_out

! read NEVPT2 energy from a BDF .out file
subroutine read_mrpt_energy_from_bdf_out(outname, itype, ref_e, corr_e, dav_e)
 implicit none
 integer :: i, fid
 integer, intent(in) :: itype   ! 1/2 for SDSPT2/FIC-NEVPT2
 character(len=51), parameter :: error_warn = 'ERROR in subroutine read_mrpt_en&
                                              &ergy_from_bdf_out: '
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8) :: rtmp(6)
 real(kind=8), intent(out) :: ref_e, corr_e, dav_e
 ! dav_e: Davidson correction energy for size-inconsistency error

 ref_e = 0d0; corr_e = 0d0; dav_e = 0d0
 if(.not. (itype==1 .or. itype==2)) then
  write(6,'(/,A)') error_warn//'currently only reading'
  write(6,'(A)') 'SDSPT2/NEVPT2 energy is supported.'
  stop
 end if

 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:11) == 'Print final') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'"Print final" not found in'
  write(6,'(A)') 'file '//TRIM(outname)
  close(fid)
  stop
 end if

 do i = 1, 4
  read(fid,'(A)') buf
 end do ! for i

 if(INDEX(buf, '=') > 0) then
  call get_dpv_after_flag(buf, '=', .true.,  ref_e) ! CASCI/CASSCF energy
 else
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:12) == 'TARGET_MCSCF') exit
  end do ! for while
  if(i /= 0) then
   write(6,'(/,A)') error_warn//'no "TARGET_MCSCF" found in'
   write(6,'(A)') 'file '//TRIM(outname)
   close(fid)
   stop
  end if
  read(buf(26:),*) ref_e
 end if

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
  write(6,'(/,A)') error_warn//'no "TARGET_XIANCI" found in'
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

! read Davidson correction and MRCISD energy from ORCA output file
subroutine read_mrcisd_e_from_orca_out(CtrType, outname, davidson_e, e)
 implicit none
 integer :: i, fid
 integer, intent(in) :: CtrType
 real(kind=8) :: casci_e, ref_weight
 real(kind=8), intent(out) :: davidson_e, e
 character(len=240), intent(in) :: outname
 character(len=240) :: buf
 logical, external :: has_davidson_q

 davidson_e = 0d0; e = 0d0; ref_weight = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 select case(CtrType)
 case(1) ! uncontracted MRCISD
  do while(.true.)
   BACKSPACE(fid)
   BACKSPACE(fid)
   read(fid,'(A)') buf
   if(buf(36:48) == 'O   R   C   A') then
    write(6,'(/,A)') 'ERROR in subroutine read_mrcisd_e_from_orca_out: failed to &
                     &read MRCISD'
    write(6,'(A)') 'energy from file '//TRIM(outname)
    close(fid)
    stop
   end if
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
  call get_dpv_after_flag(buf, '=', .false., ref_weight)
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

 case(3) ! FIC-MRCISD
  do while(.true.)
   BACKSPACE(fid)
   BACKSPACE(fid)
   read(fid,'(A)') buf
   if(buf(36:48) == 'O   R   C   A') then
    write(6,'(/,A)') 'ERROR in subroutine read_mrcisd_e_from_orca_out: failed to &
                     &read MRCISD'
    write(6,'(A)') 'energy from file '//TRIM(outname)
    close(fid)
    stop
   end if
   if(buf(1:5) == 'ETOT ') exit
  end do ! for while

  i = INDEX(buf,'...', back=.true.)
  read(buf(i+3:),*) e
  read(fid,'(A)') buf
  read(fid,'(A)') buf
  i = INDEX(buf,'...', back=.true.)
  read(buf(i+3:),*) davidson_e
 case default
  close(fid)
  write(6,'(/,A,I0)') 'ERROR in subroutine read_mrcisd_e_from_orca_out: invalid&
                      & CtrType=', CtrType
  write(6,'(A)') 'outname='//TRIM(outname)
  stop
 end select

 close(fid)
end subroutine read_mrcisd_e_from_orca_out

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
   call get_dpv_after_flag(buf, '=', .true., e)
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
   call get_dpv_after_flag(buf, ':', .false., e)
  else if(CtrType == 2) then ! ic-MRCISD
   do while(.true.)
    BACKSPACE(fid)
    BACKSPACE(fid)
    read(fid,'(A)') buf
    if(buf(23:31) == 'CI ENERGY') exit
   end do ! for while
   call get_dpv_after_flag(buf, ':', .false., e)
   read(fid,'(A)') buf
   call get_dpv_after_flag(buf, ':', .false., davidson_e)
  end if
  e = e + ptchg_e

 case('orca')
  close(fid)
  call read_mrcisd_e_from_orca_out(CtrType, outname, davidson_e, e)

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
    write(6,'(/,A)') 'ERROR in subroutine read_mrci_energy_from_output: no'
    write(6,'(A)') '"Total CI energy" found in file '//TRIM(outname)
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
    write(6,'(/,A)') 'ERROR in subroutine read_mrci_energy_from_output: no'
    write(6,'(A)') '"@ Final CI energies" found in file '//TRIM(outname)
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
   call get_dpv_after_flag(buf, '=', .true., casci_e)

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
   call get_dpv_after_flag(buf, '=', .true., e)
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

  call get_dpv_after_flag(buf, '=', .true., ref_e)

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
  ! GVB-BCCCnb correlation energy
  call get_dpv_after_flag(buf, '=', .true., corr_e(1))

  read(fid,'(A)') buf
  ! GVB-BCCCnb total energy
  call get_dpv_after_flag(buf, '=', .true., corr_e(2))
  close(fid)

 !case('nwchem')
 case default
  write(6,'(/,A)') 'ERROR in subroutine read_mrcc_energy_from_output: invalid m&
                   &rcc_prog='//TRIM(mrcc_prog)
  stop
 end select
end subroutine read_mrcc_energy_from_output

subroutine read_cc_e_from_gau_log(outname, t1diag, ref_e, tot_e)
 implicit none
 integer :: i, fid
 real(kind=8) :: ss
 real(kind=8), intent(out) :: t1diag, ref_e, tot_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 t1diag = 0d0; ref_e = 0d0; tot_e = 0d0
 call read_hf_e_and_ss_from_gau_log(outname, .false., ref_e, ss)

 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:17) == 'Wavefunction amp') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_orca_out: 'Wavefunction &
                   &amp' not found"
  write(6,'(A)') 'in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, '=')
 read(buf(i+1:),*) tot_e ! CCD/CCSD total energy

 read(fid,'(A)') buf
 if(buf(2:8) == 'T1 Diag') then
  i = INDEX(buf, '=')
  read(buf(i+1:),*) t1diag
 end if

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:14) == 'Discarding MO') exit
  if(buf(2:8)=='CCSD(T)' .or. buf(2:9)=='QCISD(T)') then
   i = INDEX(buf, '=')
   read(buf(i+1:),*) tot_e ! CCSD(T)/QCISD(T) total energy
   exit
  end if
 end do ! for while

 close(fid)
end subroutine read_cc_e_from_gau_log

subroutine read_cc_e_from_molpro_out(outname, read_t1diag, t1diag, ref_e, tot_e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: t1diag, ref_e, tot_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: read_t1diag

 t1diag = 0d0; ref_e = 0d0; tot_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(2:9) == 'Checking') then
   i = -1
   exit
  end if
  if(INDEX(buf,'HF-SCF') > 0) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_orca_out: 'Total Ener' n&
                   &ot found in file "//TRIM(outname)
  close(fid)
  stop
 end if

 read(fid,*) tot_e, ref_e

 if(read_t1diag) then
  do while(.true.)
   BACKSPACE(fid,iostat=i)
   if(i /= 0) exit
   BACKSPACE(fid,iostat=i)
   if(i /= 0) exit
   read(fid,'(A)') buf
   if(buf(2:9) == 'Checking') then
    i = -1
    exit
   end if
   if(INDEX(buf,'T1 diag') > 0) exit
  end do ! for while
  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_orca_out: 'T1 diag' not&
                   & found in file "//TRIM(outname)
   close(fid)
   stop
  end if
  i = INDEX(buf, 'T1 diagnostic:')
  read(buf(i+14:),*) t1diag
 end if

 close(fid)
end subroutine read_cc_e_from_molpro_out

subroutine read_cc_e_from_pyscf_out(outname, t1diag, ref_e, tot_e)
 implicit none
 integer :: i, fid
 real(kind=8) :: ccsd_t_corr
 real(kind=8), intent(out) :: t1diag, ref_e, tot_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical :: set_tot_e, set_ccsd_t_corr

 set_tot_e = .false.; set_ccsd_t_corr = .false.
 t1diag = 0d0; ref_e = 0d0; tot_e = 0d0; ccsd_t_corr = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')
 write(fid,'(/)',advance='no')
 i = 0

 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(buf(1:12) == 'SCF not conv') then
   i = -1; exit
  end if

  select case(buf(1:11))
  case('E(RCCSDT) =','E(UCCSDT) =')
   if(.not. set_tot_e) then
    read(buf(12:),*) tot_e
    set_tot_e = .true.
   end if
  end select

  select case(buf(1:10))
  case('RCCSD(T) c','UCCSD(T) c')
   if(.not. set_ccsd_t_corr) then
    call get_dpv_after_flag(buf, '=', .true., ccsd_t_corr)
    set_ccsd_t_corr = .true.
   end if
  case('E(RCCSD) =','E(UCCSD) =')
   if(.not. set_tot_e) then
    read(buf(11:),*) tot_e
    set_tot_e = .true.
   end if
  end select

  select case(buf(1:9))
  case('CCSD(T) c')
   if(.not. set_ccsd_t_corr) then
    call get_dpv_after_flag(buf, '=', .true., ccsd_t_corr)
    set_ccsd_t_corr = .true.
   end if
  case('E(CCSD) =')
   if(.not. set_tot_e) then
    read(buf(10:),*) tot_e
    set_tot_e = .true.
   end if
  case('T1_diag =')
   read(buf(10:),*) t1diag
  case('converged')
   exit
  end select
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_cc_e_from_pyscf_out: failed to rea&
                   &d HF and CC'
  write(6,'(A)') 'energies from file '//TRIM(outname)
  stop
 end if

 call get_dpv_after_flag(buf, '=', .true., ref_e)
 tot_e = tot_e + ccsd_t_corr
end subroutine read_cc_e_from_pyscf_out

! Read CCSD(T) energy from a PSI4 output file. Usually used for the job:
!  grad, wfn = gradient('ccsd(t)', return_wfn=True)
subroutine read_ccsd_t_e_from_psi4_grad_out(outname, RI, t1diag, ref_e, tot_e)
 implicit none
 integer :: i, fid, ntimes
 real(kind=8) :: r1, r2
 real(kind=8), intent(out) :: t1diag, ref_e, tot_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: RI

 t1diag = 0d0; ref_e = 0d0; tot_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 ! read SCF energy
 ntimes = 0
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(5:18) == 'Total Energy =') then
   ntimes = ntimes + 1
   if(ntimes == 2) exit
  end if
 end do ! for while
 read(buf(19:),*) ref_e

 if(RI) then
  ! read T1diag
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(INDEX(buf,'T1 diagnostic') > 0) exit
  end do ! for while
  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_ccsd_t_e_from_psi4_grad_out: no '&
                    &T1 diagnostic' found"
   write(6,'(A)') 'in file '//TRIM(outname)
   close(fid)
   stop
  end if
  read(fid,'(A)') buf
  i = INDEX(buf, ':')
  read(buf(i+1:),*) t1diag

  ! read CCSD(T) energy
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(INDEX(buf,'DF-CCSD(T) Tot') > 0) exit
  end do ! for while
  close(fid)

  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_ccsd_t_e_from_psi4_grad_out: no '&
                    &DF-CCSD(T) Tot'"
   write(6,'(A)') 'found in file '//TRIM(outname)
   stop
  end if
  i = INDEX(buf, ':', back=.true.)

 else ! no RI
  ! read T1diag
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(51:56) == 'T1Diag') exit
  end do ! for while

  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_ccsd_t_e_from_psi4_grad_out: no 'T&
                    &1Diag' found in"
   write(6,'(A)') 'file '//TRIM(outname)
   close(fid)
   stop
  end if

  read(fid,'(A)') ! skip 1 line
  do while(.true.)
   read(fid,'(A)') buf
   if(LEN_TRIM(buf) == 0) exit
   read(buf,*) i, r1, r2, t1diag
  end do ! for while

  ! read CCSD(T) energy
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(9:23) == 'CCSD(T) total e') exit
  end do ! for while
  close(fid)
  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_ccsd_t_e_from_psi4_grad_out: no 'C&
                    &CSD(T) total e'"
   write(6,'(A)') 'found in file '//TRIM(outname)
   stop
  end if
  i = INDEX(buf, '=', back=.true.)
 end if

 read(buf(i+1:),*) tot_e
end subroutine read_ccsd_t_e_from_psi4_grad_out

subroutine read_cc_e_from_psi4_out(outname, RI, ccd, ccsd, ccsd_t, t1diag, ref_e,&
                                   tot_e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: t1diag, ref_e, tot_e
 character(len=14), parameter :: key0(5) = ['DF-CCD Total E','* CCSD total e',&
                           'DF-CCSD Total ','* CCSD(T) tota','DF-CCSD(T) Tot']
 character(len=14) :: key
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: RI, ccd, ccsd, ccsd_t

 t1diag = 0d0; ref_e = 0d0; tot_e = 0d0

 if(ccd) then
  key = key0(1)
 else if(ccsd) then
  if(RI) then
   key = key0(3)
  else
   key = key0(2)
  end if
 else if(ccsd_t) then
  if(RI) then
   key = key0(5)
  else
   key = key0(4)
  end if
 end if

 open(newunit=fid,file=TRIM(outname),status='old',position='append')
 if(.not. ccd) then
  do while(.true.)
   BACKSPACE(fid, iostat=i)
   if(i /= 0) exit
   BACKSPACE(fid, iostat=i)
   if(i /= 0) exit
   read(fid,'(A)') buf
   if(buf(1:6) == 'memory') then
    i = -1; exit
   end if
   if(INDEX(buf,'T1 DIAG')>0 .or. INDEX(buf,'T1 diag')>0) exit
  end do ! for while
  i = INDEX(buf, '=>')
  if(i == 0) i = INDEX(buf, ':')
  read(buf(i+2:),*) t1diag
 end if

 rewind(fid)
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(INDEX(buf, key) > 0) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_psi4_out: key '"//key//&
                   "'"
  write(6,'(A)') 'not found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, '=')
 if(i == 0) i = INDEX(buf, ':')
 read(buf(i+1:),*) tot_e

 do while(.true.)
  BACKSPACE(fid, iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid, iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(1:6) == 'memory') then
   i = -1; exit
  end if
  if(INDEX(buf,'SCF E')>0 .or. INDEX(buf,'SCF e')>0) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_psi4_out: 'SCF E' or 'SC&
                   &F e' not found"
  write(6,'(A)') 'in file '//TRIM(outname)
  close(fid)
  stop
 end if

 close(fid)
 i = INDEX(buf, '=')
 if(i == 0) i = INDEX(buf, ':')
 read(buf(i+1:),*) ref_e
end subroutine read_cc_e_from_psi4_out

subroutine read_scf_e_from_gms_out(outname, scf_e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: scf_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 scf_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:7) == ' FINAL ') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_scf_e_from_gms_out: ' FINAL ' not &
                   &found in file "//TRIM(outname)
  stop
 end if

 i = INDEX(buf, 'IS')
 buf = ADJUSTL(buf(i+2:))
 read(buf,*) scf_e
end subroutine read_scf_e_from_gms_out

! read CC energy from GAMESS output
subroutine read_cc_e_from_gms_out(outname, ccd, ccsd, ccsd_t, t1diag, ref_e, &
                                  tot_e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: t1diag, ref_e, tot_e
 character(len=10), parameter :: key0(3) = ['CCD ENERGY','CCSD    EN','CCSD(T) EN']
 character(len=10) :: key
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: ccd, ccsd, ccsd_t

 t1diag = 0d0; ref_e = 0d0; tot_e = 0d0
 if(ccd) then
  key = key0(1)
 else if(ccsd) then
  key = key0(2)
 else if(ccsd_t) then
  key = key0(3)
 end if

 call read_scf_e_from_gms_out(outname, ref_e)
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 if(.not. ccd) then
  do while(.true.)
   read(fid,'(A)') buf
   if(buf(2:8) == 'T1 DIAG') exit
  end do ! for while
  i = INDEX(buf, '=')
  read(buf(i+1:),*) t1diag
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(INDEX(buf, key) > 0) exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_gms_out: key '"//key//"'"
  write(6,'(A)') 'not found in file '//TRIM(outname)
  stop
 end if

 i = INDEX(buf, ':')
 read(buf(i+1:),*) tot_e
end subroutine read_cc_e_from_gms_out

subroutine read_cc_e_from_qchem_out(outname, ccd, ccsd, ccsd_t, t1diag, ref_e, tot_e)
 implicit none
 integer :: i, na, nb, fid
 real(kind=8), intent(out) :: t1diag, ref_e, tot_e
 character(len=10) :: key
 character(len=10), parameter :: key0(3) = ['CCD total ','CCSD total','CCSD(T) to']
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: ccd, ccsd, ccsd_t

 t1diag = 0d0; ref_e = 0d0;  tot_e = 0d0
 if(ccd) then
  key = key0(1)
 else if(ccsd) then
  key = key0(2)
 else if(ccsd_t) then
  key = key0(3)
 end if

 open(newunit=fid,file=TRIM(outname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid, iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid, iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(2:8) == 'SCF ene') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_qchem_out: 'SCF ene' not&
                   & found in file "//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, '=')
 read(buf(i+1:),*) ref_e ! SCF energy

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:11) == key) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_qchem_out: key '"//key//"'"
  write(6,'(A)') 'not found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, '=')
 read(buf(i+1:),*) tot_e

 if(.not. ccd) then ! read T1diag
  read(fid,'(A)') buf
  read(fid,'(A)') buf
  i = INDEX(buf, '=')
  read(buf(i+1:),*) t1diag
  if(t1diag < 5d-3) then
   write(6,'(A)') REPEAT('-',79)
   write(6,'(A)') 'Warning from subroutine read_cc_e_from_qchem_out: T1^2 in Q-&
                  &Chem output file'
   write(6,'(A)') 'has too few digits. This does not affect the electronic ener&
                  &gy at all, but will'
   write(6,'(A)') 'lead to inaccurate T1diag value. If you care about T1diag, y&
                  &ou should use'
   write(6,'(A)') 'another CC_prog.'
   write(6,'(A)') REPEAT('-',79)
  end if
  do while(.true.)
   BACKSPACE(fid)
   BACKSPACE(fid)
   read(fid,'(A)') buf
   if(buf(2:10) == 'Alpha orb') exit
  end do ! for while
  read(fid,'(A)') buf
  read(fid,'(A)') buf
  read(buf(20:),*) i, na
  do while(.true.)
   read(fid,'(A)') buf
   if(buf(2:9) == 'Beta orb') exit
  end do ! for while
  read(fid,'(A)') buf
  read(fid,'(A)') buf
  read(buf(20:),*) i, nb
  t1diag = DSQRT(t1diag/DBLE(2*(na+nb)))
 end if

 close(fid)
end subroutine read_cc_e_from_qchem_out

! read the number of doubly occupied orbitals (frozen core orbitals not
!  included) from CFOUR output file
subroutine read_nocc_from_cfour_out(outname, nocc)
 implicit none
 integer :: i, na, nf, fid
 integer, intent(out) :: nocc
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical :: frozen

 nocc = 0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(11:26) == 'Alpha population') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_nocc_from_cfour_out: keyword 'Alph&
                   &a population' not"
  write(6,'(A)') 'found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, ':')
 read(buf(i+1:),*) na ! number of alpha spin orbitals

 frozen = .true.
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(3:15) == 'The following') exit
  if(buf(3:12) == 'Summary of') then
   frozen = .false.
   exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_nocc_from_cfour_out: keyword 'The &
                   &following' not"
  write(6,'(A)') 'found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 if(frozen) then
  read(buf(16:),*) nf ! number of frozen core orbitals
 else
  nf = 0
 end if
 close(fid)

 nocc = na - nf
end subroutine read_nocc_from_cfour_out

! read T1 diagnostic value from a CFOUR T1 binary file
subroutine read_t1diag_from_cfour_t1(nocc, t1diag)
 implicit none
 integer :: i, nvir, fid
 integer(kind=4) :: ibuf(13)
 integer, intent(in) :: nocc
! nocc: the number of doubly occupied orbitals (frozen core not included)
 real(kind=8), intent(out) :: t1diag
 real(kind=8), allocatable :: t1(:)
 logical :: alive

 t1diag = 0d0
 inquire(file='T1',exist=alive)
 if(.not. alive) then
  write(6,'(/,A)') 'ERROR in subroutine read_t1diag_from_cfour_t1: file T1 does&
                   & not exist.'
  stop
 end if

 ! determine nvir from the file size of T1
 call STAT('T1', ibuf)
 nvir = ibuf(8)/(8*nocc)

 allocate(t1(nocc*nvir), source=0d0)
 open(newunit=fid,file='T1',status='old',access='stream')
 read(fid,iostat=i) t1
 close(fid)

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_t1diag_from_cfour_t1: failed to re&
                   &ad t1 amplitudes'
  write(6,'(A)') 'from CFOUR T1 file.'
  write(6,'(A,2I6)') 'nocc, nvir=', nocc, nvir
  deallocate(t1)
  stop
 else
  t1diag = DSQRT(DOT_PRODUCT(t1,t1)/DBLE(nocc*2))
  deallocate(t1)
 end if
end subroutine read_t1diag_from_cfour_t1

subroutine read_cc_e_from_cfour_out(outname, ccd, ccsd, ccsd_t, t1diag, ref_e, &
                                    tot_e)
 implicit none
 integer :: i, k, nocc, fid
 real(kind=8) :: rtmp
 real(kind=8), intent(out) :: t1diag, ref_e, tot_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: ccd, ccsd, ccsd_t
 logical :: uhf

 t1diag = 0d0; ref_e = 0d0; tot_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(13:18) == 'E(SCF)') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_cfour_out: key 'E(SCF)' &
                   &not found in"
  write(6,'(A)') 'file '//TRIM(outname)
  close(fid)
  stop
 end if
 call get_dpv_after_flag(buf, '=', .true., ref_e) ! HF energy

 if(ccd) then
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:11) == 'Total CCD e') then
    call get_dpv_after_flag(buf, ':', .true., tot_e)
    exit
   end if
   if(buf(2:15) == 'A miracle come') then
    BACKSPACE(fid)
    BACKSPACE(fid)
    BACKSPACE(fid)
    read(fid,*) k, rtmp, tot_e
    exit
   end if
  end do ! for while
  close(fid)
  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_cfour_out: failed to lo&
                    &cate CCD energy"
   write(6,'(A)') 'in file '//TRIM(outname)
   stop
  end if
 else if(ccsd) then
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:12)=='Total CCSD e' .or. buf(3:14)=='Total CCSD e') then
    call get_dpv_after_flag(buf, ':', .true., tot_e)
    exit
   end if
  end do ! for while
  close(fid)
  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_cfour_out: 'Total CCSD &
                    &e' is not found"
   write(6,'(A)') "in file "//TRIM(outname)
   stop
  end if
 else if(ccsd_t) then
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:15) == 'Total CCSD(T) e') then
    call get_dpv_after_flag(buf, ':', .true., tot_e)
    exit
   else if(buf(3:13) == 'CCSD(T) ene') then
    call get_dpv_after_flag(buf, 'y', .true., tot_e)
    exit
   end if
  end do ! for while
  close(fid)
  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_cfour_out: neither 'Tot&
                    &al CCSD(T) e' nor"
   write(6,'(A)') "'CCSD(T) ene' is found in file "//TRIM(outname)
   stop
  end if
 end if

 call check_uhf_in_cfour_out(outname, uhf)
 if(uhf .and. (ccsd .or. ccsd_t)) then
  write(6,'(A)') REPEAT('-',79)
  write(6,'(A)') 'T1 diagnostic is currently not supported for UHF-based CC usi&
                 &ng CFOUR. T1_diag'
  write(6,'(A)') 'will be set to 0 below.'
  write(6,'(A)') REPEAT('-',79)
 end if

 if(.not. (uhf .or. ccd)) then
  call read_nocc_from_cfour_out(outname, nocc)
  call read_t1diag_from_cfour_t1(nocc, t1diag)
 end if
end subroutine read_cc_e_from_cfour_out

subroutine read_cc_e_from_molcas_out(outname, ccsd_t, t1diag, ref_e, tot_e)
 implicit none
 integer :: fid
 real(kind=8), intent(out) :: t1diag, ref_e, tot_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical :: ricd, cholesky
 logical, intent(in) :: ccsd_t

 ricd = .false.; cholesky = .false.   ! initialization

 ! determine RICD/CHOLESKY/noRI
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:4) == 'RICD') then
   ricd = .true.
   exit
  end if
  if(buf(1:8) == 'CHOLESKY') then
   cholesky = .true.
   exit
  end if
  if(buf(1:4) == '()()') exit
 end do ! for while
 close(fid)

 if(ricd .or. cholesky) then
  call read_cc_e_from_molcas_out1(outname, ccsd_t, t1diag, ref_e, tot_e)
 else
  call read_cc_e_from_molcas_out2(outname, ccsd_t, t1diag, ref_e, tot_e)
 end if
end subroutine read_cc_e_from_molcas_out

! find the number of occupied orbtials from a (Open)Molcas output file (frozen
! core orbitals not included)
subroutine read_nocc_from_molcas_out(outname, nocc)
 implicit none
 integer :: i, fid
 integer, intent(out) :: nocc
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 nocc = 0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(7:37)=='No. of occupied orbitals with a' .or. buf(2:9)=='Occupied') exit
 end do ! for while

 close(fid)

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_nocc_from_molcas_out: nocc cannot &
                   &be found in file"
  write(6,'(A)') TRIM(outname)
  stop
 end if

 i = INDEX(buf, ':')
 if(i > 0) then
  read(buf(i+1:),*) nocc
 else
  read(buf(47:),*) nocc
 end if
end subroutine read_nocc_from_molcas_out

! read T1 diagnostic value from an (Open)Molcas .RstFil unformatted file
subroutine read_t1diag_from_rstfil(proname, nocc, t1diag)
 implicit none
 integer :: i, nvir, fid
 integer(kind=4) :: ibuf(13)
 integer, intent(in) :: nocc
! nocc: the number of doubly occupied orbitals (frozen core not included)
 real(kind=8), intent(out) :: t1diag
 real(kind=8), allocatable :: t1(:)
 character(len=840) :: fullpath
 character(len=240), intent(in) :: proname

 t1diag = 0d0
 call getenv('MOLCAS_WORKDIR', fullpath)
 fullpath = TRIM(fullpath)//'/'//TRIM(proname)//'/'//TRIM(proname)//'.RstFil'

 ! Determine nvir from the file size of .RstFil. Here 40 means E1 CCSD energy,
 !  E2 CCSD energy and the number of CCSD iterations, so it is excluded.
 call STAT(TRIM(fullpath), ibuf)
 nvir = (ibuf(8)-40)/(8*nocc)

 allocate(t1(nocc*nvir), source=0d0)
 open(newunit=fid,file=TRIM(fullpath),status='old',form='unformatted')
 read(fid,iostat=i) t1
 close(fid)

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_t1diag_from_rstfil: failed to read&
                   & t1 amplitudes'
  write(6,'(A)') 'from Molcas .RstFil file.'
  write(6,'(A)') 'File prefix='//TRIM(proname)
  write(6,'(A,2I6)') 'nocc, nvir=', nocc, nvir
  deallocate(t1)
  stop
 else
  t1diag = DSQRT(DOT_PRODUCT(t1,t1)/DBLE(nocc*2))
  deallocate(t1)
 end if
end subroutine read_t1diag_from_rstfil

! read CCSD/CCSD(T) energies from a (Open)Molcas output file (for RICD)
subroutine read_cc_e_from_molcas_out1(outname, ccsd_t, t1diag, ref_e, tot_e)
 implicit none
 integer :: i, j, nocc, fid
 real(kind=8) :: r(2), corr_e
 real(kind=8), intent(out) :: t1diag, ref_e, tot_e
 character(len=240) :: buf, proname
 character(len=240), intent(in) :: outname
 logical, intent(in) :: ccsd_t

 t1diag = 0d0; ref_e = 0d0; tot_e = 0d0
 call read_nocc_from_molcas_out(outname, nocc)
 call find_specified_suffix(outname, '.out', i)
 proname = outname(1:i-1)
 call read_t1diag_from_rstfil(proname, nocc, t1diag)

 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 if(ccsd_t) then
  do while(.true.)
   BACKSPACE(fid,iostat=i)
   if(i /= 0) exit
   BACKSPACE(fid,iostat=i)
   if(i /= 0) exit
   read(fid,'(A)') buf
   if(buf(4:15) == 'CCSD(T) corr') exit
  end do ! for while

  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_molcas_out1: keyword 'CC&
                    &SD(T) corr' not"
   write(6,'(A)') 'found in file '//TRIM(outname)
   close(fid)
   stop
  end if

  i = INDEX(buf, '=')
  read(buf(i+1:),*) corr_e ! CCSD(T) corr e
 end if

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(3:16) == 'Final CCSD ene') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_molcas_out1: keyword 'Fi&
                   &nal CCSD energy' not"
  write(6,'(A)') 'found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 do i = 1, 2   ! read CCSD correlation e
  read(fid,'(A)') buf
  j = INDEX(buf, ':')
  read(buf(j+1:),*) r(i)
 end do ! for i

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(7:17) == 'Total SCF e') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_molcas_out1: keyword 'To&
                   &tal SCF e' not found"
  write(6,'(A)') 'in file '//TRIM(outname)
  stop
 end if

 read(buf(23:),*) ref_e ! SCF e

 if(ccsd_t) then ! CCSD(T)
  tot_e = ref_e + corr_e
 else            ! CCSD
  tot_e = ref_e + r(1) + r(2)
 end if
end subroutine read_cc_e_from_molcas_out1

! read CCSD/CCSD(T) energies from a (Open)Molcas output file (no RI)
subroutine read_cc_e_from_molcas_out2(outname, ccsd_t, t1diag, ref_e, tot_e)
 implicit none
 integer :: i, nocc, fid
 real(kind=8), intent(out) :: t1diag, ref_e, tot_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: ccsd_t

 t1diag = 0d0; ref_e = 0d0; tot_e = 0d0
 call read_nocc_from_molcas_out(outname, nocc)
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(7:18) == 'Total energy') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_molcas_out2: keyword 'To&
                   &tal energy' not"
  write(6,'(A)') 'found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 ! read CCSD energy
 i = INDEX(buf, ':')
 read(buf(i+1:),*) tot_e

 ! read SCF energy
 read(fid,'(A)') buf
 read(fid,'(A)') buf
 i = INDEX(buf, ':')
 read(buf(i+1:),*) ref_e

 ! read t1 diag
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:15) == 'Euclidian norm') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cc_e_from_molcas_out2: keyword 'Eu&
                   &clidian norm' not"
  write(6,'(A)') 'found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, ':')
 read(buf(i+1:),*) t1diag

 t1diag = DSQRT(t1diag*t1diag/DBLE(2*nocc))

 if(ccsd_t) then ! read CCSD(T) energy
  do while(.true.)
   read(fid,'(A)') buf
   if(buf(8:16) == 'CCSD + T3') exit
  end do ! for while
  i = INDEX(buf, '=')
  read(buf(i+1:),*) tot_e
 end if

 close(fid)
end subroutine read_cc_e_from_molcas_out2

subroutine read_mp2_e_from_mrcc_out(outname, ref_e, mp2_e)
 implicit none
 integer :: fid
 real(kind=8), intent(out) :: ref_e, mp2_e
 character(len=46), parameter :: error_warn = 'ERROR in subroutine read_mp2_e_f&
                                              &rom_mrcc_out: '
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical :: set_ref_e, set_mp2_e

 set_ref_e = .false.; set_mp2_e = .false.; ref_e = 0d0; mp2_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf

  if(buf(2:12)=='Total MP2 e' .or. buf(2:9)=='DF-MP2 e') then
   if(.not. set_mp2_e) then
    call get_dpv_after_flag(buf, ':', .true., mp2_e)
    set_mp2_e = .true.
   end if
  else if(buf(5:24) == 'FINAL HARTREE-FOCK E') then
   if(.not. set_ref_e) then
    call get_dpv_after_flag(buf, ':', .true., ref_e)
    set_ref_e = .true.
    exit
   end if
  else if(buf(2:15) == 'Executing minp') then
   close(fid)
   write(6,'(/,A)') error_warn//'failed to read HF and MP2'
   write(6,'(A)') 'energies from file '//TRIM(outname)
   stop
  end if
 end do ! for while

 close(fid)
 if(.not. (set_ref_e .and. set_mp2_e)) then
  write(6,'(/,A)') error_warn//'some of HF and MP2 energies'
  write(6,'(A)') 'is not found in file '//TRIM(outname)
  stop
 end if
end subroutine read_mp2_e_from_mrcc_out

subroutine read_cc_e_from_mrcc_out(outname, t1diag, ref_e, tot_e)
 implicit none
 integer :: fid
 real(kind=8), intent(out) :: t1diag, ref_e, tot_e
 character(len=45), parameter :: error_warn = 'ERROR in subroutine read_cc_e_fr&
                                              &om_mrcc_out: '
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical :: set_ref_e, set_tot_e

 t1diag = 0d0; ref_e = 0d0; tot_e = 0d0
 set_ref_e = .false.; set_tot_e = .false.
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf

  if(buf(2:14) == 'Total CCSDT e') then
   if(.not. set_tot_e) then
    call get_dpv_after_flag(buf, ':', .true., tot_e)
    set_tot_e = .true.
   end if
  else if(buf(2:16) == 'Total CCSD(T) e') then
   if(.not. set_tot_e) then
    call get_dpv_after_flag(buf, ':', .true., tot_e)
    set_tot_e = .true.
   end if
  else if(buf(2:13) == 'Total CCSD e') then
   if(.not. set_tot_e) then
    call get_dpv_after_flag(buf, ':', .true., tot_e)
    set_tot_e = .true.
   end if
  else if(buf(2:8) == 'T1 diag') then
   call get_dpv_after_flag(buf, ':', .true., t1diag)
  else if(buf(5:24) == 'FINAL HARTREE-FOCK E') then
   if(.not. set_ref_e) then
    call get_dpv_after_flag(buf, ':', .true., ref_e)
    set_ref_e = .true.
    exit
   end if
  else if(buf(2:15) == 'Executing minp') then
   close(fid)
   write(6,'(/,A)') error_warn//'failed to read HF and CC'
   write(6,'(A)') 'energies from file '//TRIM(outname)
   stop
  end if
 end do ! for while

 close(fid)
 if(.not. (set_ref_e .and. set_tot_e)) then
  write(6,'(/,A)') error_warn//'some of HF and CC energies'
  write(6,'(A)') 'is not found in file '//TRIM(outname)
  stop
 end if
end subroutine read_cc_e_from_mrcc_out

