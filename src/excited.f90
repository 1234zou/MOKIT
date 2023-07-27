! written by jxzou at 20210128: subroutines involving excited states wave function
! Currently only CIS/TDHF/TDDFT can be analyzed

! read NROrb, NOA, NOB from a Gaussian output file
! NFC          : number of frozen doubly occupied orbitals
! NROrb(nif_ex): number of spatial orbitals involved in excitation, nif_ex<=nif
! NOA     (noa): number of alpha occupied spin orbitals involved in excitation
! NOB     (nob): number of beta occupied spin orbitals involved in excitation
subroutine read_noa_nob_from_gau_log(logname, nfc, nif_ex, noa, nob)
 implicit none
 integer :: i, fid
 integer, intent(out) :: nfc, nif_ex, noa, nob
 character(len=240) :: buf
 character(len=240), intent(in) :: logname

 open(newunit=fid,file=TRIM(logname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:9) == 'Range of') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_noa_nob_from_gau_log: no 'Range'&
                & found in file "//TRIM(logname)
  close(fid)
  stop
 end if

 read(fid,'(A)') buf
 i = index(buf, 'NFC=')
 read(buf(i+4:),*) nfc

 read(fid,'(A)') buf
 close(fid)
 i = index(buf, '=')
 read(buf(i+1:),*) nif_ex
 buf(i:i) = ' '

 i = index(buf, '=')
 read(buf(i+1:),*) noa
 buf(i:i) = ' '

 i = index(buf, '=')
 read(buf(i+1:),*) nob
end subroutine read_noa_nob_from_gau_log

! Check whether noa and nob is consistent with input nocc, nvir.
! If consistent, return nfc, nif_ex; if not, stop
subroutine check_noa_nob_in_logname(logname, alpha, nocc, nvir, nfc, nif_ex)
 implicit none
 integer :: noa, nob, nva, nvb
 integer, intent(in) :: nocc, nvir
 integer, intent(out) :: nfc, nif_ex
 logical, intent(in) :: alpha
 character(len=240), intent(in) :: logname

 call read_noa_nob_from_gau_log(logname, nfc, nif_ex, noa, nob)
 nva = nif_ex - noa
 nvb = nif_ex - nob

 if(alpha) then ! check alpha excitation information
  if(nocc/=noa .or. nvir/=nva) then
   write(6,'(A)') 'ERROR in subroutine check_noa_nob_in_logname: nocc/=noa or&
                 & nvir/=nva.'
   write(6,'(A)') 'logname='//TRIM(logname)
   write(6,'(A,4I5)') 'nocc, noa,  nvir, nva=', nocc, noa,  nvir, nva
   stop
  end if
 else ! check beta excitation information
  if(nocc/=nob .or. nvir/=nvb) then
   write(6,'(A)') 'ERROR in subroutine check_noa_nob_in_logname: nocc/=nob or&
                 & nvir/=nvb.'
   write(6,'(A)') 'logname='//TRIM(logname)
   write(6,'(A,4I5)') 'nocc, nob,  nvir, nvb=', nocc, nob,  nvir, nvb
   stop
  end if
 end if

end subroutine check_noa_nob_in_logname

! read excitation coefficients from Gaussian output file
! This subroutine is only used for RHF-based CIS/TDHF/TDDFT
! For UHF-based CIS/TDHF/TDDFT, please use read_ex_coeff_from_gau_log2
subroutine read_ex_coeff_from_gau_log(logname, istate, nocc, nvir, exc)
 implicit none
 integer :: i, j, k, fid, nfc, nif_ex
 integer, intent(in) :: istate, nocc, nvir
! nocc: number of alpha occupied orbitals involved in excitation
! nvir: number of beta occupied orbitals involved in excitation
 real(kind=8) :: rtmp
 real(kind=8), intent(out) :: exc(nocc,nvir)
 character(len=17) :: str
 character(len=240) :: buf
 character(len=240), intent(in) :: logname

 call check_noa_nob_in_logname(logname, .true., nocc, nvir, nfc, nif_ex)

 open(newunit=fid,file=TRIM(logname),status='old',position='rewind')
 write(str,'(A13,I4)') 'Excited State', istate
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:18) == str) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_ex_coeff_from_gau_log: no '"//&
                   str//"' found in file "//TRIM(logname)
  close(fid)
  stop
 end if
 exc = 0d0

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:11)=='This state' .or. buf(2:7)=='SavETr' .or. LEN_TRIM(buf)==0) exit
  i = index(buf,'->')
  if(i == 0) i = index(buf,'<-') ! de-excitation
  if(i == 0) then
   close(fid)
   write(6,'(A)') "ERROR in subroutine read_ex_coeff_from_gau_log: no '->'&
                 & or '<-' symbol found"
   write(6,'(A)') 'in file '//TRIM(logname)
   write(6,'(A)') TRIM(buf)
   stop
  end if
  read(buf(1:i-1),*) j
  read(buf(i+2:),*) k
  j = j - nfc        ! remember to minus frozen core
  k = k - nfc - nocc ! minus number of alpha occupied orbitals
  read(buf(22:),*) rtmp
  exc(j,k) = exc(j,k) + rtmp
  ! for de-excitations, simply sum, learned from
  !  http://sobereva.com/wfnbbs/viewtopic.php?id=72
 end do ! for while

 close(fid)
end subroutine read_ex_coeff_from_gau_log

! read excitation coefficients from Gaussian output file
! This subroutine is only used for UHF-based CIS/TDHF/TDDFT
! For RHF-based CIS/TDHF/TDDFT, please use read_ex_coeff_from_gau_log
subroutine read_ex_coeff_from_gau_log2(logname, istate, nocc_a, nvir_a, nocc_b,&
                                       nvir_b, exc_a, exc_b)
 implicit none
 integer :: i, j, k, fid, nfc, nif_ex
 integer, intent(in) :: istate, nocc_a, nvir_a, nocc_b, nvir_b
! nocc: number of alpha occupied orbitals involved in excitation
! nvir: number of beta occupied orbitals involved in excitation
 real(kind=8) :: rtmp
 real(kind=8), intent(out) :: exc_a(nocc_a,nvir_a), exc_b(nocc_b,nvir_b)
 character(len=6) :: str6
 character(len=17) :: str
 character(len=240) :: buf
 character(len=240), intent(in) :: logname
 logical :: alpha

 call check_noa_nob_in_logname(logname, .true., nocc_a, nvir_a, nfc, nif_ex)
 call check_noa_nob_in_logname(logname,.false., nocc_b, nvir_b, nfc, nif_ex)

 open(newunit=fid,file=TRIM(logname),status='old',position='rewind')
 write(str,'(A13,I4)') 'Excited State', istate
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:18) == str) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_ex_coeff_from_gau_log2: no'"//&
                   str//"' found in file "//TRIM(logname)
  close(fid)
  stop
 end if
 exc_a = 0d0; exc_b = 0d0

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:11)=='This state' .or. LEN_TRIM(buf)==0) exit
  i = index(buf,'->')
  if(i == 0) i = index(buf,'<-') ! de-excitation
  if(i == 0) then
   close(fid)
   write(6,'(A)') "ERROR in subroutine read_ex_coeff_from_gau_log2: no '->'&
                 & or '<-' symbol found"
   write(6,'(A)') 'in file '//TRIM(logname)
   stop
  end if

  read(buf(1:i-1),*) str6
  call del_A_or_B(str6)
  read(str6,*) j
  if(index(str,'A') > 0) then
   alpha = .true.
  else
   alpha = .false.
  end if

  read(buf(i+2:),*) str6
  call del_A_or_B(str6)
  read(str6,*) k
  j = j - nfc        ! remember to minus frozen core
  if(alpha) then
   k = k - nfc - nocc_a ! minus number of alpha occupied orbitals
  else
   k = k - nfc - nocc_b ! minus number of beta occupied orbitals
  end if
  read(buf(22:),*) rtmp
  if(alpha) then
   exc_a(j,k) = exc_a(j,k) + rtmp
  else
   exc_b(j,k) = exc_b(j,k) + rtmp
  end if
  ! for de-excitations, simply sum, learned from
  !  http://sobereva.com/wfnbbs/viewtopic.php?id=72
 end do ! for while

 close(fid)
end subroutine read_ex_coeff_from_gau_log2

! generate NTO for restricted-type CIS/TDHF/TDDFT wave function
! T*(T^) = U(n_occ)U^ -> occupied NTOs
! (T^)*T = V(n_vir)V^ -> unoccupied NTOs
subroutine gen_nto(logname, fchname, istate, averaged)
 implicit none
 integer :: i, nfc, nif_ex, nocc, nob, nvir, ndb, nbf, nif
 integer, intent(in) :: istate
 real(kind=8), allocatable :: exc0(:,:)
 real(kind=8), allocatable :: exc(:,:), ttd(:,:), ev(:), ev2(:)
 real(kind=8), allocatable :: mo(:,:), occ_mo(:,:), vir_mo(:,:)
 character(len=240), intent(in) :: logname, fchname
 ! logname: it includes CIS excitation coefficients in .log file
 ! fchname: it includes canonical MOs in .fch file
 logical, intent(in) :: averaged
 ! True: State-Averaged NTO; False: NTO of the ground state -> a excited state

 call read_noa_nob_from_gau_log(logname, nfc, nif_ex, nocc, nob)
 ndb = nfc + nocc
 nvir = nif_ex - nocc
 allocate(exc(nocc,nvir), source=0d0)

 if(averaged) then
  write(6,'(A,I0,A)') 'State-Averaged NTO generated using ',istate,' states.'
  allocate(exc0(nocc,nvir))
  do i = 1, istate, 1
   call read_ex_coeff_from_gau_log(logname, i, nocc, nvir, exc0)
   exc = exc + exc0
  end do ! for i
  deallocate(exc0)
  exc = DSQRT(2d0)*exc/DSQRT(DBLE(istate)) ! sum of square of exc is 0.5

 else ! ground state -> a excited state
  call read_ex_coeff_from_gau_log(logname, istate, nocc, nvir, exc)
  exc = DSQRT(2d0)*exc  ! sum of square of exc is 0.5, so multiply root2
 end if

 allocate(ttd(nocc,nocc), source=0d0) ! T*(T^dagger)
 ! call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
 call dgemm('N','T',nocc,nocc,nvir,1d0,exc,nocc,exc,nocc,0d0,ttd,nocc)

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(ev(nif), source=0d0)
 call diag_get_e_and_vec(nocc, ttd, ev(nfc+1:ndb)) ! occupied NTO eigenvalues

 allocate(mo(nbf,nif), occ_mo(nbf,nocc))
 call read_mo_from_fch(fchname, nbf, nif, 'a', mo)
 occ_mo = mo(:,nfc+1:ndb) ! C'_occ = (C_occ)U
 call dgemm('N','N',nbf,nocc,nocc,1d0,occ_mo,nbf,ttd,nocc,0d0,mo(:,nfc+1:ndb),nbf)
 deallocate(occ_mo, ttd)

 allocate(ttd(nvir,nvir), source=0d0) ! (T^dagger)*T
 call dgemm('T','N',nvir,nvir,nocc,1d0,exc,nocc,exc,nocc,0d0,ttd,nvir)

 ! unoccupied NTO eigenvalues
 call diag_get_e_and_vec(nvir, ttd, ev(ndb+1:ndb+nvir))
 deallocate(exc)
 allocate(vir_mo(nbf,nvir))
 vir_mo = mo(:,ndb+1:ndb+nvir) ! C'_vir = (C_vir)U
 call dgemm('N','N',nbf,nvir,nvir,1d0,vir_mo,nbf,ttd,nvir,0d0,mo(:,ndb+1:ndb+nvir),nbf)
 deallocate(ttd)

 ! reverse the eigenvalues and virtual NTOs
 allocate(ev2(nvir))
 do i = ndb+1, ndb+nvir, 1
  ev2(i-ndb) = ev(2*ndb+nvir-i+1)
  vir_mo(:,i-ndb) = mo(:,2*ndb+nvir-i+1)
 end do ! for i
 ev(ndb+1:ndb+nvir) = ev2
 mo(:,ndb+1:ndb+nvir) = vir_mo
 deallocate(ev2, vir_mo)

 call write_mo_into_fch(fchname, nbf, nif, 'a', mo)
 call write_eigenvalues_to_fch(fchname, nif, 'a', ev, .true.)
end subroutine gen_nto

! generate AO-basis Transition Density Matrix from CI coefficients in Gaussian
! .log file
! ground state -> i-th excited state (e.g. istate=1 for S1 state)
subroutine gen_ao_tdm(logname, nbf, nif, mo, istate, tdm)
 implicit none
 integer :: nfc, nif_ex, nob, nocc, nvir
 integer, intent(in) :: nbf, nif, istate
 real(kind=8), intent(in) :: mo(nbf,nif)
 real(kind=8), intent(out) :: tdm(nbf,nbf)
 real(kind=8), allocatable :: exc(:,:)
 character(len=240), intent(in) :: logname

 call read_noa_nob_from_gau_log(logname, nfc, nif_ex, nocc, nob)
 nvir = nif_ex - nocc
 allocate(exc(nocc,nvir))
 call read_ex_coeff_from_gau_log(logname, istate, nocc, nvir, exc)
 call mo2ao_tdm(nbf, nif, mo, nocc, nvir, exc, .true., tdm)
 deallocate(exc)
end subroutine gen_ao_tdm

! generate AO-basis Transition Density Matrix from CI coefficients in Gaussian
! .log file
! ground state -> i-th excited state (e.g. istate=1 for S1 state)
subroutine gen_ao_tdm2(logname, nbf, nif, mo, istate, tdm)
 implicit none
 integer :: nfc, nif_ex, nocc_a, nocc_b, nvir_a, nvir_b
 integer, intent(in) :: nbf, nif, istate
 real(kind=8), intent(in) :: mo(nbf,nif,2)
 real(kind=8), intent(out) :: tdm(nbf,nbf,2)
 real(kind=8), allocatable :: exc_a(:,:), exc_b(:,:)
 character(len=240), intent(in) :: logname

 call read_noa_nob_from_gau_log(logname, nfc, nif_ex, nocc_a, nocc_b)
 nvir_a = nif_ex - nocc_a
 nvir_b = nif_ex - nocc_b
 allocate(exc_a(nocc_a,nvir_a), exc_b(nocc_b,nvir_b))
 call read_ex_coeff_from_gau_log2(logname, istate, nocc_a, nvir_a, nocc_b, &
                                  nvir_b, exc_a, exc_b)
 call mo2ao_tdm(nbf,nif, mo(:,:,1), nocc_a,nvir_a, exc_a, .false., tdm(:,:,1))
 call mo2ao_tdm(nbf,nif, mo(:,:,2), nocc_b,nvir_b, exc_b, .false., tdm(:,:,2))
 deallocate(exc_a, exc_b)
end subroutine gen_ao_tdm2

subroutine mo2ao_tdm(nbf, nif, mo, nae, nav, exc, total, tdm)
 implicit none
 integer, intent(in) :: nbf, nif, nae, nav
 real(kind=8) :: rtmp
 real(kind=8), intent(in) :: mo(nbf,nif), exc(nae,nav)
 real(kind=8), intent(out) :: tdm(nbf,nbf)
 real(kind=8), allocatable :: mo_ci(:,:)
 logical, intent(in) :: total

 tdm = 0d0
 allocate(mo_ci(nbf,nav), source=0d0)
 call dgemm('N','N',nbf,nav,nae,1d0,mo(:,1:nae),nbf,exc,nae,0d0,mo_ci,nbf)

 ! alpha + beta, or only one of alpha/beta
 rtmp = 1d0
 if(total) rtmp = 2d0
 call dgemm('N','T',nbf,nbf,nav,rtmp,mo_ci,nbf,mo(:,nae+1:nif),nbf,0d0,tdm,nbf)
 deallocate(mo_ci)
end subroutine mo2ao_tdm

! Compute and average the first n states of AO-based transition density matrices
subroutine average_cis_tdm(fchname, logname, n, nbf, nif, tdm)
 implicit none
 integer :: i
 integer, intent(in) :: n, nbf, nif
 real(kind=8), intent(out) :: tdm(nbf,nbf)
 real(kind=8), allocatable :: mo1(:,:), tdm1(:,:)
 real(kind=8), allocatable :: mo2(:,:,:), tdm2(:,:,:)
 character(len=240), intent(in) :: fchname, logname
 logical :: uhf

 tdm = 0d0
 call check_uhf_in_fch(fchname, uhf)

 if(uhf) then
  allocate(mo2(nbf,nif,2))
  call read_mo_from_fch(fchname, nbf, nif, 'a', mo2(:,:,1))
  call read_mo_from_fch(fchname, nbf, nif, 'b', mo2(:,:,2))
  allocate(tdm2(nbf,nbf,2), source=0d0)
  do i = 1, n, 1
   call gen_ao_tdm2(logname, nbf, nif, mo2, i, tdm2)
   tdm = tdm + tdm2(:,:,1) + tdm2(:,:,2)
  end do ! for i
  deallocate(mo2, tdm2)

 else
  allocate(mo1(nbf,nif))
  call read_mo_from_fch(fchname, nbf, nif, 'a', mo1)
  allocate(tdm1(nbf,nbf), source=0d0)
  do i = 1, n, 1
   call gen_ao_tdm(logname, nbf, nif, mo1, i, tdm1)
   tdm = tdm + tdm1
  end do ! for i
  deallocate(mo1, tdm1)
 end if

 tdm = tdm/DBLE(n) ! averaged transition density
end subroutine average_cis_tdm

subroutine del_A_or_B(str)
 implicit none
 integer :: i
 character(len=6), intent(inout) :: str

 str = ADJUSTL(str)
 i = LEN_TRIM(str)
 if(str(i:i)=='A' .or. str(i:i)=='B') str(i:i) = ' '
end subroutine del_A_or_B

! check whether UHF-type MOs are hold in a given .fch(k) file
subroutine check_uhf_in_fch(fchname, uhf)
 implicit none
 integer :: i, fid
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname
 logical, intent(out) :: uhf
!f2py intent(out) :: uhf

 uhf = .false.
 call open_file(fchname, .true., fid)

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i < 0) exit ! end-of-file
  if(i > 0) then
   write(6,'(A)') 'ERROR in subroutine check_uhf_in_fch: failed to read file&
                  & '//TRIM(fchname)
   close(fid)
   stop
  end if

  if(buf(1:7) == 'Beta MO') then
   uhf = .true.
   exit
  end if

  select case(buf(1:11))
  case('Orthonormal','Total SCF D','Mulliken Ch')
   exit
  end select
 end do ! for while

 close(fid)
end subroutine check_uhf_in_fch

! calculate the AO-based transition density matrix and oscillator strength,
! generate particle and hole NTOs, respectively
subroutine gen_nto_and_fosc_from_mo_tdm(nbf, nmo, mo, mo_tdm, ao_dip, delta_e, &
                                        ev, part_mo, hole_mo, fosc)
 implicit none
 integer :: i
 integer, intent(in) :: nbf, nmo
!f2py intent(in) :: nbf, nmo
 real(kind=8) :: td(3) ! transition dipole moment
 real(kind=8), intent(in) :: mo(nbf,nmo), mo_tdm(nmo,nmo), ao_dip(3,nbf,nbf), &
  delta_e
!f2py intent(in) :: mo, mo_tdm, ao_dip, delta_e
!f2py depend(nbf,nmo) :: mo
!f2py depend(nbf) :: ao_dip
!f2py depend(nmo) :: mo_tdm
 real(kind=8), intent(out) :: ev(nmo), part_mo(nbf,nmo), hole_mo(nbf,nmo), fosc
!f2py intent(out) :: ev, part_mo, hole_mo, fosc
!f2py depend(nbf,nmo) :: part_mo, hole_mo
!f2py depend(nmo) :: ev
 real(kind=8), allocatable :: ao_tdm(:,:), ttd(:,:)

 ev = 0d0; part_mo = 0d0; hole_mo = 0d0; fosc = 0d0
 allocate(ao_tdm(nbf,nbf))
 call calc_CXCT(nbf, nmo, mo, mo_tdm, ao_tdm)

 do i = 1, 3
  call get_ne_from_PS(nbf, ao_tdm, ao_dip(i,:,:), td(i))
 end do ! for i

 deallocate(ao_tdm)
 fosc = 2d0*DOT_PRODUCT(td, td)*delta_e/3d0

 ! T*(T^dagger)
 allocate(ttd(nmo,nmo), source=0d0)
 call dgemm('N','T',nmo,nmo,nmo, 1d0, mo_tdm, nmo, mo_tdm, nmo, 0d0, ttd, nmo)
 call diag_get_e_and_vec(nmo, ttd, ev)
 ev = ev/SUM(ev)
 ! compute particle NTO, C_NTO = CU
 call dgemm('N','N',nbf,nmo,nmo, 1d0, mo, nbf, ttd, nmo, 0d0, part_mo, nbf)

 ! (T^dagger)*T
 ttd = 0d0
 call dgemm('T','N',nmo,nmo,nmo, 1d0, mo_tdm, nmo, mo_tdm, nmo, 0d0, ttd, nmo)
 call diag_get_e_and_vec(nmo, ttd, ev)
 ev = ev/SUM(ev)
 ! compute hole NTO, C_NTO = CU
 call dgemm('N','N',nbf,nmo,nmo, 1d0, mo, nbf, ttd, nmo, 0d0, hole_mo, nbf)

 deallocate(ttd)
end subroutine gen_nto_and_fosc_from_mo_tdm

