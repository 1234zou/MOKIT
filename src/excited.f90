! written by jxzou at 20210128: subroutines involving excited states wave function
! Currently only CIS/TDHF/TDDFT can be analyzed

! read NROrb, NOA, NOB from a Gaussian output file
! NROrb(nif_ex): number of spatial orbitals involved in excitation, nif_ex<=nif
! NOA     (noa): number of alpha occupied spin orbitals
! NOB     (nob): number of beta occupied spin orbitals
subroutine read_noa_and_nob_from_gau_log(logname, nif_ex, noa, nob)
 implicit none
 integer :: i, fid
 integer, intent(out) :: nif_ex, noa, nob
 character(len=240) :: buf
 character(len=240), intent(in) :: logname

 open(newunit=fid,file=TRIM(logname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:6) == 'Range') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_noa_and_nob_from_gau_log: no 'Range'&
                & found in file "//TRIM(logname)
  close(fid)
  stop
 end if

 read(fid,'(A)') buf
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
end subroutine read_noa_and_nob_from_gau_log

! read excitation coefficients from Gaussian output file
subroutine read_ex_coeff_from_gau_log(logname, istate, nocc, nvir, exc)
 implicit none
 integer :: i, fid
 integer, intent(in) :: istate, nocc, nvir
 real(kind=8), intent(out) :: exc(nocc,nvir)
 character(len=240) :: buf
 character(len=240), intent(in) :: logname

 open(newunit=fid,file=TRIM(logname),status='old',position='rewind')
 write(str,'(A13,I4)') 'Excited State', istate
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:18) == str) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine : no '"//str//"' found in file "//&
                   TRIM(logname)
  close(fid)
  stop
 end if

 allocate(exc(nocc,nvir), source=0d0)

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:11)=='This state' .or. LEN_TRIM(buf)==0) exit
  i = index(buf,'->')
  if(i == 0) i = index(buf,'<-') ! de-excitation
  read(buf(1:i-1),*) j
  read(buf(i+2:),*) k
  k = k - nae   ! remember to minus nae
  read(buf(23:),*) rtmp
  exc(j,k) = exc(j,k) + rtmp
  ! for de-excitations, simply sum, learned from
  !  http://sobereva.com/wfnbbs/viewtopic.php?id=72
 end do ! for while

 close(fid)
end subroutine read_ex_coeff_from_gau_log

! generate AO-basis Transition Density Matrix from CI coefficients in Gaussian
! .log file
! ground state -> i-th excited state (e.g. istate=1 for S1 state)
subroutine gen_ao_tdm(logname, nbf, nif, mo, istate, tdm)
 implicit none
 integer :: i, j, k, nae, nbe, nav, nbv, fid
 integer, parameter :: iout = 6
 integer, intent(in) :: nbf, nif, istate
 real(kind=8) :: rtmp
 real(kind=8), intent(in) :: mo(nbf,nif)
 real(kind=8), intent(out) :: tdm(nbf,nbf)
 real(kind=8), allocatable :: exc(:,:)
 character(len=17) :: str
 character(len=240) :: buf
 character(len=240), intent(in) :: logname

 tdm = 0d0
 call mo2ao_tdm(nbf, nif, mo, nae, nav, exc, .true., tdm)
 deallocate(exc)
end subroutine gen_ao_tdm

subroutine mo2ao_tdm(nbf, nif, mo, nae, nav, exc, total, tdm)
 implicit none
 integer, intent(in) :: nbf, nif, nae, nav
 real(kind=8) :: rtmp
 real(kind=8), intent(in) :: mo(nbf,nif), exc(nae,nav)
 real(kind=8), intent(out) :: tdm(nbf,nbf)
 real(kind=8), allocatable :: mo_ci(:,:)
 logical, intent(in) :: total

 ! call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
 allocate(mo_ci(nbf,nav), source=0d0)
 call dgemm('N','N',nbf,nav,nae,1d0,mo(:,1:nae),nbf,exc,nae,0d0,mo_ci,nbf)

 ! alpha + beta, or only one of alpha/beta
 rtmp = 1d0
 if(total) rtmp = 2d0
 call dgemm('N','T',nbf,nbf,nav,rtmp,mo_ci,nbf,mo(:,nae+1:nif),nbf,0d0,tdm,nbf)
 deallocate(mo_ci)
 return
end subroutine mo2ao_tdm

! Compute and average the first n states of AO-based transition density matrices
subroutine average_cis_tdm(fchname, logname, n, nbf, nif, tdm)
 implicit none
 integer :: i
 integer, intent(in) :: n, nbf, nif
 real(kind=8), intent(out) :: tdm(nbf,nbf)
 real(kind=8), allocatable :: alpha_mo(:,:), beta_mo(:,:), tdm1(:,:)
 character(len=240), intent(in) :: fchname, logname
 logical :: uhf

 allocate(alpha_mo(nbf,nif))
 call read_mo_from_fch(fchname, nbf, nif, 'a', alpha_mo)

 call check_uhf_in_fch(fchname, uhf)
 if(uhf) then
  allocate(beta_mo(nbf,nif))
  call read_mo_from_fch(fchname, nbf, nif, 'b', beta_mo)
 end if

 tdm = 0d0
 allocate(tdm1(nbf,nbf), source=0d0)
 do i = 1, n, 1
  call gen_ao_tdm(logname, nbf, nif, alpha_mo, i, tdm1)
  tdm = tdm + tdm1
 end do ! for i
 deallocate(alpha_mo)

 if(uhf) then
  do i = 1, n, 1
   call gen_ao_tdm(logname, nbf, nif, beta_mo, i, tdm1)
   tdm = tdm + tdm1
  end do ! for i
  deallocate(beta_mo)
 end if
 deallocate(tdm1)

 tdm = tdm/DBLE(n) ! averaged transition density
end subroutine average_cis_tdm

! get State-Averaged NTOs from specified .fch and .log files
subroutine get_state_averaged_nto(fchname, logname, nstate)
 use util_wrapper, only: fch_u2r_wrap
 implicit none
 integer :: i, nbf, nif
 integer, intent(in) :: nstate
 character(len=240) :: fchname1, fchname2
 character(len=240), intent(in) :: fchname, logname
 real(kind=8), allocatable :: mo(:,:), tdm(:,:), S(:,:), noon(:)
 logical :: uhf

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(tdm(nbf,nbf))
 call average_cis_tdm(fchname, logname, nstate, nbf, nif, tdm)

 allocate(S(nbf,nbf))
 call get_ao_ovlp_using_fch(fchname, nbf, S)

 allocate(noon(nif), mo(nbf,nif))
 call no(nbf, nif, tdm, S, noon, mo)
 deallocate(tdm, S)

 i = index(fchname, '.fch')
 fchname1 = fchname(1:i-1)//'_SA-NTO.fch'
 fchname2 = fchname(1:i-1)//'_r.fch'
 call check_uhf_in_fch(fchname, uhf)
 if(uhf) then
  call fch_u2r_wrap(fchname)
  call copy_file(fchname2, fchname1, .true.)
 else
  call copy_file(fchname, fchname1, .false.)
 end if

 call write_mo_into_fch(fchname1, nbf, nif, 'a', mo)
 call write_eigenvalues_to_fch(fchname1, nif, 'a', noon, .true.)
 deallocate(noon, mo)
end subroutine get_state_averaged_nto

! check whether UHF-type MOs are hold in a given .fch(k) file
subroutine check_uhf_in_fch(fchname, uhf)
 implicit none
 integer :: i, fid
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
 logical, intent(out) :: uhf

 uhf = .false.
 call open_file(fchname, .true., fid)

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i < 0) exit ! end-of-file
  if(i > 0) then
   write(6,'(A)') 'ERROR in subroutine check_uhf_in_fch: failed to read file '&
                   //TRIM(fchname)
   stop
  end if

  if(buf(1:7) == 'Beta MO') then
   uhf = .true.
   exit
  end if
 end do ! for while

 close(fid)
end subroutine check_uhf_in_fch

