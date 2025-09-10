! written by jxzou at 20210128: subroutines involving excited states wave function
! TODO: implement NTO using non-orthogonal orbitals, e.g.
!       singlet RHF/RKS -> triplet ROHF/ROKS/UHF/UKS
! TODO: implement Dyson orbitals using non-orthogonal orbitals, e.g.
!       N-e singlet RHF/RKS -> 'N-1'-e doublet ROHF/ROKS/UHF/UKS

subroutine calc_fosc_using_mo_tdm_ao_dip(nbf, nmo, delta_e, mo, mo_tdm, ao_dip,&
                                         fosc)
 implicit none
 integer, intent(in) :: nbf, nmo
!f2py intent(in) :: nbf, nmo
 real(kind=8) :: td(3) ! transition dipole moment
 real(kind=8), intent(in) :: delta_e, mo(nbf,nmo), mo_tdm(nmo,nmo), ao_dip(3,nbf,nbf)
!f2py intent(in) :: delta_e, mo, mo_tdm, ao_dip
!f2py depend(nbf,nmo) :: mo
!f2py depend(nmo) :: mo_tdm
!f2py depend(nbf) :: ao_dip
 real(kind=8), intent(out) :: fosc
!f2py intent(out) :: fosc
 real(kind=8), allocatable :: ao_tdm(:,:)

 allocate(ao_tdm(nbf,nbf))
 call calc_CXCT(nbf, nmo, mo, mo_tdm, ao_tdm)
 call trace_dm_ao_ovlp(nbf, ao_tdm, ao_dip(1,:,:), td(1))
 call trace_dm_ao_ovlp(nbf, ao_tdm, ao_dip(2,:,:), td(2))
 call trace_dm_ao_ovlp(nbf, ao_tdm, ao_dip(3,:,:), td(3))
 deallocate(ao_tdm)
 fosc = 2d0*DOT_PRODUCT(td, td)*delta_e/3d0
end subroutine calc_fosc_using_mo_tdm_ao_dip

! calculate RHF-CIS MO-based density matrix using excitation coefficients
! nfc: the number of frozen core orbitals in CIS calculation
! nocc: the number of doubly occupied orbitals involved in excitations
! nvir: the number of virtual orbitals involved in excitations
! Note: |C_ia|^2 = 1 is required for the input array exc.
subroutine calc_cis_mo_dm_using_exc(nfc, nocc, nvir, exc, dm)
 implicit none
 integer :: i, j, p, q, ndb, nmo
 integer, intent(in) :: nfc, nocc, nvir
 real(kind=8) :: r_ne
 real(kind=8), parameter :: ne_thres = 1d-2
 real(kind=8), intent(in) :: exc(nocc,nvir)
 real(kind=8), intent(out) :: dm(nfc+nocc+nvir,nfc+nocc+nvir)

 dm = 0d0; ndb = nfc + nocc; nmo = ndb + nvir
 if(nfc > 0) then
  forall(i = 1:nfc) dm(i,i) = 2d0
 end if

 do i = 1, nocc, 1 ! d_ii
  j = nfc + i
  dm(j,j) = 2d0 - DOT_PRODUCT(exc(i,:),exc(i,:))
 end do ! for i

!$omp parallel do schedule(dynamic) default(shared) private(i,j)
 do i = 1, nvir, 1 ! d_aa
  j = ndb + i
  dm(j,j) = DOT_PRODUCT(exc(:,i),exc(:,i))
 end do ! for i
!$omp end parallel do

 r_ne = 0d0
 do i = 1, nmo, 1
  r_ne = r_ne + dm(i,i)
 end do ! for i

 if(DABS(r_ne - 2d0*DBLE(ndb)) > ne_thres) then
  write(6,'(/,A)') 'ERROR in subroutine calc_cis_mo_dm_using_exc: the number of&
                   & electrons calc-'
  write(6,'(A)') 'ulated from exc differs from 2*ndb.'
  write(6,'(A,F12.4,A,I0)') 'r_ne=', r_ne, ', ndb=', ndb
  stop
 end if

 do i = 1, nocc-1, 1 ! d_ij
  p = nfc + i
  do j = i+1, nocc, 1
   q = nfc + j
   dm(p,q) = -DOT_PRODUCT(exc(i,:),exc(j,:))
  end do ! for j
 end do ! for i

!$omp parallel do schedule(dynamic) default(shared) private(i,j,p,q)
 do i = 1, nvir-1, 1 ! d_ab
  p = ndb + i
  do j = i+1, nvir, 1
   q = ndb + j
   dm(p,q) = DOT_PRODUCT(exc(:,i),exc(:,j))
  end do ! for j
 end do ! for i
!$omp end parallel do

 ! RHF-CIS d_ia = 0, no need to calculate.
 call symmetrize_dmat(nmo, dm)
end subroutine calc_cis_mo_dm_using_exc

! Calculte RHF-CIS MO-based transition density matrix between two CIS excited
! states with the same spin (both singlet or both triplet).
! Note: the RHF->CIS (ground state -> singlet excited state) transition density
!       matrix is simply DSQRT(2d0)*C_ia, thus no subroutine is needed.
subroutine calc_cis_mo_tdm_using_exc(ndb, nvir, exc, tdm_occ, tdm_vir)
 implicit none
 integer :: i, j
 integer, intent(in) :: ndb, nvir
!f2py intent(in) :: ndb, nvir
 real(kind=8), intent(in) :: exc(ndb,nvir,2)
!f2py intent(in) :: exc
!f2py depend(ndb,nvir) :: exc
 real(kind=8), intent(out) :: tdm_occ(ndb,ndb), tdm_vir(nvir,nvir)
!f2py intent(out) :: tdm_occ, tdm_vir
!f2py depend(ndb) :: tdm_occ
!f2py depend(nvir) :: tdm_vir

 tdm_occ = 0d0; tdm_vir = 0d0

 do i = 1, ndb, 1 ! d_ii, i <- {C}
  tdm_occ(i,i) = -DOT_PRODUCT(exc(i,:,1), exc(i,:,2))
 end do ! for i

 do i = 1, ndb-1, 1 ! d_ij, ij <- {C}
  do j = i+1, ndb, 1
   tdm_occ(i,j) = -DOT_PRODUCT(exc(i,:,2), exc(j,:,1))
   tdm_occ(j,i) = -DOT_PRODUCT(exc(j,:,2), exc(i,:,1))
  end do ! for j
 end do ! for i

!$omp parallel do schedule(dynamic) default(shared) private(i,j)
 do i = 1, nvir, 1 ! d_aa, a <- {V}
  j = ndb + i
  tdm_vir(j,j) = DOT_PRODUCT(exc(:,i,1), exc(:,i,2))
 end do ! for i
!$omp end parallel do

!$omp parallel do schedule(dynamic) default(shared) private(i,j)
 do i = 1, nvir-1, 1 ! d_ab, ab <- {V}
  do j = i+1, nvir, 1
   tdm_vir(i,j) = DOT_PRODUCT(exc(:,i,1), exc(:,j,2))
   tdm_vir(j,i) = DOT_PRODUCT(exc(:,j,1), exc(:,i,2))
  end do ! for j
 end do ! for i
!$omp end parallel do
end subroutine calc_cis_mo_tdm_using_exc

! Calculte ROHF-based SF-CIS MO-based density matrix using excitation
!  coefficients. nopen>=2 is required.
subroutine calc_sfcis_mo_dm_using_exc(nfc, nopen, nocc, nvir, exc, dm)
 implicit none
 integer :: i, j, p, q, nval, ndb, nmo
 integer, intent(in) :: nfc, nopen, nocc, nvir
 real(kind=8) :: r_ne
 real(kind=8), parameter :: ne_thres = 1d-2
 real(kind=8), intent(in) :: exc(nocc,nvir)
 real(kind=8), intent(out) :: dm(nfc+nocc+nvir-nopen,nfc+nocc+nvir-nopen)

 dm = 0d0; nval = nocc - nopen; ndb = nfc + nval; nmo = ndb + nvir
 if(nfc > 0) then
  do i = 1, nfc, 1
   dm(i,i) = 2d0
  end do ! for i
 end if

 do i = 1, nval, 1 ! d_ii, i <- {C}
  j = nfc + i
  dm(j,j) = 2d0 - DOT_PRODUCT(exc(i,:), exc(i,:))
 end do ! for i

!$omp parallel do schedule(dynamic) default(shared) private(i,j)
 do i = 1, nvir, 1 ! d_aa, a <- {O+V}
  j = ndb + i
  dm(j,j) = DOT_PRODUCT(exc(:,i), exc(:,i))
  if(i <= nopen) then
   dm(j,j) = dm(j,j) + 1d0 - DOT_PRODUCT(exc(nval+i,:), exc(nval+i,:))
  end if
 end do ! for i
!$omp end parallel do

 r_ne = 0d0
 do i = 1, nmo, 1
  r_ne = r_ne + dm(i,i)
 end do ! for i

 if(DABS(r_ne - DBLE(2*ndb+nopen)) > ne_thres) then
  write(6,'(/,A)') 'ERROR in subroutine calc_sfcis_mo_dm_using_exc: the number &
                   &of electrons'
  write(6,'(A)') 'calculated from exc differs from 2*ndb+nopen.'
  write(6,'(A,F12.4,2(A,I0))') 'r_ne=', r_ne, ', ndb=', ndb, ', nopen=', nopen
  stop
 end if

 do i = 1, nval-1, 1 ! d_ij, ij <- {C}
  p = nfc + i
  do j = i+1, nval, 1
   q = nfc + j
   dm(p,q) = -DOT_PRODUCT(exc(i,:), exc(j,:))
  end do ! for j
 end do ! for i

!$omp parallel do schedule(dynamic) default(shared) private(i,j,p,q)
 do i = 1, nvir-1, 1 ! d_ab, ab <- {O+V}
  p = ndb + i
  do j = i+1, nvir, 1
   q = ndb + j
   dm(p,q) = DOT_PRODUCT(exc(:,i), exc(:,j))
   if(i<=nopen .and. j<=nopen) then
    dm(p,q) = dm(p,q) - DOT_PRODUCT(exc(nval+i,:), exc(nval+j,:))
   end if
  end do ! for j
 end do ! for i
!$omp end parallel do

 do i = 1, nval, 1 ! d_is, i <- {C}, s <- {O}
  p = nfc + i
  do j = 1, nopen, 1
   q = ndb + j
   dm(p,q) = -DOT_PRODUCT(exc(i,:), exc(nval+j,:))
  end do ! for j
 end do ! for i

 ! SF-CIS d_ia = 0 for i<-{C}, a<-{V}, no need to calculate.
 call symmetrize_dmat(nmo, dm)
end subroutine calc_sfcis_mo_dm_using_exc

! Calculte ROHF-based SF-CIS <S^2> expectation value using excitation
!  coefficients. nopen>=2 is required.
subroutine calc_sfcis_ssquare_using_exc(nopen, nocc, nvir, exc, ssquare)
 implicit none
 integer :: a, b
 integer, intent(in) :: nopen, nocc, nvir
!f2py intent(in) :: nopen, nocc, nvir
 real(kind=8) :: rtmp
 real(kind=8), intent(in) :: exc(nocc,nvir)
!f2py intent(in) :: exc
!f2py depend(nocc,nvir) :: exc
 real(kind=8), intent(out) :: ssquare
!f2py intent(out) :: ssquare

 if(nopen < 2) then
  write(6,'(/,A)') 'ERROR in subroutine calc_sfcis_ssquare_using_exc: nopen>=2 &
                   &is required.'
  write(6,'(A,I0)') 'But got nopen=', nopen
  stop
 end if
 write(6,'(A,3I5)') 'nopen, nocc, nvir=', nopen, nocc, nvir
 ssquare = DBLE(nopen)

 do a = 1, nopen, 1
  ssquare = ssquare - DOT_PRODUCT(exc(:,a), exc(:,a))
  rtmp = exc(nocc-nopen+a,a)
  do b = 1, nopen, 1
   ssquare = ssquare + rtmp*exc(nocc-nopen+b,b)
  end do ! for b
 end do ! for a
end subroutine calc_sfcis_ssquare_using_exc

! Calculte ROHF-based SF-CIS MO-based transition density matrix using excitation
!  coefficients. nopen>=2 is required.
subroutine calc_sfcis_mo_tdm_using_exc(nopen, nocc, nvir, exc, tdm)
 implicit none
 integer :: i, j, p, q, ndb
 integer, intent(in) :: nopen, nocc, nvir
!f2py intent(in) :: nopen, nocc, nvir
 real(kind=8) :: r
 real(kind=8), intent(in) :: exc(nocc,nvir,2)
!f2py intent(in) :: exc
!f2py depend(nocc,nvir) :: exc
 real(kind=8), intent(out) :: tdm(nocc+nvir-nopen,nocc+nvir-nopen)
!f2py intent(out) :: tdm
!f2py depend(nocc,nvir,nopen) :: tdm

 tdm = 0d0; ndb = nocc - nopen; r = 0d0

!$omp parallel do schedule(dynamic) default(shared) private(i) reduction(+:r)
 do i = 1, nvir, 1
  r = r + DOT_PRODUCT(exc(:,i,1), exc(:,i,2))
 end do ! for i
!$omp end parallel do
 r = r*2d0

 do i = 1, ndb, 1 ! d_ii, i <- {C}
  tdm(i,i) = r - DOT_PRODUCT(exc(i,:,1), exc(i,:,2))
 end do ! for i
 r = 0.5d0*r

!$omp parallel do schedule(dynamic) default(shared) private(i,j)
 do i = 1, nvir, 1 ! d_aa, a <- {O+V}
  j = ndb + i
  tdm(j,j) = DOT_PRODUCT(exc(:,i,1), exc(:,i,2))
  if(i <= nopen) then
   tdm(j,j) = tdm(j,j) + r - DOT_PRODUCT(exc(j,:,1), exc(j,:,2))
  end if
 end do ! for i
!$omp end parallel do

 do i = 1, ndb-1, 1 ! d_ij, ij <- {C}
  do j = i+1, ndb, 1
   tdm(i,j) = -DOT_PRODUCT(exc(i,:,2), exc(j,:,1))
   tdm(j,i) = -DOT_PRODUCT(exc(j,:,2), exc(i,:,1))
  end do ! for j
 end do ! for i

!$omp parallel do schedule(dynamic) default(shared) private(i,j,p,q)
 do i = 1, nvir-1, 1 ! d_ab, ab <- {O+V}
  p = ndb + i
  do j = i+1, nvir, 1
   q = ndb + j
   tdm(p,q) = DOT_PRODUCT(exc(:,i,1), exc(:,j,2))
   tdm(q,p) = DOT_PRODUCT(exc(:,j,1), exc(:,i,2))
   if(i<=nopen .and. j<=nopen) then
    tdm(p,q) = tdm(p,q) - DOT_PRODUCT(exc(p,:,2), exc(q,:,1))
    tdm(q,p) = tdm(q,p) - DOT_PRODUCT(exc(q,:,2), exc(p,:,1))
   end if
  end do ! for j
 end do ! for i
!$omp end parallel do

 do i = 1, ndb, 1 ! d_is, i <- {C}, s <- {O}
  do j = 1, nopen, 1
   p = ndb + j
   tdm(i,p) = -DOT_PRODUCT(exc(i,:,2), exc(p,:,1))
   tdm(p,i) = -DOT_PRODUCT(exc(p,:,2), exc(i,:,1))
  end do ! for j
 end do ! for i

 ! SF-CIS tdm_ia = 0 for i<-{C}, a<-{V}, no need to calculate.
end subroutine calc_sfcis_mo_tdm_using_exc

! Calculte ROHF-based MRSF-CIS MO-based density matrix using excitation
!  coefficients. nopen=2 is required since MRSF supports only the triplet
!  reference currently.
subroutine calc_mrsfcis_mo_dm_using_exc(nfc, nopen, nocc, nvir, exc, dm)
 implicit none
 integer :: i, j, p, q, nval, ndb, nmo
 integer, intent(in) :: nfc, nopen, nocc, nvir
 real(kind=8) :: r_ne
 real(kind=8), parameter :: ne_thres = 1d-2
 real(kind=8), intent(in) :: exc(nocc,nvir)
 real(kind=8), intent(out) :: dm(nfc+nocc+nvir-nopen,nfc+nocc+nvir-nopen)

 if(nopen /= 2) then
  write(6,'(/,A)') 'ERROR in subroutine calc_mrsfcis_mo_dm_using_exc: currently&
                   & only triplet'
  write(6,'(A,I0)') 'ROHF/ROKS reference is supported, i.e. nopen=2. But got no&
                    &pen=', nopen
  stop
 end if
 dm = 0d0; nval = nocc - nopen; ndb = nfc + nval; nmo = ndb + nvir

 if(nfc > 0) then
  forall(i = 1:nfc) dm(i,i) = 2d0
 end if

 do i = 1, nval, 1 ! d_ii, i <- {C}
  j = nfc + i
  dm(j,j) = 2d0 - DOT_PRODUCT(exc(i,:),exc(i,:))
 end do ! for i

 do i = 1, nval-1, 1 ! d_ij, ij <- {C}
  p = nfc + i
  do j = i+1, nval, 1
   q = nfc + j
   dm(p,q) = -DOT_PRODUCT(exc(i,:),exc(j,:))
  end do ! for j
 end do ! for i

 ! d_ss
 dm(ndb+1,ndb+1) = 1d0 + DOT_PRODUCT(exc(:,1),exc(:,1)) - &
                         DOT_PRODUCT(exc(nval+1,:),exc(nval+1,:))
 ! d_tt
 dm(ndb+2,ndb+2) = 1d0 + DOT_PRODUCT(exc(:,2),exc(:,2)) - &
                         DOT_PRODUCT(exc(nocc,:),exc(nocc,:))
 ! d_st
 dm(ndb+1,ndb+2) = (exc(nocc,1)-exc(nval+1,2))*(exc(nval+1,1)+exc(nocc,2)) + &
  DOT_PRODUCT(exc(:,1),exc(:,2)) - DOT_PRODUCT(exc(nval+1,:),exc(nocc,:))

!$omp parallel do schedule(dynamic) default(shared) private(i,j,p,q)
 do i = 1, nvir-nopen, 1 ! d_ab, ab <- {V}, including a=b
  p = nfc + nocc + i
  do j = i, nvir-nopen, 1
   q = nfc + nocc + j
   dm(p,q) = DOT_PRODUCT(exc(:,nopen+i),exc(:,nopen+j))
  end do ! for j
 end do ! for i
!$omp end parallel do

 r_ne = 0d0
 do i = 1, nmo, 1
  r_ne = r_ne + dm(i,i)
 end do ! for i

 if(DABS(r_ne - DBLE(2*ndb+nopen)) > ne_thres) then
  write(6,'(/,A)') 'ERROR in subroutine calc_mrsfcis_mo_dm_using_exc: the numbe&
                   &r of electrons'
  write(6,'(A)') 'calculated from exc differs from 2*ndb+nopen.'
  write(6,'(A,F12.4,2(A,I0))') 'r_ne=', r_ne, ', ndb=', ndb, ', nopen=', nopen
  stop
 end if

 do i = 1, nval, 1 ! d_is, i <- {C}, s <- {O}
  p = nfc + i
  dm(p,ndb+1) = exc(i,2)*exc(nval+1,2) - exc(i,1)*exc(nocc,2) - &
                DOT_PRODUCT(exc(i,:),exc(nval+1,:))
  dm(p,ndb+2) = exc(i,1)*exc(nocc,1) - exc(i,2)*exc(nval+1,1) - &
                DOT_PRODUCT(exc(i,:),exc(nocc,:))
 end do ! for i

!$omp parallel do schedule(dynamic) default(shared) private(i,p)
 do i = 1, nvir-nopen, 1 ! d_sa, s <- {O}, a <- {V}
  p = nfc + nocc + i
  dm(ndb+1,p) = exc(nocc,nopen+i)*exc(nocc,1) - exc(nval+1,nopen+i)*exc(nocc,2) + &
                DOT_PRODUCT(exc(:,1),exc(:,nopen+i))
  dm(ndb+2,p) = exc(nval+1,nopen+i)*exc(nval+1,2) - exc(nocc,nopen+i)*exc(nval+1,1) + &
                DOT_PRODUCT(exc(:,2),exc(:,nopen+i))
 end do ! for i
!$omp end parallel do

 ! MRSF-CIS d_ia = 0 for i<-{C}, a<-{V}, no need to calculate.
 call symmetrize_dmat(nmo, dm)
end subroutine calc_mrsfcis_mo_dm_using_exc

! Read NROrb, NOA, NOB from a Gaussian output file.
! NFC          : number of frozen doubly occupied orbitals
! NROrb(nif_ex): number of spatial orbitals involved in excitation, nif_ex<=nif
! NOA     (noa): number of alpha occupied spin orbitals involved in excitation
! NOB     (nob): number of beta occupied spin orbitals involved in excitation
subroutine read_noa_nob_from_gau_log(logname, nfc, nif_ex, noa, nob)
 implicit none
 integer :: i, fid
 integer, intent(out) :: nfc, nif_ex, noa, nob
!f2py intent(out) :: nfc, nif_ex, noa, nob
 character(len=240) :: buf
 character(len=240), intent(in) :: logname
!f2py intent(in) :: logname

 open(newunit=fid,file=TRIM(logname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:9) == 'Range of') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_noa_nob_from_gau_log: no 'Range' f&
                   &ound in file "//TRIM(logname)
  close(fid)
  stop
 end if

 read(fid,'(A)') buf
 i = INDEX(buf, 'NFC=')
 read(buf(i+4:),*) nfc

 read(fid,'(A)') buf
 close(fid)
 i = INDEX(buf, '=')
 read(buf(i+1:),*) nif_ex
 buf(i:i) = ' '

 i = INDEX(buf, '=')
 read(buf(i+1:),*) noa
 buf(i:i) = ' '

 i = INDEX(buf, '=')
 read(buf(i+1:),*) nob
end subroutine read_noa_nob_from_gau_log

! Check whether noa and nob is consistent with input nocc, nvir.
! If consistent, return nfc, nif_ex; if not, stop
subroutine check_noa_nob_in_logname(logname, alpha, nocc, nvir, nfc, nif_ex)
 implicit none
 integer :: noa, nob, nva, nvb
 integer, intent(in) :: nocc, nvir
!f2py intent(in) :: nocc, nvir
 integer, intent(out) :: nfc, nif_ex
!f2py intent(out) :: nfc, nif_ex
 logical, intent(in) :: alpha
!f2py intent(in) :: alpha
 character(len=240), intent(in) :: logname
!f2py intent(in) :: logname

 call read_noa_nob_from_gau_log(logname, nfc, nif_ex, noa, nob)
 nva = nif_ex - noa
 nvb = nif_ex - nob

 if(alpha) then ! check alpha excitation information
  if(nocc/=noa .or. nvir/=nva) then
   write(6,'(/,A)') 'ERROR in subroutine check_noa_nob_in_logname: nocc/=noa or&
                    & nvir/=nva.'
   write(6,'(A)') 'logname='//TRIM(logname)
   write(6,'(A,4I5)') 'nocc, noa,  nvir, nva=', nocc, noa,  nvir, nva
   stop
  end if
 else ! check beta excitation information
  if(nocc/=nob .or. nvir/=nvb) then
   write(6,'(/,A)') 'ERROR in subroutine check_noa_nob_in_logname: nocc/=nob or&
                    & nvir/=nvb.'
   write(6,'(A)') 'logname='//TRIM(logname)
   write(6,'(A,4I5)') 'nocc, nob,  nvir, nvb=', nocc, nob,  nvir, nvb
   stop
  end if
 end if
end subroutine check_noa_nob_in_logname

! Read CIS/TDHF/TDDFT excitation coefficients (also called amplitudes) from a
!  specified Gaussian output file. This subroutine is only used for RHF-based
!  CIS/TDHF/TDDFT. For UHF-based ones, please use read_ex_coeff_from_gau_log2
!  below.
! Note: the sum of square of elements in exc is 0.5. There are two explanations
!  for this number: 1) these are only Alpha CI coefficients; 2) these are scaled
!  by sqrt(2) for an unknown reason by Gaussian developers.
subroutine read_ex_coeff_from_gau_log(logname, istate, nocc, nvir, exc)
 implicit none
 integer :: i, j, k, fid, nfc, nif_ex
 integer, intent(in) :: istate, nocc, nvir
!f2py intent(in) :: istate, nocc, nvir
! nocc: number of alpha occupied orbitals involved in excitation
! nvir: number of beta occupied orbitals involved in excitation
 real(kind=8) :: rtmp
 real(kind=8), intent(out) :: exc(nocc,nvir)
!f2py intent(out) :: exc
!f2py depend(nocc,nvir) :: exc
 real(kind=8), parameter :: diff_thres = 0.01d0
 character(len=17) :: str
 character(len=240) :: buf
 character(len=240), intent(in) :: logname
!f2py intent(in) :: logname

 exc = 0d0
 call check_noa_nob_in_logname(logname, .true., nocc, nvir, nfc, nif_ex)

 write(str,'(A13,I4)') 'Excited State', istate
 open(newunit=fid,file=TRIM(logname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:18) == str) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_ex_coeff_from_gau_log: no '"//&
                    str//"' found"
  write(6,'(A)') 'in file '//TRIM(logname)
  close(fid)
  stop
 end if

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:11)=='This state' .or. buf(2:7)=='SavETr' .or. LEN_TRIM(buf)==0) exit
  i = INDEX(buf,'->')
  if(i == 0) i = INDEX(buf,'<-') ! de-excitation
  if(i == 0) then
   close(fid)
   write(6,'(/,A)') "ERROR in subroutine read_ex_coeff_from_gau_log: no '->' or&
                    & '<-' symbol"
   write(6,'(A)') 'found in file '//TRIM(logname)
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
 rtmp = 0d0
 do i = 1, nvir, 1
  rtmp = rtmp + DOT_PRODUCT(exc(:,i), exc(:,i))
 end do ! for i

 ! check the difference with 0.5
 if(DABS(rtmp - 0.5d0) > diff_thres) then
  write(6,'(/,A)') 'ERROR in subroutine read_ex_coeff_from_gau_log: the sum of &
                   &|C_ia|^2 differs'
  write(6,'(A)') 'from 0.5. Did you forget to write IOP(9/40=4) in gjf?'
  stop
 end if
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
  i = INDEX(buf,'->')
  if(i == 0) i = INDEX(buf,'<-') ! de-excitation
  if(i == 0) then
   close(fid)
   write(6,'(A)') "ERROR in subroutine read_ex_coeff_from_gau_log2: no '->'&
                  & or '<-' symbol"
   write(6,'(A)') 'found in file '//TRIM(logname)
   stop
  end if

  read(buf(1:i-1),*) str6
  call del_a_or_b(str6)
  read(str6,*) j
  if(index(str,'A') > 0) then
   alpha = .true.
  else
   alpha = .false.
  end if

  read(buf(i+2:),*) str6
  call del_a_or_b(str6)
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

! Read SF-/MRSF- CIS/TD excitation coefficients (also called amplitudes) from a
!  specified GAMESS output file (.gms).
! Note:
! 1) istate=0 means the ground state 'STATE #   1' in the .gms file.
! 2) nval is the number of valence doubly occupied orbitals in nocc. Frozen
!    core orbitals are neither considered in nval, nor in nocc.
subroutine read_sf_ex_coeff_from_gms_gms(gmsname, istate, nval, nocc, nvir, exc)
 implicit none
 integer :: i, j, k, fid
 integer, intent(in) :: istate, nval, nocc, nvir
!f2py intent(in) :: istate, nval, nocc, nvir
 real(kind=8) :: rtmp
 real(kind=8), intent(out) :: exc(nocc,nvir)
!f2py intent(out) :: exc
!f2py depend(nocc,nvir) :: exc
 real(kind=8), parameter :: diff_thres = 0.01d0
 character(len=2) :: str
 character(len=11) :: key
 character(len=240) :: buf
 character(len=240), intent(in) :: gmsname
!f2py intent(in) :: gmsname

 exc = 0d0
 write(key,'(A7,I4)') 'STATE #', istate+1
 open(newunit=fid,file=TRIM(gmsname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:12) == key) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_sf_ex_coeff_from_gms_gms: no '"//&
                    key//"'"
  write(6,'(A)') 'is found in file '//TRIM(gmsname)
  close(fid)
  stop
 end if

 do i = 1, 4 ! skip 4 lines
  read(fid,'(A)') buf
 end do ! for i

 do while(.true.)
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
  read(buf,*) k, rtmp, i, str, j
  exc(i,j-nval) = rtmp
 end do ! for while

 close(fid)
 rtmp = 0d0
 do i = 1, nvir, 1
  rtmp = rtmp + DOT_PRODUCT(exc(:,i),exc(:,i))
 end do ! for i

 if(DABS(rtmp - 1d0) > diff_thres) then
  write(6,'(/,A)') 'ERROR in subroutine read_sf_ex_coeff_from_gms_gms: the sum &
                   &of |C_ia|^2'
  write(6,'(A)') 'differs from 1.0. Did you forget to modify GAMESS code so as &
                 &to print more'
  write(6,'(A)') 'excitation coefficients?'
  stop
 end if
end subroutine read_sf_ex_coeff_from_gms_gms

! Generate CIS/TD MO-based density matrix using excitation coefficients from a
!  Gaussian output file. The reference can only be RHF/RKS.
! Note:
! 1) Iop(9/40=4) is recommended when generating the log/out file.
! 2) the calculation of CIS/TD MO-based dm requires only the excitation
!  coefficients (also called amplitudes sometimes). The MO coefficients are not
!  needed. For CIS/TD AO-based dm, see subroutine gen_cis_ao_dm_from_fch_and_log
!  below.
subroutine gen_cis_mo_dm_from_gau_log(logname, istate, averaged, nif, mo_dm)
 implicit none
 integer :: i, nfc, nif_ex, nocc, nob, ndb, nvir
 integer, intent(in) :: istate, nif
!f2py intent(in) :: istate, nif
 real(kind=8), intent(out) :: mo_dm(nif,nif)
!f2py intent(out) :: mo_dm
!f2py depend(nif) :: mo_dm
 real(kind=8), allocatable :: mo_dm0(:,:), exc(:,:)
 character(len=240), intent(in) :: logname
!f2py intent(in) :: logname
 logical, intent(in) :: averaged
!f2py intent(in) :: averaged
 ! True: generate MO-based averaged density matrix of all 0~i states
 ! False: only generate MO-based density matrix of the i-th excited state

 if(istate < 1) then
  write(6,'(/,A)') 'ERROR in subroutine gen_cis_mo_dm_from_gau_log: istate>=1 i&
                   &s required.'
  write(6,'(A)') 'logname='//TRIM(logname)
  write(6,'(A,I0)') 'istate=', istate
  stop
 end if

 call read_noa_nob_from_gau_log(logname, nfc, nif_ex, nocc, nob)
 if(nif-nfc-nif_ex /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine gen_cis_mo_dm_from_gau_log: inconsisten&
                   &cy detected in'
  write(6,'(A)') 'file '//TRIM(logname)
  write(6,'(A,3I6)') 'nif, nfc, nif_ex=', nif, nfc, nif_ex
  stop
 end if

 ndb = nfc + nocc; nvir = nif_ex - nocc
 allocate(exc(nocc,nvir))

 if(averaged) then
  write(6,'(A,I0,A)') 'State-averaged NO generated using states 0~',istate,'.'
  call init_rhf_occ_mat(ndb, nif, mo_dm) ! RHF occupation numbers
  allocate(mo_dm0(nif,nif))
  do i = 1, istate, 1
   call read_ex_coeff_from_gau_log(logname, i, nocc, nvir, exc)
   exc = DSQRT(2d0)*exc  ! sum of square of exc is 0.5, so multiply root2
   call calc_cis_mo_dm_using_exc(nfc, nocc, nvir, exc, mo_dm0)
   mo_dm = mo_dm + mo_dm0
  end do ! for i
  deallocate(mo_dm0)
  mo_dm = mo_dm/DBLE(istate+1)

 else ! a specific excited state
  call read_ex_coeff_from_gau_log(logname, istate, nocc, nvir, exc)
  exc = DSQRT(2d0)*exc  ! sum of square of exc is 0.5, so multiply root2
  call calc_cis_mo_dm_using_exc(nfc, nocc, nvir, exc, mo_dm)
 end if

 deallocate(exc)
 ! The diagonalization of mo_dm leads to CIS/TD natural orbital occupation
 ! numbers, but we will not do it here.
end subroutine gen_cis_mo_dm_from_gau_log

! Generate SF-CIS/TD MO-based density matrix using excitation coefficients
!  from a GAMESS output file (.gms). The spin multiplicity of the ROHF/ROKS
!  reference must be >= 3.
subroutine gen_sfcis_mo_dm_from_gms_gms(gmsname, averaged, istate, nfc, nval, &
                                        nif, mo_dm)
 implicit none
 integer :: i, mult, nopen, nocc, nvir
 integer, intent(in) :: istate, nfc, nval, nif
!f2py intent(in) :: istate, nfc, nval, nif
 real(kind=8), intent(out) :: mo_dm(nif,nif)
!f2py intent(out) :: mo_dm
!f2py depend(nif) :: mo_dm
 real(kind=8), allocatable :: mo_dm0(:,:), exc(:,:)
 character(len=240), intent(in) :: gmsname
!f2py intent(in) :: gmsname
 logical, intent(in) :: averaged
!f2py intent(in) :: averaged

 if(istate < 0) then
  write(6,'(/,A)') 'ERROR in subroutine gen_sfcis_mo_dm_from_gms_gms: istate>=0&
                   & is required.'
  write(6,'(A)') 'gmsname='//TRIM(gmsname)
  write(6,'(A,I0)') 'istate=', istate
  stop
 end if

 call read_mult_from_gms_gms(gmsname, mult)
 if(mult < 3) then
  write(6,'(/,A)') 'ERROR in subroutine gen_sfcis_mo_dm_from_gms_gms: mult>=3 i&
                   &s required for '//TRIM(gmsname)
  write(6,'(A,I0)') 'But detected mult=', mult
  stop
 end if

 nopen = mult - 1; nocc = nval + nopen; nvir = nif - nfc - nval
 allocate(exc(nocc,nvir))

 if(averaged) then
  write(6,'(A,I0,A)') 'State-averaged SF-CIS NO generated using states 0~', &
                      istate, '.'
  mo_dm = 0d0
  allocate(mo_dm0(nif,nif))
  do i = 0, istate, 1
   call read_sf_ex_coeff_from_gms_gms(gmsname, i, nval, nocc, nvir, exc)
   call calc_sfcis_mo_dm_using_exc(nfc, nopen, nocc, nvir, exc, mo_dm0)
   mo_dm = mo_dm + mo_dm0
  end do ! for i
  deallocate(mo_dm0)
  mo_dm = mo_dm/DBLE(istate+1)
 else ! a specific electronic state
  call read_sf_ex_coeff_from_gms_gms(gmsname, istate, nval, nocc, nvir, exc)
  call calc_sfcis_mo_dm_using_exc(nfc, nopen, nocc, nvir, exc, mo_dm)
 end if

 deallocate(exc)
end subroutine gen_sfcis_mo_dm_from_gms_gms

! Generate MRSF-CIS/TD MO-based density matrix using excitation coefficients
!  from a GAMESS output file (.gms). The spin multiplicity of the ROHF/ROKS
!  reference can only be 3.
subroutine gen_mrsfcis_mo_dm_from_gms_gms(gmsname, averaged, istate, nfc, nval,&
                                          nif, mo_dm)
 implicit none
 integer :: i, mult, nopen, nocc, nvir
 integer, intent(in) :: istate, nfc, nval, nif
!f2py intent(in) :: istate, nfc, nval, nif
 real(kind=8), intent(out) :: mo_dm(nif,nif)
!f2py intent(out) :: mo_dm
!f2py depend(nif) :: mo_dm
 real(kind=8), allocatable :: mo_dm0(:,:), exc(:,:)
 character(len=240), intent(in) :: gmsname
!f2py intent(in) :: gmsname
 logical, intent(in) :: averaged
!f2py intent(in) :: averaged

 if(istate < 0) then
  write(6,'(/,A)') 'ERROR in subroutine gen_mrsfcis_mo_dm_from_gms_gms: istate>&
                   &=0 is required.'
  write(6,'(A)') 'gmsname='//TRIM(gmsname)
  write(6,'(A,I0)') 'istate=', istate
  stop
 end if

 call read_mult_from_gms_gms(gmsname, mult)
 if(mult /= 3) then
  write(6,'(/,A)') 'ERROR in subroutine gen_mrsfcis_mo_dm_from_gms_gms: MRSF ca&
                   &n only be applied'
  write(6,'(A)') 'in the triplet case currently. But a different spin multiplic&
                 &ity is detected in'
  write(6,'(A)') 'file '//TRIM(gmsname)
  write(6,'(A,I0)') 'mult=', mult
  stop
 end if

 nopen = mult - 1; nocc = nval + nopen; nvir = nif - nfc - nval
 allocate(exc(nocc,nvir))

 if(averaged) then
  write(6,'(A,I0,A)') 'State-averaged MRSF-CIS NO generated using states 0~', &
                      istate, '.'
  mo_dm = 0d0
  allocate(mo_dm0(nif,nif))
  do i = 0, istate, 1
   call read_sf_ex_coeff_from_gms_gms(gmsname, i, nval, nocc, nvir, exc)
   call calc_mrsfcis_mo_dm_using_exc(nfc, nopen, nocc, nvir, exc, mo_dm0)
   mo_dm = mo_dm + mo_dm0
  end do ! for i
  deallocate(mo_dm0)
  mo_dm = mo_dm/DBLE(istate+1)
 else ! a specific electronic state
  call read_sf_ex_coeff_from_gms_gms(gmsname, istate, nval, nocc, nvir, exc)
  call calc_mrsfcis_mo_dm_using_exc(nfc, nopen, nocc, nvir, exc, mo_dm)
 end if

 deallocate(exc)
end subroutine gen_mrsfcis_mo_dm_from_gms_gms

! Generate MRSF-CIS/TD State-averaged Mixed-spin MO-based density matrix using
!  excitation coefficients from two GAMESS output files (.gms). The spin multi-
!  plicity of the ROHF/ROKS reference must be 3.
! GAMESS MRSF either calculate singlets or triplets at one time. So singlets/
!  triplets excitation coefficients are stored/printed in two .gms files.
! gmsname1/gmsname2: singlet/triplet states
subroutine gen_mrsfcis_ave_mo_dm_from_gms_gms(gmsname1, gmsname2, istate1, &
                                              istate2, nfc, nval, nif, mo_dm)
 implicit none
 integer :: i, mult(2), nocc, nvir
 integer, intent(in) :: istate1, istate2, nfc, nval, nif
!f2py intent(in) :: istate1, istate2, nfc, nval, nif
 real(kind=8), intent(out) :: mo_dm(nif,nif)
!f2py intent(out) :: mo_dm
!f2py depend(nif) :: mo_dm
 real(kind=8), allocatable :: mo_dm0(:,:), exc(:,:)
 character(len=240), intent(in) :: gmsname1, gmsname2
!f2py intent(in) :: gmsname1, gmsname2

 if(istate1<0 .or. istate2<0) then
  write(6,'(/,A)') 'ERROR in subroutine gen_mrsfcis_ave_mo_dm_from_gms_gms: ist&
                   &ate>=0 is required.'
  write(6,'(A)') 'gmsname1='//TRIM(gmsname1)//', gmsname2='//TRIM(gmsname2)
  write(6,'(A,2I8)') 'istate1, istate2=', istate1, istate2
  stop
 end if

 call read_mult_from_gms_gms(gmsname1, mult(1))
 call read_mult_from_gms_gms(gmsname2, mult(2))
 if(ANY(mult /= 3)) then
  write(6,'(/,A)') 'ERROR in subroutine gen_mrsfcis_ave_mo_dm_from_gms_gms: MRS&
                   &F can only be'
  write(6,'(A)') 'applied for triplet reference. Please check files:'
  write(6,'(A)') 'gmsname1='//TRIM(gmsname1)//', gmsname2='//TRIM(gmsname2)
  write(6,'(A,2I8)') 'mult=', mult
  stop
 end if

 write(6,'(A)') 'State-averaged mixed-spin MRSF-CIS NO generated using'
 write(6,'(A,I0,A)') 'states 0~',istate1, ' from file '//TRIM(gmsname1)
 write(6,'(A,I0,A)') 'states 0~',istate2, ' from file '//TRIM(gmsname2)

 nocc = nval + 2; nvir = nif - nfc - nval; mo_dm = 0d0
 allocate(exc(nocc,nvir), mo_dm0(nif,nif))

 do i = 0, istate1, 1
  call read_sf_ex_coeff_from_gms_gms(gmsname1, i, nval, nocc, nvir, exc)
  call calc_mrsfcis_mo_dm_using_exc(nfc, 2, nocc, nvir, exc, mo_dm0)
  mo_dm = mo_dm + mo_dm0
 end do ! for i

 do i = 0, istate2, 1
  call read_sf_ex_coeff_from_gms_gms(gmsname2, i, nval, nocc, nvir, exc)
  call calc_mrsfcis_mo_dm_using_exc(nfc, 2, nocc, nvir, exc, mo_dm0)
  mo_dm = mo_dm + mo_dm0
 end do ! for i

 deallocate(exc, mo_dm0)
 mo_dm = mo_dm/DBLE(istate1+istate2+2)
end subroutine gen_mrsfcis_ave_mo_dm_from_gms_gms

! Generate RHF/RKS-based CIS/TD (state-averaged) NOs using Gaussian
!  .fch(k) and .log files. Canonical MOs must be included in fchname. (SA-)NOs
!  and correponding density will be saved into xxx_NO.fch.
! Note: Iop(9/40=4) is recommended when generating the log/out file.
subroutine gen_cis_no_from_fch_and_log(fchname, logname, istate, averaged)
 implicit none
 integer :: i, nbf, nif
 integer, intent(in) :: istate
!f2py intent(in) :: istate
 real(kind=8), allocatable :: ao_dm(:,:)
 real(kind=8), allocatable :: mo(:,:), mo_dm(:,:)
 character(len=240) :: no_fch
 character(len=240), intent(in) :: fchname, logname
!f2py intent(in) :: fchname, logname
 logical, intent(in) :: averaged
!f2py intent(in) :: averaged

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(mo_dm(nif,nif))
 call gen_cis_mo_dm_from_gau_log(logname, istate, averaged, nif, mo_dm)
 allocate(mo(nbf,nif))
 call read_mo_from_fch(fchname, nbf, nif, 'a', mo)

 ! P = Cn(C^T), where n is mo_dm below and usually not diagonal.
 allocate(ao_dm(nbf,nbf))
 call calc_cnct(nbf, nif, mo, mo_dm, ao_dm)
 deallocate(mo, mo_dm)

 call find_specified_suffix(fchname, '.fch', i)
 no_fch = fchname(1:i-1)//'_NO.fch'
 call copy_file(fchname, no_fch, .false.)

 call write_dm_into_fch(no_fch, .true., nbf, ao_dm)
 deallocate(ao_dm)
 call gen_no_using_dm_in_fch(no_fch, 1)
end subroutine gen_cis_no_from_fch_and_log

! Generate ROHF/ROKS-based SF-/MRSF- CIS/TD (state-averaged) NOs using Gaussian
! .fch(k) and GAMESS .gms files. Canonical MOs must be included in fchname.
! (SA-)NOs and correponding density will be saved into xxx_NO.fch. We cannot
! read nfc and nval from .fch and .gms, so these two variables are provided as
! input.
subroutine gen_sfcis_no_from_fch_and_gms(fchname, gmsname, istate, nfc, nval, &
                                         mrsf, averaged)
 implicit none
 integer :: i, nbf, nif
 integer, intent(in) :: istate, nfc, nval
!f2py intent(in) :: istate, nfc, nval
 real(kind=8), allocatable :: ao_dm(:,:)
 real(kind=8), allocatable :: mo(:,:), mo_dm(:,:)
 character(len=240) :: no_fch
 character(len=240), intent(in) :: fchname, gmsname
!f2py intent(in) :: fchname, gmsname
 logical, intent(in) :: mrsf, averaged
!f2py intent(in) :: mrsf, averaged

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(mo_dm(nif,nif))

 if(mrsf) then
  call gen_mrsfcis_mo_dm_from_gms_gms(gmsname, averaged, istate, nfc, nval, &
                                      nif, mo_dm)
 else
  call gen_sfcis_mo_dm_from_gms_gms(gmsname, averaged, istate, nfc, nval, nif, &
                                    mo_dm)
 end if

 allocate(mo(nbf,nif))
 call read_mo_from_fch(fchname, nbf, nif, 'a', mo)

 ! P = Cn(C^T), where n is mo_dm below and usually not diagonal.
 allocate(ao_dm(nbf,nbf))
 call calc_cnct(nbf, nif, mo, mo_dm, ao_dm)
 deallocate(mo, mo_dm)

 call find_specified_suffix(fchname, '.fch', i)
 no_fch = fchname(1:i-1)//'_NO.fch'
 call copy_file(fchname, no_fch, .false.)

 call write_dm_into_fch(no_fch, .true., nbf, ao_dm)
 deallocate(ao_dm)
 call gen_no_using_dm_in_fch(no_fch, 1)
end subroutine gen_sfcis_no_from_fch_and_gms

! Generate ROHF/ROKS-based MRSF-CIS/TD state-averaged mixed-spin NOs using
!  Gaussian .fch(k) and two GAMESS .gms files. For two states with the same
!  spin, please use subroutine gen_sfcis_no_from_fch_and_gms above.
subroutine gen_mrsfcis_sa_no_from_fch_and_gms(fchname, gmsname1, gmsname2, &
                                              istate1, istate2, nfc, nval)
 implicit none
 integer :: i, nbf, nif
 integer, intent(in) :: istate1, istate2, nfc, nval
!f2py intent(in) :: istate1, istate2, nfc, nval
 real(kind=8), allocatable :: ao_dm(:,:)
 real(kind=8), allocatable :: mo(:,:), mo_dm(:,:)
 character(len=240) :: no_fch
 character(len=240), intent(in) :: fchname, gmsname1, gmsname2
!f2py intent(in) :: fchname, gmsname1, gmsname2

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(mo_dm(nif,nif))
 call gen_mrsfcis_ave_mo_dm_from_gms_gms(gmsname1, gmsname2, istate1, istate2, &
                                         nfc, nval, nif, mo_dm)
 allocate(mo(nbf,nif))
 call read_mo_from_fch(fchname, nbf, nif, 'a', mo)

 ! P = Cn(C^T), where n is mo_dm below and usually not diagonal.
 allocate(ao_dm(nbf,nbf))
 call calc_cnct(nbf, nif, mo, mo_dm, ao_dm)
 deallocate(mo, mo_dm)

 call find_specified_suffix(fchname, '.fch', i)
 no_fch = fchname(1:i-1)//'_NO.fch'
 call copy_file(fchname, no_fch, .false.)

 call write_dm_into_fch(no_fch, .true., nbf, ao_dm)
 deallocate(ao_dm)
 call gen_no_using_dm_in_fch(no_fch, 1)
end subroutine gen_mrsfcis_sa_no_from_fch_and_gms

! Generate ground state -> singlet excited state (state-averaged) NTO for RHF/
!  RKS- CIS/TDHF/TDDFT wave function. Canonical MOs must be kept in fchname. NTOs
!  will be punched into file xxx_NTO.fch
! Note: For averaged=.True., logname must include the 1~`istate`-th singlet exc-
!  ited state results. For averaged=.False., logname must include the `istate`-th
!  singlet excited state result. Do not use this subroutine for S0->Tn, since
!  that transition density is exactly zero for RHF-CIS.
! T*(T^) = U(n_occ)(U^T) -> occupied NTOs
! (T^)*T = V(n_vir)(V^T) -> unoccupied NTOs
subroutine gen_cis_ge_nto_from_fch_and_log(fchname, logname, istate, averaged)
 implicit none
 integer :: i, nfc, nif_ex, nocc, nob, nvir, ndb, nbf, nif
 integer, intent(in) :: istate
!f2py inteny(in) :: istate
 real(kind=8), allocatable :: exc(:,:), exc0(:,:), ev(:), mo(:,:), occ_mo(:,:),&
  vir_mo(:,:), u(:,:), vt(:,:), s(:)
 character(len=240) :: nto_fch
 character(len=240), intent(in) :: fchname, logname
!f2py inteny(in) :: fchname, logname
 ! logname: it includes CIS excitation coefficients in .log file
 ! fchname: it includes canonical MOs in .fch file
 logical, intent(in) :: averaged
!f2py inteny(in) :: averaged
 ! True: State-Averaged NTO; False: NTO of the ground state -> an excited state

 call read_noa_nob_from_gau_log(logname, nfc, nif_ex, nocc, nob)
 ndb = nfc + nocc
 nvir = nif_ex - nocc
 allocate(exc(nocc,nvir))

 if(averaged) then
  write(6,'(A,I0,A)') 'State-averaged NTO generated using ', istate, &
                      ' excited states.'
  exc = 0d0
  allocate(exc0(nocc,nvir))
  ! The RHF->CIS (ground state -> singlet excited state) transition density
  ! matrix is simply DSQRT(2d0)*C_ia for spatial orbital basis (or C_ia for
  ! spin orbital basis), so the averaged transition density matrix can be
  ! calculated by simply summing different exc0.
  do i = 1, istate, 1
   call read_ex_coeff_from_gau_log(logname, i, nocc, nvir, exc0)
   exc = exc + exc0
  end do ! for i
  deallocate(exc0)
  exc = DSQRT(2d0)*exc/DSQRT(DBLE(istate))
 else ! ground state -> an excited state
  call read_ex_coeff_from_gau_log(logname, istate, nocc, nvir, exc)
  exc = DSQRT(2d0)*exc ! ! sum of square of exc(:,:) is 0.5, so multiply sqrt(2)
 end if

 ! The alpha/beta spin NTOs have identical MO coefficients and singular values.
 ! Only one set of them is calculated. In such a set of spin NTOs, the sum of
 ! square of singular values is 1. If you consider the spatial orbital basis,
 ! the sum would be 2. Nowadays people usually use spin NTOs when they talk about
 ! NTOs. So we do not multiply exc by sqrt(2) below.
 allocate(u(nocc,nocc), vt(nvir,nvir), s(nocc))
 call do_svd(nocc, nvir, exc, u, vt, s)
 deallocate(exc)
 do i = 1, nocc, 1
  s(i) = s(i)*s(i)
 end do ! for i

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(ev(nif), source=0d0)
 allocate(vir_mo(nocc,nocc))

 do i = 1, nocc, 1
  ev(nfc+i) = s(nocc-i+1)
  ev(ndb+i) = s(i)
  vir_mo(:,i) = u(:,nocc-i+1)
 end do ! for i

 deallocate(s)
 u = vir_mo
 deallocate(vir_mo)

 allocate(mo(nbf,nif))
 call read_mo_from_fch(fchname, nbf, nif, 'a', mo)
 allocate(occ_mo(nbf,nocc), source=0d0)
 ! C'_occ = (C_occ)U
 call dgemm('N', 'N', nbf, nocc, nocc, 1d0, mo(:,nfc+1:ndb), nbf, u, nocc, &
            0d0, occ_mo, nbf)
 mo(:,nfc+1:ndb) = occ_mo
 deallocate(occ_mo, u)

 allocate(vir_mo(nbf,nvir), source=0d0)
 ! C'_vir = (C_vir)V
 call dgemm('N', 'T', nbf, nvir, nvir, 1d0, mo(:,ndb+1:ndb+nvir), nbf, vt, &
            nvir, 0d0, vir_mo, nbf)
 mo(:,ndb+1:ndb+nvir) = vir_mo
 deallocate(vir_mo, vt)

 call find_specified_suffix(fchname, '.fch', i)
 nto_fch = fchname(1:i-1)//'_NTO.fch'
 call copy_file(fchname, nto_fch, .false.)

 call write_mo_into_fch(nto_fch, nbf, nif, 'a', mo)
 deallocate(mo)
 call write_eigenvalues_to_fch(nto_fch, nif, 'a', ev, .true.)
 deallocate(ev)
end subroutine gen_cis_ge_nto_from_fch_and_log

subroutine gen_cis_ee_nto_from_fch_and_log(fchname, logname, istate1, istate2)
 implicit none
 integer :: i, nfc, nif_ex, nocc, nob, nvir, ndb, nbf, nif
 integer, intent(in) :: istate1, istate2
!f2py intent(in) :: istate1, istate2
 real(kind=8), allocatable :: exc(:,:,:), tdm_occ(:,:), tdm_vir(:,:), ev(:), &
  u(:,:), vt(:,:), s(:), mo(:,:), occ_mo(:,:), vir_mo(:,:)
 character(len=240) :: nto_fch
 character(len=240), intent(in) :: fchname, logname
!f2py intent(in) :: fchname, logname

 call read_noa_nob_from_gau_log(logname, nfc, nif_ex, nocc, nob)
 ndb = nfc + nocc
 nvir = nif_ex - nocc
 allocate(exc(nocc,nvir,2))
 call read_ex_coeff_from_gau_log(logname, istate1, nocc, nvir, exc(:,:,1))
 call read_ex_coeff_from_gau_log(logname, istate2, nocc, nvir, exc(:,:,2))
 !exc = DSQRT(2d0)*exc  ! sum of square of exc(:,:,i) is 0.5, so multiply sqrt(2)

 allocate(tdm_occ(nocc,nocc), tdm_vir(nvir,nvir))
 call calc_cis_mo_tdm_using_exc(nocc, nvir, exc, tdm_occ, tdm_vir)
 deallocate(exc)

 allocate(u(nocc,nocc), vt(nocc,nocc), s(nocc))
 call do_svd(nocc, nocc, tdm_occ, u, vt, s)
 deallocate(tdm_occ, vt)
 do i = 1, nocc, 1
  s(i) = s(i)*s(i)
 end do ! for i

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(ev(nif), source=0d0)
 allocate(vir_mo(nocc,nocc))

 do i = 1, nocc, 1
  ev(nfc+i) = s(nocc-i+1)
  ev(ndb+i) = s(i)
  vir_mo(:,i) = u(:,nocc-i+1)
 end do ! for i

 deallocate(s)
 u = vir_mo
 deallocate(vir_mo)

 allocate(mo(nbf,nif))
 call read_mo_from_fch(fchname, nbf, nif, 'a', mo)
 allocate(occ_mo(nbf,nocc), source=0d0)
 ! C'_occ = (C_occ)U
 call dgemm('N', 'N', nbf, nocc, nocc, 1d0, mo(:,nfc+1:ndb), nbf, u, nocc, &
            0d0, occ_mo, nbf)
 mo(:,nfc+1:ndb) = occ_mo
 deallocate(occ_mo, u)

 allocate(u(nvir,nvir), vt(nvir,nvir), s(nvir))
 call do_svd(nvir, nvir, tdm_vir, u, vt, s)
 deallocate(tdm_vir, u, s)
 allocate(vir_mo(nbf,nvir), source=0d0)
 ! C'_vir = (C_vir)V
 call dgemm('N', 'T', nbf, nvir, nvir, 1d0, mo(:,ndb+1:ndb+nvir), nbf, vt, &
            nvir, 0d0, vir_mo, nbf)
 mo(:,ndb+1:ndb+nvir) = vir_mo
 deallocate(vir_mo, vt)

 call find_specified_suffix(fchname, '.fch', i)
 nto_fch = fchname(1:i-1)//'_NTO.fch'
 call copy_file(fchname, nto_fch, .false.)

 call write_mo_into_fch(nto_fch, nbf, nif, 'a', mo)
 deallocate(mo)
 call write_eigenvalues_to_fch(nto_fch, nif, 'a', ev, .true.)
 deallocate(ev)
end subroutine gen_cis_ee_nto_from_fch_and_log

subroutine del_a_or_b(str)
 implicit none
 integer :: i
 character(len=6), intent(inout) :: str

 str = ADJUSTL(str)
 i = LEN_TRIM(str)
 if(str(i:i)=='A' .or. str(i:i)=='B') str(i:i) = ' '
end subroutine del_a_or_b

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
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i < 0) exit ! end-of-file
  if(i > 0) then
   write(6,'(/,A)') 'ERROR in subroutine check_uhf_in_fch: failed to read file&
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
 integer, intent(in) :: nbf, nmo
!f2py intent(in) :: nbf, nmo
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
 real(kind=8), allocatable :: ttd(:,:)

 ev = 0d0; part_mo = 0d0; hole_mo = 0d0
 call calc_fosc_using_mo_tdm_ao_dip(nbf, nmo, delta_e, mo, mo_tdm, ao_dip, fosc)

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

