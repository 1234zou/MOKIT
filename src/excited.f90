! written by jxzou at 20210128: subroutines involving excited states wave function

! Currently only CIS/TDHF/TDDFT can be analyzed

subroutine ao2mo_tdm(nbf, nif, mo, nae, nav, exc, tdm, total)
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
end subroutine ao2mo_tdm

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
 open(newunit=fid,file=TRIM(logname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:6) == 'Range') exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine gen_tdm1: no 'Range' found in&
                   & file "//TRIM(logname)
  close(fid)
  stop
 end if

 read(fid,'(A)') buf
 i = index(buf,'NAE=')
 read(buf(i+4:),*) nae
 nav = nif - nae

 i = index(buf,'NBE=')
 read(buf(i+4:),*) nbe
 nbv = nif - nbe

 write(str,'(A13,I4)') 'Excited State', istate
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:18) == str) exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine gen_tdm1: no '"//str//"' found in&
                   & file "//TRIM(logname)
  close(fid)
  stop
 end if

 allocate(exc(nae,nav), source=0d0)

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

 call ao2mo_tdm(nbf, nif, mo, nae, nav, exc, tdm, .true.)
 deallocate(exc)
 return
end subroutine gen_ao_tdm

