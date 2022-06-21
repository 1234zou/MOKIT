! written by jxzou at 20220517: perform CIS/TDHF after HF calculations

! perform CIS/TDHF
subroutine do_cis()
 use mol, only: nacto, nacta, nactb, nacte
 use mr_keyword, only: ist, eist, excited, cis_prog, gjfname, nstate, nproc, &
  mem, hf_fch, casnofch, gau_path, casscf, TDHF
 use util_wrapper, only: unfchk, formchk
 implicit none
 integer :: i, system
 character(len=240) :: ntofch, cisgjf, cischk, cisfch, cislog

 if(.not. excited) return
 if(eist == 0) then
  return
 else
  if(eist /= 1) then
   write(6,'(/,A)') 'ERROR in subroutine do_cis: eist/=1.'
   stop
  end if
  if(ist == 5) return ! readno
 end if

 select case(cis_prog)
 case('gaussian')
  i = index(gjfname, '.gjf', back=.true.)
  cisgjf = gjfname(1:i-1)//'_CIS.gjf'
  cischk = gjfname(1:i-1)//'_CIS.chk'
  cisfch = gjfname(1:i-1)//'_CIS.fch'
  cislog = gjfname(1:i-1)//'_CIS.log'
  ntofch = gjfname(1:i-1)//'_CIS_NTO.fch'
  hf_fch = gjfname(1:i-1)//'_rhf.fch'

  call prt_cis_gjf(cisgjf, nstate, nproc, mem, TDHF)
  call unfchk(hf_fch, cischk)

  write(6,'(A)') '$'//TRIM(gau_path)//' '//TRIM(cisgjf)
  i = system(TRIM(gau_path)//' '//TRIM(cisgjf))
  if(i /= 0) then
   write(6,'(/,A)') 'ERROR in subroutine do_cis: Gaussian CIS job failed.'
   write(6,'(A)') 'Please open file '//TRIM(cisgjf)//' and check.'
   stop
  end if
  call read_cis_energies_from_gau_log(cislog)
  call formchk(cischk, cisfch)
  call delete_file(cischk)
  hf_fch = ntofch ! update hf_fch
  call copy_file(cisfch, hf_fch, .false.)
  call gen_nto(cislog, hf_fch, nstate, .true.)

 case default
  write(6,'(A)') 'ERROR in subroutine do_cis: currently only CIS_prog=Gaussian&
                & is supported currently.'
  stop
 end select

 casscf = .true. ! activate CASSCF for SA-CASSCF calculation
 call get_active_space_for_SACAS(hf_fch, nacto, nacta, nactb)
 nacte = nacta + nactb
end subroutine do_cis

! print Gaussian gjf file for CIS/TDHF
subroutine prt_cis_gjf(cisgjf, nstate, nproc, mem, TDHF)
 implicit none
 integer :: i, fid
 integer, intent(in) :: nstate, nproc, mem
 character(len=240) :: cischk
 character(len=240), intent(in) :: cisgjf
 logical, intent(in) :: TDHF

 i = index(cisgjf, '.gjf', back=.true.)
 cischk = cisgjf(1:i-1)//'.chk'

 open(newunit=fid,file=TRIM(cisgjf),status='replace')
 write(fid,'(A)') '%chk='//TRIM(cischk)
 write(fid,'(A,I0,A)') '%mem=',mem,'GB'
 write(fid,'(A,I0)') '%nprocshared=',nproc
 ! we calculate two more roots than desired
 if(TDHF) then
  write(fid,'(A,I0,A,/)') '#p HF TD(nstates=',nstate+2,') chkbasis nosymm int=n&
                          &obasistransform guess=read geom=allcheck iop(9/40=5)'
 else
  write(fid,'(A,I0,A,/)') '#p CIS(nstates=',nstate+2,') chkbasis nosymm int=nob&
                          &asistransform guess=read geom=allcheck iop(9/40=5)'
 end if
 close(fid)
end subroutine prt_cis_gjf

subroutine read_cis_energies_from_gau_log(logname)
 implicit none
 integer :: i, mult, fid
 real(kind=8) :: e
 character(len=240) :: buf
 character(len=240), intent(in) :: logname

 open(newunit=fid,file=TRIM(logname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:19) == 'Excited state symm') exit
 end do ! for while

 if(i /= 0) then
  close(fid)
  write(6,'(A)') "ERROR in subroutine read_cis_energies_from_gau_log: no 'Excit&
                 &ed state symm' found in file "//TRIM(logname)
  stop
 end if

 write(6,'(/,A)') 'CIS energies from Gaussian log:'
 BACKSPACE(fid)
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:7) == 'SavETr') exit
  if(buf(2:19) == 'Excited state symm') then
   read(fid,'(A)') buf
   write(6,'(A)') TRIM(buf)
  end if
 end do ! for while

 close(fid)
end subroutine read_cis_energies_from_gau_log

! determine the size of active space for subsequent SA-CASSCF calculation
subroutine get_active_space_for_SACAS(fchname, nacto, nacta, nactb)
 implicit none
 integer :: nbf, nif, na, nb
 integer, intent(out) :: nacto, nacta, nactb
 real(kind=8), parameter :: nto_thres = 1d-2
 real(kind=8), allocatable :: ev(:)
 character(len=240), intent(in) :: fchname
 ! currently assume the fchname contains state-avaraged NTO

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(ev(nif))
 call read_eigenvalues_from_fch(fchname, nif, 'a', ev)

 call read_na_and_nb_from_fch(fchname, na, nb)
 nacto = COUNT(ev(1:na) > nto_thres)
 deallocate(ev)

 nacta = nacto
 nactb = nacto
 nacto = 2*nacto
end subroutine get_active_space_for_SACAS

