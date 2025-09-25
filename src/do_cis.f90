! written by jxzou at 20220517: perform CIS/TDHF after HF calculations
! TODO: comparisons between various initial guess orbitals for SA-CAS
! TODO: develop spin-adapted RO-CIS methods in MOKIT, and replace MRSF-CIS by
!       spin-adapted RO-CIS methods of MOKIT.

! perform CIS/TDHF
subroutine do_cis()
 use phys_cons, only: au2ev
 use mol, only: mult, nbf, nif, nopen, nacto, nacta, nactb, nacte, uhf_ssquare
 use mr_keyword, only: mem, nproc, ist, eist, excited, sa_nto, tdhf, readrhf, &
  readuhf, nstate, dkh2_or_x2c, mixed_spin, casscf, cis_prog, gjfname, hf_fch,&
  gau_path
 use util_wrapper, only: unfchk, formchk, dat2fch_wrap
 implicit none
 integer :: i, na, nb, RENAME
 real(kind=8), parameter :: ssquare_thres = 0.5d0
 real(kind=8), allocatable :: cis_e(:), cis_ssquare(:), cis_fosc(:), ex_e(:)
 character(len=10) :: tmp_hf_prog
 character(len=240) :: nto_fch, cisno_fch, cis_gjf, cis_chk, cis_fch, cis_log,&
  uno_fch, pyname, outname, rohf_fch, datname, gmsname

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

 if(mult /= 1) then
  write(6,'(/,A)') 'ERROR in subroutine do_cis: only singlet is supported curre&
                   &ntly.'
  write(6,'(A,I0)') 'mult=', mult
  stop
 end if

 call read_nbf_and_nif_from_fch(hf_fch, nbf, nif)
 call read_na_and_nb_from_fch(hf_fch, na, nb)
 nopen = na - nb

 if(readrhf) then
  write(6,'(A)') 'Use R(O)HF .fch file provided by the user.'
  uhf_ssquare = 0d0
 else if(readuhf) then
  write(6,'(A)') 'Use UHF .fch file provided by the user.'
  call read_ssquare_from_fch(hf_fch, uhf_ssquare)
 end if
! TODO: for PySCF-generated UHF .fch file, we need to add <S^2> into it.

 allocate(cis_e(0:nstate), cis_ssquare(0:nstate), cis_fosc(nstate), ex_e(nstate))

 if(uhf_ssquare < ssquare_thres) then
  select case(TRIM(cis_prog))
  case('gaussian')
   call find_specified_suffix(gjfname, '.gjf', i)
   cis_gjf = gjfname(1:i-1)//'_CIS.gjf'
   cis_chk = gjfname(1:i-1)//'_CIS.chk'
   cis_fch = gjfname(1:i-1)//'_CIS.fch'
   cis_log = gjfname(1:i-1)//'_CIS.log'
   nto_fch = gjfname(1:i-1)//'_CIS_NTO.fch'
   cisno_fch = gjfname(1:i-1)//'_CIS_NO.fch'
   hf_fch = gjfname(1:i-1)//'_rhf.fch'
   call prt_cis_gjf(cis_gjf, mem, nproc, nstate, dkh2_or_x2c,mixed_spin,tdhf)
   call unfchk(hf_fch, cis_chk)
   call submit_gau_job(gau_path, cis_gjf, .true.)
   call read_cis_e_from_gau_log(cis_log, nstate, cis_e, cis_ssquare, cis_fosc)
   ex_e = (cis_e(1:nstate) - cis_e(0))*au2ev
   if(tdhf) then
    write(6,'(/,A)') 'TDHF energies from Gaussian:'
   else
    write(6,'(/,A)') 'CIS energies from Gaussian:'
   end if
   do i = 1, nstate, 1
    write(6,'(A,I4,A,F14.6,A,F6.2,A,F6.4,A,F8.3)') 'Excited State', i, ': E=',&
    cis_e(i),' a.u. (',ex_e(i),' eV), f=',cis_fosc(i),', <S**2>=',cis_ssquare(i)
   end do ! for i
   call formchk(cis_chk, cis_fch)
   if(sa_nto) then
    call gen_cis_ge_nto_from_fch_and_log(cis_fch, cis_log, nstate, .true.)
    hf_fch = nto_fch ! update hf_fch
   else
    call gen_cis_no_from_fch_and_log(cis_fch, cis_log, nstate, .true.)
    hf_fch = cisno_fch ! update hf_fch
   end if
   call delete_files(2, [cis_chk, cis_fch])
  case default
   write(6,'(/,A)') 'ERROR in subroutine do_cis: currently only CIS_prog=Gaussi&
                    &an is supported'
   write(6,'(A)') 'currently.'
   stop
  end select
 else
  ! perform MRSF-CIS calculations for an open-shell singlet species
  call find_specified_suffix(hf_fch, '.fch', i)
  uno_fch = hf_fch(1:i-1)//'_uno.fch'
  cis_fch = hf_fch(1:i-1)//'_t.fch'    ! temporarily use
  nto_fch = hf_fch(1:i-1)//'_t_NO.fch' ! temporarily use
  rohf_fch = hf_fch(1:i-1)//'_T_rohf.fch'
  datname = hf_fch(1:i-1)//'_T_rohf.dat'
  gmsname = hf_fch(1:i-1)//'_T_rohf.gms'
  pyname = hf_fch(1:i-1)//'.py'
  outname = hf_fch(1:i-1)//'.out'
  call find_specified_suffix(gjfname, '.gjf', i)
  cisno_fch = gjfname(1:i-1)//'_MRSFCIS_NO.fch'

  ! generate UNOs for singlet UHF
  call prt_uno_pyscf_script(hf_fch)
  call submit_pyscf_job(pyname, .true.)

  ! Perform a triplet ROHF calculation using singlet UNOs. Gaussian ROHF does
  ! not support stable/stable=opt, always use PySCF.
  tmp_hf_prog = 'pyscf'
  call do_rohf_using_uno(tmp_hf_prog, mem, nproc, dkh2_or_x2c, uno_fch)

  ! Perform MRSF-CIS based on the triplet ROHF
  call do_sf_cis_using_rohf(mem, nproc, nstate, .true., .false., rohf_fch)

  !call read_sf_e_from_gms_gms(gmsname, nstate, cis_e, cis_ssquare)
  call read_mrsf_e_from_gms_gms(gmsname, nstate, cis_e, cis_ssquare, cis_fosc)
  ex_e = (cis_e(1:nstate) - cis_e(0))*au2ev
  write(6,'(/,A)') 'MRSF-CIS energies from GAMESS:'
  write(6,'(A,F14.6,A,24X,A,F8.3)') 'Ground State    0: E=', cis_e(0), ' a.u.',&
                                    '<S**2>=', cis_ssquare(0)
  do i = 1, nstate, 1
   write(6,'(A,I4,A,F14.6,A,F6.2,A,F6.4,A,F8.3)') 'Excited State', i, ': E=',&
   cis_e(i),' a.u. (',ex_e(i),' eV), f=',cis_fosc(i),', <S**2>=',cis_ssquare(i)
  end do ! for i

  ! The GAMESS triplet ROHF MOs should be used to generate MRSF SA-NOs
  call sys_copy_file(TRIM(uno_fch), TRIM(cis_fch), .false.)
  call dat2fch_wrap(datname, cis_fch)
  call delete_files(3, [datname, pyname, outname])
  call gen_sfcis_no_from_fch_and_gms(cis_fch, gmsname, nstate, 0, nb-1, &
                                     .true., .true.)
  call delete_file(TRIM(cis_fch))
  i = RENAME(TRIM(nto_fch), TRIM(cisno_fch))
  hf_fch = cisno_fch ! update hf_fch
 end if

 deallocate(cis_e, cis_ssquare, cis_fosc, ex_e)
 casscf = .true. ! activate CASSCF for the subsequent SA-CAS calculation
 call get_active_space_for_sacas(hf_fch, sa_nto, nacto, nacta, nactb)
 nacte = nacta + nactb
end subroutine do_cis

! read expectation value <S^2> from a specified .fch(k) file
subroutine read_ssquare_from_fch(fchname, ssquare)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: ssquare
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname

 ssquare = 0d0
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:5) == 'S**2 ') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_ssquare_from_fch: no 'S**2 ' found&
                   & in file"
  write(6,'(A)') TRIM(fchname)
  stop
 end if

 read(buf(45:),*) ssquare
end subroutine read_ssquare_from_fch

! print Gaussian gjf file for CIS/TDHF
subroutine prt_cis_gjf(cis_gjf, mem, nproc, nstate, dkh2, mixed, tdhf)
 implicit none
 integer :: i, fid
 integer, intent(in) :: mem, nproc, nstate
 character(len=240) :: cis_chk
 character(len=240), intent(in) :: cis_gjf
 logical, intent(in) :: dkh2, mixed, tdhf

 call find_specified_suffix(cis_gjf, '.gjf', i)
 cis_chk = cis_gjf(1:i-1)//'.chk'
 open(newunit=fid,file=TRIM(cis_gjf),status='replace')

 write(fid,'(A)') '%chk='//TRIM(cis_chk)
 write(fid,'(A,I0,A)') '%mem=',mem,'GB'
 write(fid,'(A,I0)') '%nprocshared=',nproc
 write(fid,'(A)',advance='no') '#p RHF'
 if(tdhf) then
  write(fid,'(A)',advance='no') ' TD'
 else
  write(fid,'(A)',advance='no') ' CIS'
 end if
 write(fid,'(A)',advance='no') '(nstates='
 ! we calculate two more roots than nstate
 if(mixed) then
  write(fid,'(I0,A)',advance='no') (nstate+3)/2,',50-50'
 else
  write(fid,'(I0)',advance='no') nstate+2
 end if

 write(fid,'(A)',advance='no') ') chkbasis nosymm int'
 if(dkh2) then
  write(fid,'(A)',advance='no') '(DKH2,nobasistransform)'
 else
  write(fid,'(A)',advance='no') '=nobasistransform'
 end if
 write(fid,'(A,/)') ' guess=read geom=allcheck iop(9/40=4)'
 close(fid)
end subroutine prt_cis_gjf

subroutine read_cis_e_from_gau_log(logname, nstate, cis_e, ssquare, fosc)
 use phys_cons, only: au2ev
 implicit none
 integer :: i, j, k, fid
 integer, intent(in) :: nstate
 real(kind=8), intent(out) :: cis_e(0:nstate), ssquare(0:nstate), fosc(nstate)
 character(len=17) :: key
 character(len=240) :: buf
 character(len=240), intent(in) :: logname
 logical, allocatable :: found(:)

 cis_e = 0d0; ssquare = 0d0; fosc = 0d0
 open(newunit=fid,file=TRIM(logname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:9) == 'SCF Done') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cis_e_from_gau_log: no 'SCF&
                   & Done' found"
  write(6,'(A)') 'in file '//TRIM(logname)
  close(fid)
  stop
 end if

 i = INDEX(buf, '=')
 read(buf(i+1:),*) cis_e(0) ! HF energy in a.u.

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:19) == 'Excited state symm') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_cis_e_from_gau_log: no 'Exc&
                   &ited state symm'"
  write(6,'(A)') 'found in file '//TRIM(logname)
  close(fid)
  stop
 end if

 write(key,'(A,I4)') 'Excited State', k+1
 allocate(found(nstate))
 found = .false.
 k = 0

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:7) == 'SavETr') exit
  if(buf(2:18) == key) then
   k = k + 1
   found(k) = .true.
   i = INDEX(buf, 'Sym'); j = INDEX(buf, 'eV')
   read(buf(i+3:j-1),*) cis_e(k) ! in eV
   i = INDEX(buf, 'f=')
   read(buf(i+2:),*) fosc(k)
   i = INDEX(buf, '<S**2>=')
   read(buf(i+7:),*) ssquare(k)
   if(k == nstate) exit
   write(key,'(A,I4)') 'Excited State', k+1
   do while(.true.)
    read(fid,'(A)') buf
    if(LEN_TRIM(buf) == 0) exit
   end do ! for while
  end if
 end do ! for while

 close(fid)
 if(ANY(found .eqv. .false.)) then
  write(6,'(/,A)') 'ERROR in subroutine read_cis_e_from_gau_log: some ex&
                   &cited state is'
  write(6,'(A)') 'not located in file '//TRIM(logname)
  write(6,'(40L2)') found
  deallocate(found)
  stop
 end if

 deallocate(found)
 ! convert excitation energies to (absolute) electronic energies
 cis_e(1:nstate) = cis_e(1:nstate)/au2ev + cis_e(0)
end subroutine read_cis_e_from_gau_log

! print/create a PySCF script for generating UNOs
subroutine prt_uno_pyscf_script(uhf_fch)
 implicit none
 integer :: i, fid
 character(len=240) :: pyname, uno_fch0, uno_fch
 character(len=240), intent(in) :: uhf_fch

 call find_specified_suffix(uhf_fch, '.fch', i)
 uno_fch0 = uhf_fch(1:i-1)//'_UNO.fch'
 uno_fch = uhf_fch(1:i-1)//'_uno.fch'
 pyname = uhf_fch(1:i-1)//'.py'

 open(newunit=fid,file=TRIM(pyname),status='replace')
 write(fid,'(A)') 'from mokit.lib.gaussian import uno'
 write(fid,'(A)') 'from os import rename'
 write(fid,'(A)') "uhf_fch = '"//TRIM(uhf_fch)//"'"
 write(fid,'(A)') "uno_fch0 = '"//TRIM(uno_fch0)//"'"
 write(fid,'(A)') "uno_fch = '"//TRIM(uno_fch)//"'"
 write(fid,'(A)') 'uno(uhf_fch)'
 write(fid,'(A)') 'rename(uno_fch0, uno_fch)'
 close(fid)
end subroutine prt_uno_pyscf_script

! modify/print a PySCF ROHF script
subroutine prt_rohf_pyscf_script(mem, nproc, pyname)
 implicit none
 integer :: i, fid, fid1, RENAME
 integer, intent(in) :: mem, nproc
 character(len=240) :: buf, pyname1
 character(len=240), intent(in) :: pyname

 call find_specified_suffix(pyname, '.py', i)
 pyname1 = pyname(1:i-1)//'.t'
 open(newunit=fid,file=TRIM(pyname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(pyname1),status='replace')

 do while(.true.)
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while
 write(fid1,'(A)') 'from mokit.lib.py2fch import py2fch'
 write(fid1,'(A,/)') 'from mokit.lib.stability import rohf_stable_opt_internal'

 read(fid,'(A)') buf
 if(buf(1:15) == 'lib.num_threads') then
  write(fid1,'(A,I0,A)') 'lib.num_threads(',nproc,')'
 else
  write(6,'(/,A)') "ERROR in subroutine prt_rohf_pyscf_script: 'lib.num_threads&
                   &' is not located"
  write(6,'(A)') 'at the target line.'
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:13) == 'mf.max_memory') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while
 write(fid1,'(A,I0,A)') 'mf.max_memory= ', mem*1000, ' # MB'

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:18) == '#dm = mf.make_rdm1') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 write(fid1,'(A)') 'dm = mf.make_rdm1()'
 write(fid1,'(A)') 'mf.max_cycle = 128'
 write(fid1,'(A)') 'mf.kernel(dm0=dm)'
 write(fid1,'(A)') 'if mf.converged is False:'
 write(fid1,'(A)') '  mf = mf.newton()'
 write(fid1,'(A)') '  mf.kernel()'
 write(fid1,'(A)') 'mf = rohf_stable_opt_internal(mf)'
 write(fid1,'(A)') "py2fch(hf_fch, nbf, nif, mf.mo_coeff, 'a', mf.mo_energy, Fa&
                   &lse, True)"
 close(fid1)
 i = RENAME(TRIM(pyname1), TRIM(pyname))
end subroutine prt_rohf_pyscf_script

! Perform high spin ROHF calculation at spin S+1 using UNO at spin S.
! For example, triplet ROHF using UNO of open-shell singlet
subroutine do_rohf_using_uno(hf_prog, mem, nproc, dkh2, uno_fch)
 use mr_keyword, only: gau_path
 use util_wrapper, only: bas_fch2py_wrap, unfchk, formchk
 implicit none
 integer :: i, charge, mult, fid
 integer, intent(in) :: mem, nproc
 character(len=10), intent(in) :: hf_prog
 character(len=240) :: rohf_inp, rohf_chk, rohf_fch
 character(len=240), intent(in) :: uno_fch
 logical, intent(in) :: dkh2

 call read_charge_and_mult_from_fch(uno_fch, charge, mult)
 call find_specified_suffix(uno_fch, '_uno', i)

 select case(TRIM(hf_prog))
 case('pyscf')
  rohf_inp = uno_fch(1:i-1)//'_T_rohf.py'
  rohf_fch = uno_fch(1:i-1)//'_T_rohf.fch'
  call copy_file(uno_fch, rohf_fch, .false.)
  call modify_charge_and_mult_in_fch(rohf_fch, charge, mult+2)
  call bas_fch2py_wrap(rohf_fch, .false.)
  call prt_rohf_pyscf_script(mem, nproc, rohf_inp)
  call submit_pyscf_job(rohf_inp, .true.)
 case('gaussian')
  rohf_inp = uno_fch(1:i-1)//'_T_rohf.gjf'
  rohf_chk = uno_fch(1:i-1)//'_T_rohf.chk'
  rohf_fch = uno_fch(1:i-1)//'_T_rohf.fch'
  open(newunit=fid,file=TRIM(rohf_inp),status='replace')
  write(fid,'(A)') '%chk='//TRIM(rohf_chk)
  write(fid,'(A,I0,A)') '%mem=', mem, 'GB'
  write(fid,'(A,I0)') '%nprocshared=', nproc
  write(fid,'(A)',advance='no') '#p ROHF chkbasis nosymm'
  if(dkh2) then
   write(fid,'(A)',advance='no') ' int(nobasistransform,DKH2)'
  else
   write(fid,'(A)',advance='no') ' int=nobasistransform'
  end if
  write(fid,'(A)') ' scf(maxcycle=512) guess=read geom=check'
  write(fid,'(/,A,/)') 'ROHF file generated by MOKIT'
  write(fid,'(I0,1X,I0,/)') charge, mult+2
  close(fid)
  call unfchk(uno_fch, rohf_chk)
  call submit_gau_job(gau_path, rohf_inp, .true.)
  call formchk(rohf_chk)
  call delete_file(TRIM(rohf_chk))
 case default
  write(6,'(/,A)') 'ERROR in subroutine do_rohf_using_uno: invalid HF_prog='//&
                    TRIM(hf_prog)
  write(6,'(A)') 'Only PySCF/Gaussian are supported in this step.'
  stop
 end select
end subroutine do_rohf_using_uno

! Perform SF-/MRSF-CIS calculation using high-spin ROHF wfn.
! nstate: the number of excited states. E.g. S0/S1/S2 means nstate=2.
! SF-/MRSF-CIS will compute nstate+3 roots.
subroutine do_sf_cis_using_rohf(mem, nproc, nstate, mrsf, triplet, rohf_fch)
 use mr_keyword, only: check_gms_path, gms_path, gms_scr_path, gms_dat_path
 implicit none
 integer :: i, mult, fid, fid1, SYSTEM, RENAME
 integer, intent(in) :: mem, nproc, nstate
 character(len=240) :: buf, tmpf, inpname, gmsname
 character(len=240), intent(in) :: rohf_fch
 logical, intent(in) :: mrsf, triplet
 ! currently only the triplet is allowed for MRSF

 mult = 1
 if(mrsf) then
  i = SYSTEM('fch2inp '//TRIM(rohf_fch)//' -mrsfcis')
  if(triplet) mult = 3
 else
  i = SYSTEM('fch2inp '//TRIM(rohf_fch)//' -sfcis')
 end if

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine do_sf_cis_using_rohf: failed to call ut&
                   &ility fch2inp.'
  write(6,'(A)') 'rohf_fch='//TRIM(rohf_fch)
  stop
 end if

 call find_specified_suffix(rohf_fch, '.fch', i)
 inpname = rohf_fch(1:i-1)//'.inp'
 gmsname = rohf_fch(1:i-1)//'.gms'
 tmpf = rohf_fch(1:i-1)//'.t'

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(tmpf),status='replace')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:8) == '$SYSTEM') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 i = FLOOR(125d0*DBLE(mem)/DBLE(nproc))
 write(fid1,'(1X,A,I0,A)') '$SYSTEM MWORDS=', i, ' $END'

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:6)=='$CIS ' .or. buf(2:7)=='$TDDFT') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(mrsf) then
  write(fid1,'(2(A,I0),A)') ' $TDDFT NSTATE=',nstate+4,' MULT=',mult,' $END'
  ! MRSF-CIS is used here, so NRAD and NLEB are not needed
 else
  write(fid1,'(A,I0,A)') ' $CIS NSTATE=', nstate+4, ' $END'
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(tmpf), TRIM(inpname))

 call check_gms_path()
 call submit_gms_job(gms_path, gms_scr_path, gms_dat_path, inpname, nproc)
end subroutine do_sf_cis_using_rohf

! Read electronic energies of SF methods from a specified GAMESS .gms file.
! SF-CIS/TD in GAMESS has no oscillator strength, so there is no fosc below.
subroutine read_sf_e_from_gms_gms(gmsname, nstate, e, ssquare)
 implicit none
 integer :: i, j, k, fid
 integer, intent(in) :: nstate
 real(kind=8) :: rtmp
 real(kind=8), intent(out) :: e(0:nstate), ssquare(0:nstate)
 character(len=2) :: str
 character(len=27), parameter :: key = 'SUMMARY OF SPIN-FLIP RESULT'
 character(len=240) :: buf
 character(len=240), intent(in) :: gmsname
 logical, allocatable :: found(:)

 e = 0d0; ssquare = 0d0
 open(newunit=fid,file=TRIM(gmsname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(11:37) == key) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_sf_e_from_gms_gms: '"//key//"'"
  write(6,'(A)') 'is not found in file '//TRIM(gmsname)
  close(fid)
  stop
 end if

 do i = 1, 5   ! skip 5 lines
  read(fid,'(A)') buf
 end do ! for i
 allocate(found(0:nstate))
 found = .false.

 do i = 1, nstate+2, 1
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
  read(buf,*) k
  if(k>0 .and. k<=nstate+1) then
   found(k-1) = .true.
   read(buf,*) j, str, e(k-1), rtmp, ssquare(k-1)
  end if
 end do ! for k

 close(fid)
 if(ANY(found .eqv. .false.)) then
  write(6,'(/,A)') 'ERROR in subroutine read_sf_e_from_gms_gms: some state is n&
                   &ot found in file'
  write(6,'(A)') TRIM(gmsname)
  write(6,'(40L2)') found
  deallocate(found)
  stop
 end if

 deallocate(found)
end subroutine read_sf_e_from_gms_gms

! read electronic energies of MRSF methods from a specified GAMESS .gms gile
subroutine read_mrsf_e_from_gms_gms(gmsname, nstate, e, ssquare, fosc)
 implicit none
 integer :: i, j, k, fid
 integer, intent(in) :: nstate
 real(kind=8) :: rtmp
 real(kind=8), intent(out) :: e(0:nstate), ssquare(0:nstate), fosc(nstate)
 character(len=2) :: str
 character(len=27), parameter :: key = 'SUMMARY OF MRSF-DFT RESULTS'
 character(len=240) :: buf
 character(len=240), intent(in) :: gmsname
 logical, allocatable :: found(:)

 e = 0d0; ssquare = 0d0; fosc = 0d0
 open(newunit=fid,file=TRIM(gmsname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(buf(2:14) == 'SCF IS UNCONV')then
   write(6,'(/,A)') 'ERROR in subroutine read_mrsf_e_from_gms_gms: ROHF SCF con&
                    &vergence failure.'
   write(6,'(A)') 'Please open and check file '//TRIM(gmsname)
   close(fid)
   stop
  end if
  if(buf(2:8) == 'ITER EX') then
   write(6,'(/,A)') "ERROR in subroutine read_mrsf_e_from_gms_gms: '"//key//"'"
   write(6,'(A)') 'not found in file '//TRIM(gmsname)
   close(fid)
   stop
  end if
  if(buf(11:37) == key) exit
 end do ! for while

 do i = 1, 6   ! skip 6 lines
  read(fid,'(A)') buf
 end do ! for i
 allocate(found(0:nstate))
 found = .false.

 do i = 1, nstate+2, 1
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
  read(buf,*) k
  if(k>0 .and. k<=nstate+1) then
   found(k-1) = .true.
   read(buf,*) j, str, e(k-1), rtmp, ssquare(k-1)
  end if
 end do ! for k

 if(ANY(found .eqv. .false.)) then
  write(6,'(/,A)') 'ERROR in subroutine read_mrsf_e_from_gms_gms: some state is&
                   & not found in file'
  write(6,'(A)') TRIM(gmsname)
  write(6,'(40L2)') found
  deallocate(found)
  close(fid)
  stop
 end if
 deallocate(found)

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:11) == 'TRANSITION') exit
 end do ! for while
 read(fid,'(A)') buf
 read(fid,'(A)') buf

 do i = 1, nstate, 1
  read(fid,'(A)') buf
  k = LEN_TRIM(buf)
  j = INDEX(buf(1:k), ' ', back=.true.)
  read(buf(j+1:k),*) fosc(i)
 end do ! for i

 close(fid)
end subroutine read_mrsf_e_from_gms_gms

! determine the size of active space for the subsequent SA-CASSCF calculation
subroutine get_active_space_for_sacas(fchname, sa_nto, nacto, nacta, nactb)
 implicit none
 integer :: i, nbf, nif, na, nb
 integer, intent(out) :: nacto, nacta, nactb
 real(kind=8), parameter :: nto_thres = 0.01d0
 real(kind=8), parameter :: no_thres = 0.02d0
 real(kind=8), allocatable :: ev(:)
 character(len=240), intent(in) :: fchname
 logical, intent(in) :: sa_nto

 nacto = 0; nacta = 0; nactb = 0
 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(ev(nif))
 call read_eigenvalues_from_fch(fchname, nif, 'a', ev)
 call read_na_and_nb_from_fch(fchname, na, nb)

 if(sa_nto) then
  ! NTO eigenvales always occur in pairwise
  nacto = COUNT(ev(1:na) > nto_thres)
  deallocate(ev)
  nacta = nacto
  nactb = nacto
  nacto = 2*nacto
 else
  do i = 1, nif, 1
   if(ev(i)>=no_thres .and. ev(i)<=2d0-no_thres) then
    nacto = nacto + 1
    if(i <= na) nacta = nacta + 1
    if(i <= nb) nactb = nactb + 1
   end if
  end do ! for i
 end if
end subroutine get_active_space_for_sacas

