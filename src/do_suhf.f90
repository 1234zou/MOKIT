! written by jeanwsr at Jan 2025: subroutines for SUHF

subroutine do_suhf()
 use mr_keyword, only: mo_rhf, ist, hf_fch, bgchg, chgname, nskip_uno, &
  loc_asrot, suhf_prog
 use mol, only: chem_core, ecp_core
 use util_wrapper, only: bas_fch2py_wrap
 implicit none
 integer :: i, hf_type
 real(kind=8) :: unpaired_e
 real(kind=8) :: e, ssquare
 character(len=24) :: data_string = ' '
 character(len=240) :: proname, pyname, outname, fchname, fchname_suno

 if(ist /= 7) return  ! no need for this subroutine
 write(6,'(//,A)') 'Enter subroutine do_suhf...'

 ! calculate the number of core orbitals from array core_orb
 call calc_ncore(hf_fch, chem_core, ecp_core)
 write(6,'(3(A,I0))') 'chem_core=', chem_core, ', ecp_core=', ecp_core, &
                      ', Nskip_UNO=', nskip_uno

 i = INDEX(hf_fch, '.fch', back=.true.)
 proname = hf_fch(1:i-1)
 write(6,'(A)') 'SUHF using program '//TRIM(suhf_prog)

 if(mo_rhf) then ! currently impossible for suhf
!   if(ist == 6) then
!    write(6,'(A)') 'One set of MOs: GVB/STO-6G -> MO projection -> Rotate has be&
!                   &en invoked.'
!   else
!    write(6,'(A)') 'One set of MOs: invoke RHF virtual MO projection -> localiza&
!                   &tion -> paring.'
!    pyname = TRIM(proname)//'_proj_loc_pair.py'
!    call prt_rhf_proj_script_into_py(pyname)
!    call prt_auto_pair_script_into_py(pyname)
!    call submit_pyscf_job(pyname, .true.)
!   end if

 else

  fchname = hf_fch(1:i-1)//'_suhf.fch'
  fchname_suno = hf_fch(1:i-1)//'_suno.fch'

  if(loc_asrot) then
    write(6,'(A)')' invoke SUNO associated rotation.'
    pyname = TRIM(proname)//'_suno_asrot.py'
    outname = TRIM(proname)//'_suno_asrot.out'
  else
    pyname = TRIM(proname)//'_suhf.py'
    outname = TRIM(proname)//'_suhf.out'
  end if
  call bas_fch2py_wrap(hf_fch, .false., pyname)
  call prt_suhf_script_into_py(pyname)
  if(loc_asrot) call prt_assoc_rot_script_into_py(pyname, .true.)
  
  !if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(pyname))
  call submit_pyscf_job(pyname, .true.)
  hf_type = 3 ! although there's no other cases, let's keep 
            !this type for future change.
  call read_suhf_e_and_ss_from_exscf_out(outname, hf_type, e, ssquare)
  ! ssquare is the deformed spin square, not the pure one
  ! currently not necessary to print it
  write(6,'(/,A,F18.8,1X,A)') 'E(SUHF) = ', e,'a.u.'
  call calc_unpaired_from_fch(fchname_suno, 1, .false., unpaired_e)

 end if

 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine do_suhf at '//TRIM(data_string)

end subroutine do_suhf

subroutine prt_suhf_script_into_py(pyname)
 use mr_keyword, only: nproc, hf_fch, uno_thres
 implicit none
 integer :: i, fid1, fid2, RENAME
 character(len=240) :: buf, pyname1, suhf_fch, suno_fch, suno_out
 character(len=240), intent(in) :: pyname

 buf = ' '
 call find_specified_suffix(pyname, '.py', i)
 pyname1 = pyname(1:i-1)//'.t'

 call find_specified_suffix(hf_fch, '.fch', i)
 suhf_fch = hf_fch(1:i-1)//'_suhf.fch'
 suno_fch = hf_fch(1:i-1)//'_suno.fch'
 suno_out = hf_fch(1:i-1)//'_suno.txt'

 open(newunit=fid1,file=TRIM(pyname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(pyname1),status='replace')

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine prt_suhf_script_into_py: end-of-file det&
                   &ected.'
  write(6,'(A)') 'File may be problematic: '//TRIM(pyname)
  close(fid1)
  close(fid2,status='delete')
  stop
 end if

 write(fid2,'(A)') 'import numpy as np'
 write(fid2,'(A)') 'import os'
 write(fid2,'(A)') 'from mokit.lib.py2fch import py2fch'
 !write(fid2,'(A)') 'from mokit.lib.uno import uno'
 write(fid2,'(A)') 'from mokit.lib.gaussian import mo_fch2py'
 write(fid2,'(A)') 'from mokit.lib.rwwfn import read_nbf_and_nif_from_fch'
 !                  '&, read_na_and_nb_from_fch, \'
 !write(fid2,'(A)') ' construct_vir'
 write(fid2,'(A)') 'from mokit.lib.util import get_idx_from_noon'
 write(fid2,'(A)') 'from pyphf import suscf'
 write(fid2,'(A,I0)') 'nproc = ', nproc
 write(fid2,'(A,/)') 'lib.num_threads(nproc)'

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(1:15) == 'lib.num_threads') cycle
  if(buf(1:9) == 'mf.kernel') exit
  write(fid2,'(A)') TRIM(buf)
 end do

 write(fid2,'(/,A)') '# read MOs from .fch(k) file'
 write(fid2,'(A)') 'mf.kernel()' ! run 1 cycle because suhf needs mo_occ
 write(fid2,'(A)') "hf_fch = '"//TRIM(hf_fch)//"'"
 write(fid2,'(A)') 'mf.mo_coeff = mo_fch2py(hf_fch)'
 write(fid2,'(A)') 'nbf, nif = read_nbf_and_nif_from_fch(hf_fch)'
 write(fid2,'(A)') '# read done'
 write(fid2,'(/,A)') '# check if input MOs are orthonormal'
 write(fid2,'(A)') "S = mol.intor_symmetric('int1e_ovlp')"
 write(fid2,'(A)') 'check_orthonormal(nbf, nif, mf.mo_coeff[0], S)'
 write(fid2,'(A)') 'check_orthonormal(nbf, nif, mf.mo_coeff[1], S)'

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(TRIM(pyname1), TRIM(pyname))

 open(newunit=fid1,file=TRIM(pyname),status='old',position='append')
 write(fid1,'(/,A)') '# do SUHF'
 write(fid1,'(A)') 'smf = suscf.SUHF(mf)'
 write(fid1,'(A)') 'smf.kernel()'
 write(fid1,'(A)') 'sumo = smf.mo_reg'
 write(fid1,'(A)') 'suno = smf.natorb[2]'
 write(fid1,'(A)') 'sunoon = smf.natocc[2]'

 !write(fid1,'(/,A)') '# save the SUHF MO into .fch file'
 ! not used now, skip

 write(fid1,'(/,A)') '# save the SUNO into .fch file'
 write(fid1,'(A)') "suno_fch = '"//TRIM(suno_fch)//"'"
 write(fid1,'(A)') "suno_out = '"//TRIM(suno_out)//"'"
 write(fid1,'(A)') 'nopen = mol.spin'
 ! compute idx for next step asrot
 write(fid1,'(A,F12.10)') 'uno_thres = ', uno_thres
 write(fid1,'(A)') 'idx = get_idx_from_noon(suno_out, sunoon, nbf, nif, nopen, &
                   &uno_thres)'
 write(fid1,'(A)') 'na, nb = mf.nelec'
 write(fid1,'(A)') 'mf.mo_coeff = (suno, suno)'
!  write(fid1,'(/,A)') '# transform UHF canonical orbitals to UNO'
!  write(fid1,'(A)') 'mo_a = mf.mo_coeff[0].copy()'
!  write(fid1,'(A)') 'na, nb = read_na_and_nb_from_fch(hf_fch)'
!  write(fid1,'(A,E12.5,A)') 'idx, noon, alpha_coeff = uno(nbf,nif,na,nb,mf.mo_co&
!                            &eff[0], mf.mo_coeff[1],S,',uno_thres,')'
!  write(fid1,'(A)') 'alpha_coeff = construct_vir(nbf, nif, idx[1], alpha_coeff, S)'
!  write(fid1,'(A)') 'mf.mo_coeff = (alpha_coeff, alpha_coeff)'
!  write(fid1,'(A)') '# done transform'

 !write(fid1,'(A)') "os.system('fch_u2r "//TRIM(hf_fch)//' '//TRIM(uno_fch)//"')"
 !write(fid1,'(A)') 'sleep(1) # in some node, py2fch begins when fch_u2r unfinished'
 ! Thanks to the suggestion of Kalinite. The problem in the above two lines
 ! does not exist now.
 write(fid1,'(A)') "with os.popen('fch_u2r '+hf_fch+' '+suno_fch) as run:"
 write(fid1,'(A)') '  null = run.read()'
 write(fid1,'(A)') "py2fch(suno_fch,nbf,nif,suno,'a',sunoon,True,True)"
 write(fid1,'(A)') '# save done'
 close(fid1)
end subroutine prt_suhf_script_into_py

! read HF electronic energy from a PySCF .out file
subroutine read_suhf_e_and_ss_from_exscf_out(outname, wfn_type, e, ss)
 implicit none
 integer :: i, j, fid
 integer, intent(in) :: wfn_type
 real(kind=8), intent(out) :: e, ss
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 e = 0d0; ss = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')
 write(fid,'(/)',advance='no')
 ! add a blank line, in case 'Final E(SUHF)' is already in the last line

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(1:13) == 'Final E(SUHF)') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_suhf_e_and_ss_from_exscf_out: no &
                   &'Final E(SUHF)' found"
  write(6,'(A)') 'in file '//TRIM(outname)
  close(fid)
  stop
 end if

 j = INDEX(buf, '=')
 read(buf(j+1:),*) e
 !close(fid)

 select case(wfn_type)
!  case(1,2) ! R(O)HF
!   i = INDEX(outname, '.out', back=.true.)
!   inpname = outname(1:i-1)//'.py'
!   call read_mult_from_pyscf_inp(inpname, mult)
!   ss = DBLE((mult-1))*0.5d0
!   ss = DBLE(ss*(ss+1))
 case(3)   ! SUHF. 
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:8) == 'deformed') exit
  end do ! for while
  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_suhf_e_and_ss_from_exscf_out: no &
                   &'deformed' found"
   write(6,'(A)') 'in file '//TRIM(outname)
   close(fid)
   stop
  end if
  j = INDEX(buf, '=')
  read(buf(j+1:),*) ss
  close(fid)
 case default
  write(6,'(/,A,I0)') 'ERROR in subroutine read_suhf_e_and_ss_from_exscf_out: i&
                      &nvalid wfn_type=', wfn_type
  stop
 end select
end subroutine read_suhf_e_and_ss_from_exscf_out

