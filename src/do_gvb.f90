! written by jxzou at 20210729: move subroutine do_gvb from file automr.f90 to
! this new file, and rename it as do_gvb_gms
! updated by jxzou at 20210915: add excludeXH section

! perform GVB computation (only in Strategy 1,3) using GAMESS/QChem/Gaussian
subroutine do_gvb()
 use mr_keyword, only: gvb, gvb_prog, ist, hf_fch, mo_rhf, npair_wish, onlyXH,&
  excludeXH, localm, LocDocc
 use mol, only: nbf, nif, ndb, npair, nopen, lin_dep, nacta, nactb, nacte, &
  nacto, npair0
 use util_wrapper, only: gvb_exclude_XH_A_wrap
 implicit none
 integer :: i, nvir
 character(len=24) :: data_string = ' '
 character(len=240) :: proname, proname1, pair_fch, sort_fch, inpname, datname,&
                       gmsname
 logical :: pm_loc

 if(.not. gvb) return
 write(6,'(//,A)') 'Enter subroutine do_gvb...'

 select case(ist)
 case(1,3,6)
 case default
  write(6,'(/,A,I0)') 'ERROR in subroutine do_gvb: ist=', ist
  write(6,'(A)') 'Only ist=1,3,6 is supported currently.'
  stop
 end select

 i = INDEX(hf_fch, '.fch', back=.true.)
 proname = hf_fch(1:i-1)

 ! In RHF virtual MO projection, it will generate a file uno.out additionally
 call read_npair_from_uno_out(nbf, nif, ndb, npair, nopen, lin_dep)
 write(6,'(2(A,L1))') 'LocDocc=', LocDocc, ', Lin_dep=', lin_dep

 if(npair_wish>0 .and. npair_wish/=npair) then
  write(6,'(/,2(A,I0),A)') 'Warning: AutoMR recommends GVB(',npair,'), but user&
                           & specifies GVB(',npair_wish,'). Try to fulfill...'
  if(npair_wish < npair) then
   ndb = ndb + npair - npair_wish
   npair = npair_wish
   nvir = nif - ndb - 2*npair - nopen
   write(6,'(A)') 'OK, fulfilled. You are supposed to check GVB orbitals after &
                  &converged.'
   write(6,'(/,A)') 'Updated number of various orbitals:'
   write(6,'(4(A,I0))') 'doubly_occ=',ndb,', npair=',npair,', nopen=',nopen, &
                        ', nvir=', nvir
  else if(npair_wish > npair) then
   write(6,'(/,A)') 'ERROR in subroutine do_gvb: too many pairs specified. Cann&
                    &ot be fulfilled.'
   stop
  end if
 end if

 if(mo_rhf) then
  if(ist == 6) then
   pair_fch = hf_fch
  else ! ist /= 6
   pair_fch = TRIM(proname)//'_proj_loc_pair.fch'
  end if
 else ! not mo_rhf
  pair_fch = TRIM(proname)//'_uno_asrot.fch'
 end if

 select case(TRIM(gvb_prog))
 case('gamess')
  call do_gvb_gms(proname, pair_fch, .false.)
 case('qchem')
  call do_gvb_qchem(proname, pair_fch)
 case('gaussian')
  call do_gvb_gau(proname, pair_fch)
 case default
  write(6,'(/,A)') 'ERROR in subroutine do_gvb: invalid GVB_prog='//TRIM(gvb_prog)
  write(6,'(A)') 'Currently supported programs: GAMESS, QChem, Gaussian.'
  stop
 end select

 if(mo_rhf) then ! paired LMOs obtained from RHF virtual projection
  write(proname1,'(A,I0)') TRIM(proname)//'_proj_loc_pair2gvb',npair
 else ! paired LMOs obtained from associated rotation of UNOs
  write(proname1,'(A,I0)') TRIM(proname)//'_uno_asrot2gvb',npair
 end if

 sort_fch = TRIM(proname1)//'_s.fch'
 select case(TRIM(localm))
 case('pm')
  pm_loc = .true.
 case('boys')
  pm_loc = .false.
 case default
  write(6,'(/,A)') 'ERROR in subroutine do_gvb: invalid orbital localization me&
                   &thod LocalM='//TRIM(localm)
  stop
 end select

 ! localize singly occ (if any)
 if(nopen > 0) call localize_singly_occ_orb(sort_fch, pm_loc)

 ! localize doubly occ if required
 if(LocDocc) then
  write(6,'(A)') 'Perform orbital localization on doubly occupied orbitals...'
  call localize_orb(sort_fch, 1, ndb, pm_loc)
 end if

 ! exclude X-H bonds with little multi-reference characters from the GVB active
 ! space
 if(excludeXH) then
  datname = TRIM(proname1)//'.dat'
  gmsname = TRIM(proname1)//'.gms'
  call gvb_exclude_XH_A_wrap(datname, gmsname, onlyXH, inpname)
  call get_npair_from_inpname(inpname, i)
  ndb = ndb + npair - i
  npair = i
  call do_gvb_gms(inpname, sort_fch, .true.)
 end if

 ! determine the number of active orbitals/electrons in subsequent CAS/DMRG
 ! computations
 nacta = npair0 + nopen
 nactb = npair0
 nacte = nacta + nactb
 nacto = nacte

 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine do_gvb at '//TRIM(data_string)
end subroutine do_gvb

! perform GVB computation (only in Strategy 1,3) using GAMESS
subroutine do_gvb_gms(proname, pair_fch, name_determined)
 use mr_keyword, only: ist, mem, nproc, gms_path, gms_scr_path, mo_rhf, &
  datname, bgchg, chgname, cart, check_gms_path, GVB_conv, fcgvb
 use mol, only: nbf, nif, ndb, nopen, npair, npair0, gvb_e
 use util_wrapper, only: fch2inp_wrap
 implicit none
 integer :: i, SYSTEM, RENAME
 real(kind=8) :: unpaired_e
 character(len=240) :: inpname, gmsname
 character(len=240), intent(in) :: proname, pair_fch
 character(len=480) :: longbuf = ' '
 logical, intent(in) :: name_determined
 ! True: using the provided proname as inpname
 ! False: create inpname, gmsname, datname according to proname

 call check_gms_path()

 if(name_determined) then
  inpname = proname
  i = INDEX(inpname, '.inp', back=.true.)
  datname = inpname(1:i-1)//'.dat'
  gmsname = inpname(1:i-1)//'.gms'
 else ! not determined
  if(mo_rhf) then ! paired LMOs obtained from RHF virtual projection
   call fch2inp_wrap(pair_fch, .true., npair, nopen)
   if(ist == 6) then
    write(inpname,'(A,I0,A)') TRIM(proname)//'2gvb',npair,'.inp'
    write(gmsname,'(A,I0,A)') TRIM(proname)//'2gvb',npair,'.gms'
    write(datname,'(A,I0,A)') TRIM(proname)//'2gvb',npair,'.dat'
    i = RENAME(TRIM(proname)//'.inp', inpname)
   else
    write(inpname,'(A,I0,A)') TRIM(proname)//'_proj_loc_pair2gvb',npair,'.inp'
    write(gmsname,'(A,I0,A)') TRIM(proname)//'_proj_loc_pair2gvb',npair,'.gms'
    write(datname,'(A,I0,A)') TRIM(proname)//'_proj_loc_pair2gvb',npair,'.dat'
    i = RENAME(TRIM(proname)//'_proj_loc_pair.inp', inpname)
   end if
  else ! paired LMOs obtained from associated rotation of UNOs
   call fch2inp_wrap(pair_fch, .true., npair, nopen)
   write(inpname,'(A,I0,A)') TRIM(proname)//'_uno_asrot2gvb',npair,'.inp'
   write(gmsname,'(A,I0,A)') TRIM(proname)//'_uno_asrot2gvb',npair,'.gms'
   write(datname,'(A,I0,A)') TRIM(proname)//'_uno_asrot2gvb',npair,'.dat'
   i = RENAME(TRIM(proname)//'_uno_asrot.inp', inpname)
  end if
 end if

 call modify_memory_in_gms_inp(inpname, mem, nproc)
 call add_gvb_conv(inpname, GVB_conv)
 if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))

 ! perform GVB with all doubly occupied orbitals frozen
! call add_frz2gms_inp(inpname)
! call submit_gms_fzgvb_job(gms_path, gms_scr_path, inpname, nproc)

 ! modify the input file for GVB with all doubly occupied orbitals frozen
 if(fcgvb) then
  write(6,'(A)') 'Remark: FcGVB=.T. GVB with all doubly occupied orbitals frozen.'
  call add_frz2gms_inp(inpname)
 end if

 ! perform GVB
 call submit_gms_job(gms_path, gms_scr_path, inpname, nproc)

 if(name_determined) write(6,'(A)') 'After excluding inactive X-H pairs from th&
                                    &e original GVB:'
 call read_gvb_energy_from_gms(gmsname, gvb_e)
 write(6,'(/,A,F18.8,1X,A4)') 'E(GVB) = ', gvb_e, 'a.u.'

 ! sort the GVB pairs by CI coefficients of the 1st NOs
 if(cart) then ! Cartesian functions
  write(longbuf,'(A,5(1X,I0))') 'gvb_sort_pairs '//TRIM(datname),nbf,nif,ndb,nopen,npair
 else          ! spherical harmonic functions
  call read_cart_nbf_from_dat(datname, i)
  write(longbuf,'(A,5(1X,I0))') 'gvb_sort_pairs '//TRIM(datname),i,nif,ndb,nopen,npair
 end if

 i = SYSTEM(TRIM(longbuf))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine do_gvb_gms: failed to call utility gvb_&
                   &sort_pairs.'
  write(6,'(A)') 'Did you delete it or forget to compile it?'
  write(6,'(A)') 'Or maybe there is some unexpected error.'
  stop
 end if

 ! generate corresponding .fch file from _s.dat file
 i = INDEX(datname, '.dat')
 inpname = datname(1:i-1)//'_s.fch'
 datname = datname(1:i-1)//'_s.dat'
 call copy_file(pair_fch, inpname, .false.)
 write(longbuf,'(2(A,I0),A)') 'dat2fch '//TRIM(datname)//' '//TRIM(inpname)//&
                              ' -gvb ', npair, ' -open ', nopen, ' >/dev/null'
 i = SYSTEM(TRIM(longbuf))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine do_gvb_gms: failed to call utility dat2&
                   &fch.'
  write(6,'(A)') 'Did you delete it or forget to compile it?'
  stop
 end if

 ! extract NOONs from the above .dat file and print them into .fch file
 write(longbuf,'(A,3(1X,I0),A5)') 'extract_noon2fch '//TRIM(datname)//' '//&
                     TRIM(inpname), ndb+1, ndb+nopen+2*npair, nopen, ' -gau'
 i = SYSTEM(TRIM(longbuf))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine do_gvb_gms: failed to call utility&
                  & extract_noon2fch.'
  write(6,'(A)') 'Did you delete it or forget to compile it?'
  stop
 end if

 ! update Total SCF Density in .fch(k) file
 call update_density_using_no_and_on(inpname)

 ! calculate odd/unpaired electron number
 call calc_unpaired_from_fch(inpname, 2, .false., unpaired_e)

 ! find npair0: the number of active pairs (|C2| > 0.1)
 call find_npair0_from_dat(datname, npair, npair0)
end subroutine do_gvb_gms

! perform GVB computation (only in Strategy 1,3) using QChem
! TODO: change to variational PP. Current PP uses projection. The NO and NOONs
! of variational PP remain to be checked.
subroutine do_gvb_qchem(proname, pair_fch)
 use mr_keyword, only: mem, nproc, mo_rhf, bgchg, chgname, datname
 use mol, only: nopen, npair, npair0, gvb_e
 use util_wrapper, only: fch2qchem_wrap, fch2inp_wrap
 implicit none
 integer :: i, SYSTEM, RENAME
 real(kind=8) :: unpaired_e
 real(kind=8), allocatable :: coeff(:,:)
 character(len=240) :: buf, inpname, outname, fchname0, fchname1, fchname, &
                       pre_fch, gms_inp
 character(len=240), intent(in) :: proname, pair_fch
 character(len=500) :: longbuf
 logical :: alive

 if(mo_rhf) then ! paired LMOs obtained from RHF virtual projection
  pre_fch = TRIM(proname)//'_proj_loc_pair.fch'
  write(buf,'(A,I0)') TRIM(proname)//'_proj_loc_pair2gvb',npair
 else ! paired LMOs obtained from associated rotation of UNOs
  pre_fch = TRIM(proname)//'_uno_asrot.fch'
  write(buf,'(A,I0)') TRIM(proname)//'_uno_asrot2gvb',npair
 end if

 inpname = TRIM(buf)//'.in'
 outname = TRIM(buf)//'.out'
 fchname0= TRIM(buf)//'.0.FChk'
 fchname1= TRIM(buf)//'.FChk'
 fchname = TRIM(buf)//'_s.fch'
 datname = TRIM(buf)//'_s.dat'
 gms_inp = TRIM(buf)//'_s.inp'

 call fch2qchem_wrap(pair_fch, npair, inpname)
 ! Note: nopen is determined automatically in fch2qchem

 if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
 call modify_memory_in_qchem_inp(inpname, mem)
 call submit_qchem_job(inpname, nproc)
 call delete_file(TRIM(fchname0)) ! O.K. even if fchname0 does not exist
 call copy_file(pre_fch, fchname, .false.)

 ! for Q-Chem 5.3.2, only *.fchk file, so rename fchname1
 ! (thanks to bug report from Dr. Qi Ou)
 inquire(file=TRIM(fchname1),exist=alive)
 if(.not. alive) then
  fchname1 = TRIM(buf)//'.fchk'
  inquire(file=TRIM(fchname1),exist=alive)
  if(.not. alive) then
   write(6,'(/,A)') 'ERROR in subroutine do_gvb_qchem: fchname1 does not exist.'
   write(6,'(A)') 'This may be a bug or unexpected version of Q-Chem.'
   write(6,'(A)') 'fchname1 = '//TRIM(fchname1)
   stop
  end if
 end if
 call copy_orb_and_den_in_fch(fchname1, fchname, .false.)
 call delete_file(TRIM(fchname1))

 ! copy NOONs from output file into fch file
 longbuf = 'extract_noon2fch '//TRIM(outname)//' '//TRIM(fchname)
 write(6,'(A)') '$'//TRIM(longbuf)
 i = SYSTEM(TRIM(longbuf))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine do_gvb_qchem: failed to call utility ex&
                   &tract_noon2fch.'
  write(6,'(A)') 'Did you delete it or forget to compile it?'
  stop
 end if

 call fch2inp_wrap(fchname, .true., npair, nopen)
 i = RENAME(TRIM(gms_inp), TRIM(datname))
 allocate(coeff(2,npair))
 call read_pair_coeff_from_qchem_out(outname, npair, coeff)
 call write_pair_coeff_into_gms_inp(datname, npair, coeff)
 call determine_npair0_from_pair_coeff(npair, coeff, npair0)
 deallocate(coeff)

 call read_gvb_e_from_qchem_out(outname, gvb_e)
 write(6,'(/,A,F18.8,1X,A4)') 'E(GVB) = ', gvb_e, 'a.u.'

 call calc_unpaired_from_fch(fchname, 2, .false., unpaired_e)

 write(6,'(A)') REPEAT('-',79)
 write(6,'(A)') 'Remark: the GVB-PP in Q-Chem uses coupled-cluster-like formula&
                &e, see DOI:10.1016/'
 write(6,'(A)') '0301-0104(95)00321-5. It is a non-variational method, so its e&
                &nergy usually'
 write(6,'(A)') 'differs slightly with that calculated by GVB-PP in GAMESS/Gaus&
                &sian. This tiny'
 write(6,'(A)') 'difference does not affect the subsequent calculations.'
 write(6,'(A)') REPEAT('-',79)
end subroutine do_gvb_qchem

! perform GVB computation (only in Strategy 1,3) using Gaussian
subroutine do_gvb_gau(proname, pair_fch)
 use mr_keyword, only: mem, nproc, gau_path, mo_rhf, bgchg, chgname, cart,&
  datname
 use mol, only: nbf, nif, ndb, nopen, npair, npair0, gvb_e
 use util_wrapper, only: unfchk, formchk, fch2inp_wrap
 implicit none
 integer :: i, SYSTEM, RENAME
 real(kind=8) :: unpaired_e
 real(kind=8), allocatable :: coeff(:,:) ! GVB pair coeff
 character(len=240) :: buf, chkname, fchname, gjfname, logname, inpname
 character(len=480) :: longbuf
 character(len=240), intent(in) :: proname, pair_fch

 call check_exe_exist(gau_path)

 if(mo_rhf) then ! paired LMOs obtained from RHF virtual projection
  write(buf,'(A,I0)') TRIM(proname)//'_proj_loc_pair2gvb',npair
 else ! paired LMOs obtained from associated rotation of UNOs
  write(buf,'(A,I0)') TRIM(proname)//'_uno_asrot2gvb',npair
 end if
 chkname = TRIM(buf)//'.chk'
 fchname = TRIM(buf)//'.fch'
 gjfname = TRIM(buf)//'.gjf'
 logname = TRIM(buf)//'.log'
 inpname = TRIM(buf)//'.inp'
 datname = TRIM(buf)//'.dat'
 call unfchk(pair_fch, chkname)
 call prt_gvb_gau_inp(gjfname, mem, nproc, npair)
 if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(gjfname))

 buf = TRIM(gau_path)//' '//TRIM(gjfname)
 write(6,'(A)') '$'//TRIM(buf)
 i = SYSTEM(TRIM(buf))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine do_gvb_gau: Gaussian GVB job failed.'
  write(6,'(A)') 'Please open file '//TRIM(logname)//' and check why.'
  stop
 end if

 call read_gvb_e_from_gau_log(logname, gvb_e)
 write(6,'(/,A,F18.8,1X,A4)') 'E(GVB) = ', gvb_e, 'a.u.'

 call formchk(chkname, fchname)
 call delete_file(chkname)
 call fch2inp_wrap(fchname, .true., npair, nopen)

 i = RENAME(TRIM(inpname), TRIM(datname))
 allocate(coeff(2,npair))
 call read_pair_coeff_from_gau_log(logname, npair, coeff)
 call write_pair_coeff_into_gms_inp(datname, npair, coeff)
 call determine_npair0_from_pair_coeff(npair, coeff, npair0)
 deallocate(coeff)

 ! sort the GVB pairs by CI coefficients of the 1st NOs
 if(cart) then ! Cartesian functions
  write(longbuf,'(A,5(1X,I0))') 'gvb_sort_pairs '//TRIM(datname),nbf,nif,ndb,nopen,npair
 else          ! spherical harmonic functions
  call read_cart_nbf_from_dat(datname, i)
  write(longbuf,'(A,5(1X,I0))') 'gvb_sort_pairs '//TRIM(datname),i,nif,ndb,nopen,npair
 end if
 write(6,'(A)') '$'//TRIM(longbuf)
 i = SYSTEM(TRIM(longbuf))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine do_gvb_gau: failed to call utility&
                  & gvb_sort_pairs.'
  write(6,'(A)') 'Did you delete it or forget to compile it?'
  write(6,'(A)') 'If neither, there is some unexpected error.'
  stop
 end if

 ! generate corresponding .fch file from _s.dat file
 i = INDEX(datname, '.dat')
 inpname = datname(1:i-1)//'_s.fch'
 datname = datname(1:i-1)//'_s.dat'
 call copy_file(fchname, inpname, .false.)
 write(longbuf,'(2(A,I0))') 'dat2fch '//TRIM(datname)//' '//TRIM(inpname)//' -gvb ',&
                             npair, ' -open ', nopen
 write(6,'(A)') '$'//TRIM(longbuf)
 i = SYSTEM(TRIM(longbuf))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine do_gvb_gau: failed to call utility&
                  & dat2fch.'
  write(6,'(A)') 'Did you delete it or forget to compile it?'
  stop
 end if

 ! extract NOONs from the above .dat file and print them into .fch file
 write(longbuf,'(A,3(1X,I0),A5)') 'extract_noon2fch '//TRIM(datname)//' '//&
                     TRIM(inpname), ndb+1, ndb+nopen+2*npair, nopen, ' -gau'
 write(6,'(A)') '$'//TRIM(longbuf)
 i = SYSTEM(TRIM(longbuf))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine do_gvb_gau: failed to call utility&
                   & extract_noon2fch.'
  write(6,'(A)') 'Did you delete it or forget to compile it?'
  write(6,'(A)') 'If neither, there is some unexpected error.'
  stop
 end if

 call calc_unpaired_from_fch(inpname, 2, .false., unpaired_e)
end subroutine do_gvb_gau

! create Gaussian GVB input file
subroutine prt_gvb_gau_inp(gjfname, mem, nproc, npair)
 implicit none
 integer :: i, fid
 integer, intent(in) :: mem, nproc, npair
 character(len=2), allocatable :: pair(:)
 character(len=240) :: chkname
 character(len=240), intent(in) :: gjfname

 i = INDEX(gjfname, '.gjf', back=.true.)
 chkname = gjfname(1:i-1)//'.chk'

 open(newunit=fid,file=TRIM(gjfname),status='replace')
 write(fid,'(A)') '%chk='//TRIM(chkname)
 write(fid,'(A,I0,A)') '%mem=',mem,'GB'
 write(fid,'(A,I0)') '%nprocshared=',nproc
 write(fid,'(A,I0,A,/)') '#p GVB(',npair,') chkbasis nosymm int=nobasistransfor&
                         &m guess=read geom=allcheck scf(maxcycle=300)'
 allocate(pair(npair))
 pair = ' 2'
 write(fid,'(30A2)') (pair(i),i=1,npair)
 write(fid,'(/)',advance='no')
 close(fid)
 deallocate(pair)
end subroutine prt_gvb_gau_inp

! read GVB energy from a Q-Chem output file
subroutine read_gvb_e_from_qchem_out(outname, gvb_e)
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
  write(6,'(/,A)') "ERROR in subroutine read_gvb_e_from_qchem_out: no 'GVB Con&
                   &verge' found in file "//TRIM(outname)
  stop
 end if

 read(buf(16:),*) i, gvb_e
end subroutine read_gvb_e_from_qchem_out

! read GVB energy from a Gaussian output file
subroutine read_gvb_e_from_gau_log(logname, gvb_e)
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
  write(6,'(/,A)') "ERROR in subroutine read_gvb_e_from_gau_log: no 'TOTAL&
                     & ENERGY' found in file "//TRIM(logname)
  stop
 end if

 read(buf(34:),*) gvb_e
end subroutine read_gvb_e_from_gau_log

subroutine read_pair_coeff_from_qchem_out(outname, npair, coeff)
 implicit none
 integer :: i, j, fid, ndb
 integer, intent(in) :: npair
 real(kind=8) :: r1, r2, noon(npair)
 real(kind=8), intent(out) :: coeff(2,npair)
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 coeff = 0d0; buf = ' '
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:7) == 'Alpha:') exit
 end do ! for while

 read(buf(8:),*) ndb

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(4:14) == 'Orbital occ') exit
 end do ! for while

 do i = 1, ndb+2, 1
  read(fid,'(A)') buf
 end do ! for i

 noon = 0
 do i = 1, npair, 1
  read(fid,*) j, r1, r2
  noon(i) = r1 + r2
 end do ! for i

 close(fid)
 forall(i=1:npair)
  coeff(1,i) = DSQRT(noon(i)*0.5d0)
  coeff(2,i) = -DSQRT(1d0 - noon(i)*0.5d0)
 end forall
end subroutine read_pair_coeff_from_qchem_out

! read GVB pair coefficients from a Gaussian output file
subroutine read_pair_coeff_from_gau_log(logname, npair, coeff)
 implicit none
 integer :: i, j, k, fid
 integer, intent(in) :: npair
 real(kind=8), allocatable :: coeff0(:,:)
 real(kind=8), intent(out) :: coeff(2,npair)
 character(len=240) :: buf
 character(len=240), intent(in) :: logname

 coeff = 0d0
 open(newunit=fid,file=TRIM(logname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(2:15) == 'Separated pair') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_pair_coeff_from_gau_log: no&
                    & 'Separated pair' found in file "//TRIM(logname)
  close(fid)
  stop
 end if

 read(fid,'(A)') buf
 do i = 1, npair, 1
  read(fid,'(A)') buf
  j = INDEX(buf, '(')
  read(buf(19:j-1),*) coeff(1,i)
  j = INDEX(buf, ')')
  k = INDEX(buf, '(', back=.true.)
  read(buf(j+1:k-1),*) coeff(2,i)
 end do ! for i

 close(fid)
 ! reverse the order of pair coeff
 allocate(coeff0(2,npair))
 forall(i = 1:npair) coeff0(:,npair-i+1) = coeff(:,i)
 coeff = coeff0
 deallocate(coeff0)
end subroutine read_pair_coeff_from_gau_log

! write the GVB pair coefficients into a given GAMESS .inp/.dat file
subroutine write_pair_coeff_into_gms_inp(datname, npair, coeff)
 implicit none
 integer :: i, fid, fid1, RENAME
 integer, intent(in) :: npair
 real(kind=8), intent(in) :: coeff(2,npair)
 character(len=240) :: buf, datname1
 character(len=240), intent(in) :: datname

 datname1 = TRIM(datname)//'.t'
 open(newunit=fid,file=TRIM(datname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(datname1),status='replace')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:5) == '$SCF') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine write_pair_coeff_into_gms_inp: no&
                     & '$SCF' found in file "//TRIM(datname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if
 i = INDEX(buf, '$END')
 write(fid1,'(A)') buf(1:i-1)

 do i = 1, npair, 1
  write(fid1,'(A,I3,2(A,F18.14))') '   CICOEF(',2*i-1,')=',coeff(1,i),',',coeff(2,i)
 end do ! for i
 write(fid1,'(A)') ' $END'

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(datname1), TRIM(datname))
end subroutine write_pair_coeff_into_gms_inp

! determine npair0 from the GVB pair coefficients
subroutine determine_npair0_from_pair_coeff(npair, coeff, npair0)
 implicit none
 integer :: i
 integer, intent(in) :: npair
 integer, intent(out) :: npair0
 real(kind=8), allocatable :: coeff0(:,:)
 real(kind=8), intent(in) :: coeff(2,npair)

 allocate(coeff0(2,npair))

 do i = 1, npair, 1
  if(DABS(coeff(2,i)) > DABS(coeff(1,i))) then
   coeff0(1,i) = DABS(coeff(2,i))
   coeff0(2,i) = -DABS(coeff(1,i))
  else
   if(coeff(2,i) < 0d0) then
    coeff0(:,i) = coeff(:,i)
   else
    coeff0(:,i) = -coeff(:,i)
   end if
  end if
 end do ! for i

 npair0 = COUNT(coeff0(2,:) < -0.1d0)
 deallocate(coeff0)
end subroutine determine_npair0_from_pair_coeff

! find the number of pairs from a given inpname
subroutine get_npair_from_inpname(inpname, npair)
 implicit none
 integer :: i, j
 integer, intent(out) :: npair
 character(len=240), intent(in) :: inpname

 npair = 0
 i = INDEX(inpname, 'gvb', back=.true.)
 j = INDEX(inpname, '.inp', back=.true.)
 if(i==0 .or. j==0) then
  write(6,'(A)') 'ERROR in subroutine get_npair_from_inpname: invalid filename&
                 &='//TRIM(inpname)
  stop
 end if
 read(inpname(i+3:j-1),*) npair
end subroutine get_npair_from_inpname

! add/modify CONV value in a GAMESS .inp file
! This is GVB SCF density matrix convergence criterion/threshold
subroutine add_gvb_conv(inpname, GVB_conv)
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=240) :: buf, inpname1
 character(len=4), intent(in) :: GVB_conv
 character(len=240), intent(in) :: inpname

 if(GVB_conv=='1d-5' .or. GVB_conv=='1D-5') return
 i =  index(inpname, '.inp', back=.true.)
 inpname1 = inpname(1:i-1)//'.t'

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:5) == '$SCF') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 i = INDEX(buf, '$END')
 if(i > 0) then
  buf = buf(1:i-1)//'CONV='//GVB_conv//' $END'
 else
  buf = TRIM(buf)//' CONV='//GVB_conv
 end if

 if(LEN_TRIM(buf) > 79) then
  i = INDEX(buf(1:79), ' ', back=.true.)
  write(fid1,'(A)') buf(1:i-1)
  write(fid1,'(A)') ' '//TRIM(buf(i+1:))
 else
  write(fid1,'(A)') TRIM(buf)
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine add_gvb_conv

! add keywords of freezing all doubly orbitals for a GAMESS .inp file
subroutine add_frz2gms_inp(inpname)
 implicit none
 integer :: i, ncore, fid, fid1, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 i = INDEX(inpname, '.inp', back=.true.)
 inpname1 = inpname(1:i-1)//'_f.inp'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid,'(A)') buf
  write(fid1,'(A)') TRIM(buf)
  if(buf(2:5) == '$SCF') exit
 end do ! for while

 i = INDEX(buf, 'NCO')
 if(i == 0) then
  close(fid)
  close(fid1, status='delete')
  write(6,'(/,A)') "ERROR in subroutine add_frz2gms_inp: no 'NCO' found in file"
  write(6,'(A)') trim(inpname)
  stop
 end if
 read(buf(i+4:),*) ncore

 if(INDEX(buf,'$END') == 0) then
  do while(.true.)
   read(fid,'(A)') buf
   write(fid1,'(A)') TRIM(buf)
   if(INDEX(buf,'$END') > 0) exit
  end do ! for while
 end if

 if(ncore > 0) then
  write(fid1,'(A)',advance='no') ' $MOFRZ FRZ=.T. IFRZ(1)='
  if(ncore > 1) then
   write(fid1,'(19(I0,A1))') (i,',',i=1,ncore)
  else ! ncore == 1
   write(fid1,'(A)') '1'
  end if
  write(fid1,'(A)') ' $END'

  read(fid,'(A)') buf
  if(buf(2:7) == '$MOFRZ') then
   if(INDEX(buf,'$END') == 0) then
    do while(.true.)
     read(fid,'(A)') buf
     if(buf(2:5) == '$END') exit
    end do ! for while
   end if
  else
   write(fid1,'(A)') TRIM(buf)
  end if
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine add_frz2gms_inp

