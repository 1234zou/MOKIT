! written by jxzou at 20210729: move subroutine do_gvb from file automr.f90 to
! this new file, and rename it as do_gvb_gms

! perform GVB computation (only in Strategy 1,3) using GAMESS/QChem/Gaussian
subroutine do_gvb()
 use print_id, only: iout
 use mr_keyword, only: gvb, gvb_prog, ist, hf_fch, mo_rhf, npair_wish
 use mol, only: nbf, nif, ndb, npair, nopen, lin_dep
 implicit none
 integer :: i
 character(len=24) :: data_string = ' '
 character(len=240) :: proname, pair_fch

 if(.not. gvb) return
 write(iout,'(//,A)') 'Enter subroutine do_gvb...'
 if(.not. (ist==1 .or. ist==3)) then
  write(iout,'(/,A,I0)') 'ERROR in subroutine do_gvb: ist=', ist
  write(iout,'(A)') 'Only ist=1 or 3 is supported currently.'
  stop
 end if

 i = index(hf_fch, '.fch', back=.true.)
 proname = hf_fch(1:i-1)

 ! In RHF virtual MO projection, it will generate a file uno.out additionally
 call read_npair_from_uno_out(nbf, nif, ndb, npair, nopen, lin_dep)

 if(mo_rhf) then ! paired LMOs obtained from RHF virtual projection
  if(npair_wish>0 .and. npair_wish/=npair) then
   write(iout,'(A)') 'ERROR in subroutine do_gvb: npair_wish cannot be assigned in this strategy.'
   stop
  end if
  pair_fch = TRIM(proname)//'_proj_loc_pair.fch'

 else ! paired LMOs obtained from associated rotation of UNOs
  if(npair_wish>0 .and. npair_wish/=npair) then
   write(iout,'(2(A,I0),A)') 'Warning: AutoMR recommends GVB(',npair,'), but&
    & user specifies GVB(',npair_wish,'). Try to fulfill...'
   if(npair_wish < npair) then
    ndb = ndb + npair - npair_wish
    npair = npair_wish
    write(iout,'(A)') 'OK, fulfilled.'
   else if(npair_wish > npair) then
    write(iout,'(A)') 'ERROR in subroutine do_gvb: too large pairs specified.&
                     & Cannot fulfilled.'
    stop
   end if
  end if

  pair_fch = TRIM(proname)//'_uno_asrot.fch'
 end if

 select case(TRIM(gvb_prog))
 case('gamess')
  call do_gvb_gms(proname, pair_fch)
 case('qchem')
  call do_gvb_qchem(proname, pair_fch)
 case('gaussian')
  call do_gvb_gau(proname, pair_fch)
 case default
  write(iout,'(/,A)') 'ERROR in subroutine do_gvb: invalid GVB_prog='//&
                       TRIM(gvb_prog)
  write(iout,'(A)') 'Currently supported programs: GAMESS, QChem, Gaussian.'
  stop
 end select

 call fdate(data_string)
 write(iout,'(A)') 'Leave subroutine do_gvb at '//TRIM(data_string)
 return
end subroutine do_gvb

! perform GVB computation (only in Strategy 1,3) using GAMESS
subroutine do_gvb_gms(proname, pair_fch)
 use print_id, only: iout
 use mr_keyword, only: mem, nproc, gms_path, gms_scr_path, mo_rhf, &
  ist, datname, bgchg, chgname, cart, check_gms_path
 use mol, only: nbf, nif, ndb, nopen, npair, gvb_e, nacto, nacta, &
  nactb, nacte, npair0
 implicit none
 integer :: i, j, system, RENAME
 character(len=240) :: buf, chkname, inpname, gmsname
 character(len=240), intent(in) :: proname, pair_fch
 character(len=300) :: longbuf = ' '

 call check_gms_path()

 if(mo_rhf) then ! paired LMOs obtained from RHF virtual projection
  write(buf,'(2(A,I0))') 'fch2inp '//TRIM(pair_fch)//' -gvb ',npair,' -open ',nopen
  write(iout,'(A)') '$'//TRIM(buf)
  i = system(TRIM(buf))
  write(inpname,'(A,I0,A)') TRIM(proname)//'_proj_loc_pair2gvb',npair,'.inp'
  write(gmsname,'(A,I0,A)') TRIM(proname)//'_proj_loc_pair2gvb',npair,'.gms'
  write(datname,'(A,I0,A)') TRIM(proname)//'_proj_loc_pair2gvb',npair,'.dat'
  i = RENAME(TRIM(proname)//'_proj_loc_pair.inp', inpname)

 else ! paired LMOs obtained from associated rotation of UNOs
  write(buf,'(2(A,I0))') 'fch2inp '//TRIM(pair_fch)//' -gvb ',npair,' -open ',nopen
  write(iout,'(A)') '$'//TRIM(buf)
  i = system(TRIM(buf))
  write(inpname,'(A,I0,A)') TRIM(proname)//'_uno_asrot2gvb',npair,'.inp'
  write(gmsname,'(A,I0,A)') TRIM(proname)//'_uno_asrot2gvb',npair,'.gms'
  write(datname,'(A,I0,A)') TRIM(proname)//'_uno_asrot2gvb',npair,'.dat'
  i = RENAME(TRIM(proname)//'_uno_asrot.inp', inpname)
 end if

 call modify_memory_in_gms_inp(inpname, mem, nproc)

 ! call GAMESS to do GVB computations (delete .dat file first, if any)
 buf = TRIM(gms_scr_path)//'/'//TRIM(datname)
 call delete_file(buf)
 if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))

 write(longbuf,'(A,I0,A)') TRIM(inpname)//' 01 ',nproc,' >'//TRIM(gmsname)//" 2>&1"
! write(iout,'(A)') '$$GMS '//TRIM(longbuf)
 i = system(TRIM(gms_path)//' '//TRIM(longbuf))

 call read_gvb_energy_from_gms(gmsname, gvb_e)
 write(iout,'(/,A,F18.8,1X,A4)') 'E(GVB) = ', gvb_e, 'a.u.'

 ! move the .dat file into current directory
 i = system('mv '//TRIM(gms_scr_path)//'/'//TRIM(datname)//' .')
 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine do_gvb: fail to move file. Possibly&
                   & wrong gms_scr_path.'
  write(iout,'(A)') 'gms_scr_path='//TRIM(gms_scr_path)
  stop
 end if

 ! sort the GVB pairs by CI coefficients of the 1st NOs
 if(cart) then ! Cartesian functions
  write(longbuf,'(A,5(1X,I0))') 'gvb_sort_pairs '//TRIM(datname),nbf,nif,ndb,nopen,npair
 else          ! spherical harmonic functions
  call read_nbf_from_dat(datname, i)
  write(longbuf,'(A,5(1X,I0))') 'gvb_sort_pairs '//TRIM(datname),i,nif,ndb,nopen,npair
 end if
 i = system(TRIM(longbuf))
 if(i /= 0) then
  write(iout,'(/,A)') 'ERROR in subroutine do_gvb: failed to call utility gvb_sort_pairs.'
  write(iout,'(A)') 'Did you delete it or forget to compile it?'
  write(iout,'(A)') 'Or maybe there is some unexpected error.'
  stop
 end if

 ! generate corresponding .fch file from _s.dat file
 i = index(datname, '.dat')
 inpname = datname(1:i-1)//'_s.fch'
 datname = datname(1:i-1)//'_s.dat'
 call copy_file(pair_fch, inpname, .false.)
 write(longbuf,'(2(A,I0))') 'dat2fch '//TRIM(datname)//' '//TRIM(inpname)//' -gvb ',&
                             npair, ' -open ', nopen
 i = system(TRIM(longbuf))
 if(i /= 0) then
  write(iout,'(/,A)') 'ERROR in subroutine do_gvb: failed to call utility dat2fch.'
  write(iout,'(A)') 'Did you delete it or forget to compile it?'
  stop
 end if

 ! extract NOONs from the above .dat file and print them into .fch file
 write(longbuf,'(A,3(1X,I0),A5)') 'extract_noon2fch '//TRIM(datname)//' '//&
                     TRIM(inpname), ndb+1, ndb+nopen+2*npair, nopen, ' -gau'
 i = system(TRIM(longbuf))
 if(i /= 0) then
  write(iout,'(/,A)') 'ERROR in subroutine do_gvb: failed to call utility extract_noon2fch.'
  write(iout,'(A)') 'Did you delete it or forget to compile it?'
  stop
 end if

 ! update Total SCF Density in .fch(k) file
 call update_density_using_no_and_on(inpname)

 ! find npair0: the number of active pairs (|C2| > 0.1)
 call find_npair0_from_dat(datname, npair, npair0)

 ! determine the number of orbitals/electrons in following CAS/DMRG computations
 nacta = npair0 + nopen
 nactb = npair0
 nacte = nacta + nactb
 nacto = nacte

 return
end subroutine do_gvb_gms

! perform GVB computation (only in Strategy 1,3) using QChem
subroutine do_gvb_qchem(proname, pair_fch)
 use print_id, only: iout
 implicit none
 character(len=240), intent(in) :: proname, pair_fch

 write(iout,'(A)') 'ERROR in subroutine do_gvb_qchem: implmentation not&
                  & finished yet.'
 stop
 return
end subroutine do_gvb_qchem

! perform GVB computation (only in Strategy 1,3) using Gaussian
subroutine do_gvb_gau(proname, pair_fch)
 use print_id, only: iout
 implicit none
 character(len=240), intent(in) :: proname, pair_fch

 return
end subroutine do_gvb_gau

