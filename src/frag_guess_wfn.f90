! written by jxzou at 20201121: perform fragment-guess wavefunction calculations
! updated by jxzou at 20210409: add support of MOROKUMA-EDA and GKS-EDA
! updated by jxzou at 20210414: automatically add $DFTTYP and $PCM into GAMESS .inp file
! updated by jxzou at 20210417: construct supermolecule MOs from direct sum of fragments MOs

! The fragment-guess wavefunction calculation in Gaussian has the following
! shortcomings:
! (1) cannot ensure the wavefunction stability of each fragment. This often
!     makes the calculations of transition-metal molecules take more SCF cycles,
!     since the initial guess is not as good as we expect.
! (2) do not support mixed methods, e.g. cannot use RHF/RDFT for some fragments
!     and UHF/UDFT for other fragments
! This utility aims to overcome these shortcomings.

module frag_info
 implicit none
 integer :: nfrag0 ! No. of fragments in .gjf file
 integer :: nfrag  ! No. of fragments (including total system) to be calculated
 ! For none-EDA and Morokuma-EDA calculations, nfrag = nfrag0
 ! For GKS-EDA calculations, nfrag = 2*nfrag0 + 1

 type :: frag
  integer :: charge = 0
  integer :: mult = 1
  integer :: natom = 0
  integer :: wfn_type = 0            ! 0/1/2/3 for undetermined/RHF/ROHF/UHF
  integer, allocatable :: atm_map(:) ! map to parent system, i.e. supermolecule
  real(kind=8) :: e = 0d0            ! electronic energy
  real(kind=8) :: ssquare = 0d0      ! S(S+1)
  real(kind=8), allocatable :: coor(:,:)   ! size natom
  character(len=2), allocatable :: elem(:) ! size natom
  character(len=240) :: fname = ' '        ! gjf filename
  logical :: pos = .true.            ! .True./.False. for alpha/beta total spin
  logical :: noiter = .false.        ! .True./.False. for skipping SCF or not
  logical, allocatable :: ghost(:)   ! .True./.False. for ghost atoms or not, size natom
 end type frag

 type(frag), allocatable :: frags(:)
contains

 subroutine copy_frag(frag1, frag2)
  implicit none
  type(frag), intent(in) :: frag1
  type(frag), intent(out) :: frag2

  frag2%charge = frag1%charge
  frag2%mult = frag1%mult
  frag2%natom = frag1%natom
  frag2%wfn_type = frag1%wfn_type
  frag2%atm_map = frag1%atm_map
  frag2%e = frag1%e
  frag2%ssquare = frag1%ssquare
  frag2%coor = frag1%coor
  frag2%elem = frag1%elem
  frag2%fname = frag1%fname
  frag2%pos = frag1%pos
  frag2%noiter = frag1%noiter
  frag2%ghost = frag1%ghost
  return
 end subroutine copy_frag

end module frag_info

module theory_level
 implicit none
 integer :: mem, nproc
 integer :: wfn_type = 0   ! 0/1/2/3 for undetermined/RHF/ROHF/UHF
 integer :: eda_type = 0   ! 0/1/2/3/4 for none/Morokuma-EDA/LMO-EDA/GKS-EDA/SAPT
 integer :: disp_type = 0  ! 0/1/2 for no dispersion correction/GD3/GD3BJ
 character(len=9) :: method = ' '
 character(len=21) :: basis = ' '
 character(len=40) :: scrf  = ' '
 character(len=20) :: solvent = ' '
 logical :: sph = .true.  ! .True./.False. for spherical harmonic/Cartesian functions
end module theory_level

program main
 use print_id, only: iout
 implicit none
 integer :: i
 character(len=240) :: gau_path, gjfname

 i = iargc()
 if(i /= 1) then
  write(iout,'(/,A)') ' ERROR in subroutine frag_guess_wfn: wrong command line&
                      & arguments!'
  write(iout,'(/,A)') " Example 1 (in bash): frag_guess_wfn water_dimer.gjf >& water_dimer.out &"
  write(iout,'(A,/)') " Example 2 (in dash): frag_guess_wfn water_dimer.gjf >water_dimer.out 2>&1 &"
  stop
 end if

 gjfname = ' '
 call getarg(1, gjfname)
 if(index(gjfname,'.gjf') == 0) then
  write(iout,'(/,A)') "ERROR in subroutine frag_guess_wfn: '.gjf' suffix not&
                     & found in filename "//TRIM(gjfname)
  stop
 end if
 call require_file_exist(gjfname)

 call get_gau_path(gau_path)
 call frag_guess_wfn(gau_path, gjfname)
 call gen_inp_of_frags()
 stop
end program main 

! perform fragment-guess wavefunction calculations, one by one fragment
subroutine frag_guess_wfn(gau_path, gjfname)
 use print_id, only: iout
 use util_wrapper, only: formchk, unfchk
 use frag_info, only: nfrag0, nfrag, frag, frags, copy_frag
 use theory_level
 implicit none
 integer :: i, j, k, m, fid, ifrag, iatom
 integer :: charge, mult, natom, maxfrag
 integer, allocatable :: cm(:), nuc(:)
 real(kind=8), allocatable :: coor(:,:), tmp_coor(:,:)
 character(len=2), allocatable :: elem(:)
 character(len=240) :: buf, chkname, fchname, logname, basname
 character(len=240), intent(in) :: gau_path, gjfname
 character(len=1200) :: longbuf
 logical :: guess_read
 type(frag) :: tmp_frag1, tmp_frag2

 buf = ' '; longbuf = ' '
 call read_eda_type_from_gjf(gjfname, eda_type)
 call read_disp_ver_from_gjf(gjfname, disp_type)
 if(eda_type==4 .and. disp_type>0) then
  write(iout,'(/,A)') 'ERROR in subroutine frag_guess_wfn: dispersion correction&
                     & cannot be used in SAPT.'
  write(iout,'(A)') 'It may be supported in DFT-SAPT, but this method is not&
                   & supported for frag_guess_wfn.'
  stop
 end if
 call check_sph_in_gjf(gjfname, sph)

 if(sph .and. eda_type==1) then
  write(iout,'(A)') 'Warning in subroutine frag_guess_wfn: spherical harmonic&
                   & functions (5D 7F) cannot be'
  write(iout,'(A)') 'used in Morokuma-EDA method in GAMESS. Automatically switch&
                   & to Cartesian functions (6D 10F).'
  sph = .false.
 end if

 call read_mem_and_nproc_from_gjf(gjfname, mem, nproc)
 write(iout,'(A,I0,A)') '%mem=', mem, 'MB'
 write(iout,'(A,I0)') '%nprocshared=', nproc
 call read_natom_from_gjf(gjfname, natom) ! read the total number of atoms
 write(iout,'(A,I0)') 'natom=', natom

 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:1) == '#') exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine frag_guess_wfn: incomplete file '//TRIM(gjfname)
  write(iout,'(A)') "Failed to locate the Route Section '#' in .gjf file."
  close(fid)
  stop
 end if

 longbuf = TRIM(buf)
 do while(.true.)
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
  longbuf = TRIM(longbuf)//' '//TRIM(buf)
 end do ! for while

 call lower(longbuf)
 call read_scrf_string_from_buf(longbuf, scrf)
 if(LEN_TRIM(scrf)>0 .and. eda_type==4) then
  write(iout,'(/,A)') 'ERROR in subroutine frag_guess_wfn: SCRF detected.'
  write(iout,'(A)') 'Implicit solvent model is not supported in SAPT computations.'
  stop
 end if
 call determine_solvent_from_gau2gms(scrf, solvent)
 call read_nfrag_from_buf(longbuf, nfrag0)

 select case(eda_type)
 case(0) ! fragment guess wfn
  nfrag = nfrag0 + 1
  maxfrag = 998
 case(1) ! MOROKUMA-EDA
  nfrag = nfrag0 + 1
  maxfrag = 998
 case(2) ! LMO-EDA
  nfrag = 2*nfrag0 + 1
  maxfrag = 499
 case(3) ! GKS-EDA
  nfrag = 2*nfrag0 + 1
  maxfrag = 499
 case(4) ! SAPT
  if(nfrag0 /= 2) then
   write(iout,'(A,I0)') 'ERROR in subroutine frag_guess_wfn: invalid nfrag0=',nfrag0
   write(iout,'(A)') 'Only two fragments/monomers are allowed in SAPT.'
   stop
  end if
  nfrag = 3
  maxfrag = 3
 case default
  write(iout,'(A,I0)') 'ERROR in subroutine frag_guess_wfn: invalid eda_type=',eda_type
  write(iout,'(A,I0)') 'Only 0,1,2,3,4 are allowed.'
  close(fid)
  stop
 end select

 if(nfrag0 > maxfrag) then
  write(iout,'(A)') 'ERROR in subroutine frag_guess_wfn: nfrag0 > maxfrag, too large!'
  write(iout,'(2(A,I0))') 'nfrag0=', nfrag0, ', maxfrag=', maxfrag
  close(fid)
  stop
 end if

 write(iout,'(A,I0)') 'nfrag0=', nfrag0
 write(iout,'(A,I0)') 'nfrag=', nfrag

 call read_method_and_basis_from_buf(longbuf, method, basis, wfn_type)
 write(iout,'(A,I0)') 'method='//TRIM(method)//', basis='//TRIM(basis)//&
                     &', wfn_type=', wfn_type

 if(index(basis,'gen') > 0) then
  close(fid)
  call record_gen_basis_in_gjf(gjfname, basname)
  open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
  do while(.true.)
   read(fid,'(A)') buf
   if(LEN_TRIM(buf) == 0) exit
  end do ! for while
 else
  basname = ' '
 end if

 do while(.true.)
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
 end do ! for while
 allocate(cm(2+2*nfrag0),source=0)

 read(fid,*,iostat=i) cm
 if(i /= 0) then
  write(iout,'(/,A)') 'ERROR in subroutine frag_guess_wfn: incomplete charges and&
                     & multiplicities detected'
  write(iout,'(A)') 'in file '//TRIM(gjfname)
  close(fid)
  stop
 end if

 charge = cm(1); mult = cm(2)

 if(mult/=1 .and. eda_type==1) then
  write(iout,'(A)') 'ERROR in subroutine frag_guess_wfn: Morokuma-EDA can only&
                   & be applied to RHF. But'
  write(iout,'(A)') 'the total spin is not singlet.'
  close(fid)
  stop
 end if

 ! compare the spin and method name, determine RHF/ROHF/UHF-type wfn
 if(mult == 1) then
  if(wfn_type == 0) then
   wfn_type = 1 ! RHF
  else if(wfn_type == 2) then
   write(iout,'(A)') 'ERROR in subroutine frag_guess_wfn: total spin is singlet.&
                    & But you specify ROHF/RODFT.'
   close(fid)
   stop
  end if
 else ! mult /= 1
  if(wfn_type == 0) then
   wfn_type = 3 ! UHF
  else
   if(wfn_type == 1) then
    write(iout,'(A)') 'ERROR in subroutine frag_guess_wfn: total spin is non-singlet.&
                     & But you specify RHF/RDFT.'
    stop
   end if
  end if
 end if

 if(eda_type==1 .and. wfn_type/=1) then
  write(iout,'(A)') 'ERROR in subroutine frag_guess_wfn: Morokuma-EDA can only&
                   & be applied to RHF. But'
  write(iout,'(A,I0)') 'RHF-type wavefunction is not specified. wfn_type=',wfn_type
  close(fid)
  stop
 end if

 ! Initially, the restricted/unrestricted type of wfn of each fragment is
 !  inherited from the method of the total system.
 ! This type may be updated when reading Cartesian coordinates
 allocate(frags(nfrag))
 frags(:)%wfn_type = wfn_type

 forall(i=1:nfrag0)
  frags(i)%charge = cm(2*i+1)
  frags(i)%mult = cm(2*i+2)
 end forall
 deallocate(cm)

 i = SUM(frags(1:nfrag0)%charge)
 if(i /= charge) then
  write(iout,'(/,A)') 'ERROR in subroutine frag_guess_wfn: sum of fragment&
                     & charges is not equal to total charge.'
  write(iout,'(2(A,I0))') 'Total charge=', charge, ', sum(frag_charges)=', i
  write(iout,'(A)') 'Wrong charges in file '//TRIM(gjfname)
  deallocate(frags)
  close(fid)
  stop
 end if
 if(ANY(frags(:)%mult==0) .or. ANY(frags(:)%mult==-1)) then
  write(iout,'(/,A)') 'ERROR in subroutine frag_guess_wfn: the spin multiplicity&
                     & of some fragment is 0 or -1.'
  write(iout,'(A)') 'This is impossible. Check your file '//TRIM(gjfname)
  deallocate(frags)
  close(fid)
  stop
 end if
 j = 0 ! ne(alpha) - ne(beta)
 do i = 1, nfrag0, 1
  k = frags(i)%mult
  if(k > 0) then
   j = j + k - 1
  else
   j = j + k + 1
  end if
 end do ! for i
 if(mult-1 /= j) then
  write(iout,'(/,A)') 'ERROR in subroutine frag_guess_wfn: number of unpaired&
                     & electrons calculated from'
  write(iout,'(A)') 'fragments spin multiplicities is not equal to that calcu&
                    &lated from the total spin'
  write(iout,'(A)') 'multiplicity. Some spin multiplicity must be wrong.'
  write(iout,'(A)') 'Check your spin in file '//TRIM(gjfname)
  deallocate(frags)
  close(fid)
  stop
 end if

 ! now make all spin of fragment positive, move negative(if any) into pos
 do i = 1, nfrag0, 1
  k = frags(i)%mult
  if(k < 0) then
   frags(i)%mult = -k
   frags(i)%pos = .false.
  end if
 end do ! for i

 allocate(coor(3,natom), source=0d0)
 allocate(elem(natom))
 elem = '  '

 do i = 1, natom, 1
  read(fid,'(A)',iostat=j) buf
  if(j/=0 .or. LEN_TRIM(buf)==0) exit

  j = index(buf,'(')
  k = index(buf,'=')
  m = index(buf,')')
  if(j*k*m == 0) then
   write(iout,'(/,A)') 'ERROR in subroutine frag_guess_wfn: wrong format in&
                      & file '//TRIM(gjfname)//'.'
   write(iout,'(A)') 'Please read examples in Section 5.3.2 ~ 5.3.4 of MOKIT&
                    & manual and learn how to'
   write(iout,'(A)') 'write a valid input file.'
   stop
  end if
  read(buf(1:j-1),*) elem(i)
  read(buf(k+1:m-1),*) ifrag

  if(i==1 .and. eda_type/=0 .and. ifrag/=1) then
   write(iout,'(/,A)') 'ERROR in subroutine frag_guess_wfn: error definition of fragment number.'
   write(iout,'(A)') 'EDA input file required. This task requires definition&
                    & of fragment number'
   write(iout,'(A)') 'monomer1, monomer2, monomer3,...'
   stop
  end if

  if(ifrag > nfrag0) then
   write(iout,'(A)') 'ERROR in subroutine frag_guess_wfn: the number of fragments&
                    & detected in Cartesian coordinates is'
   write(iout,'(A)') 'larger than that found in charges and spin multiplicities.'
   write(iout,'(2(A,I0))') 'ifrag = ', ifrag, ', nfrag0=', nfrag0
   stop
  end if

  read(buf(m+1:),*) coor(1:3,i)
  frags(ifrag)%natom = frags(ifrag)%natom + 1
  iatom = frags(ifrag)%natom

  if(iatom == 1) then
   allocate(frags(ifrag)%coor(3,1))
   frags(ifrag)%coor(:,1) = coor(:,i)

   frags(ifrag)%elem = [elem(i)]
   frags(ifrag)%atm_map = [i]
  else
   tmp_coor = frags(ifrag)%coor
   deallocate(frags(ifrag)%coor)
   allocate(frags(ifrag)%coor(3,iatom))
   frags(ifrag)%coor(:,1:iatom-1) = tmp_coor
   deallocate(tmp_coor)
   frags(ifrag)%coor(:,iatom) = coor(:,i)

   frags(ifrag)%elem = [frags(ifrag)%elem, elem(i)]
   frags(ifrag)%atm_map =[frags(ifrag)%atm_map, [i]]
  end if

  j = LEN_TRIM(buf)
  if(buf(j:j)=='r' .or. buf(j:j)=='R') then
   if(frags(ifrag)%mult == 1) then
    frags(ifrag)%wfn_type = 1 ! RHF
   else
    frags(ifrag)%wfn_type = 2 ! ROHF
   end if
  end if

  if(wfn_type==1 .and. frags(ifrag)%wfn_type/=1) then
   close(fid)
   write(iout,'(/,A)') 'ERROR in subroutine frag_guess_wfn: conflicted type of&
                      & wfn between the total system and some fragment.'
   write(iout,'(A)') 'Restricted closed shell wfn is specified for the total&
                    & system. But'
   write(iout,'(A)') 'restricted open shell or unrestricted wfn is specified&
                    & for some fragment.'
   write(iout,'(3(A,I0))') 'natom=', natom, ', i=', i, ', ifrag=', ifrag
   stop
  end if
 end do ! for i

 close(fid)

 do i = 1, nfrag0, 1
  iatom = frags(i)%natom
  if(iatom == 0) then
   write(iout,'(/,A,I0,A)') 'ERROR in subroutine frag_guess_wfn: the ',i,'-th&
    & fragment has 0 atom.'
   write(iout,'(A)') 'Please check your specification of fragments in file '//TRIM(gjfname)
   stop
  end if
  allocate(frags(i)%ghost(iatom))
  frags(i)%ghost = .false.
 end do ! for i

 ! this is the super molecule/total system
 frags(nfrag)%natom = natom
 allocate(frags(nfrag)%elem(natom), frags(nfrag)%coor(3,natom), nuc(natom), &
          frags(nfrag)%atm_map(natom))
 forall(i = 1:natom) frags(nfrag)%atm_map(i) = i
 call read_elem_and_coor_from_gjf(gjfname, natom, frags(nfrag)%elem, nuc, &
  frags(nfrag)%coor, frags(nfrag)%charge, frags(nfrag)%mult)
 deallocate(nuc)
 iatom = frags(nfrag)%natom
 allocate(frags(nfrag)%ghost(iatom))
 frags(nfrag)%ghost = .false.
 frags(nfrag)%noiter = .true.

 if(eda_type == 3) then ! GKS-EDA
  do i = 1, nfrag0, 1
   call gen_extend_bas_frag(frags(i), frags(nfrag), frags(i+nfrag0))
  end do ! for i
 else if(eda_type == 4) then ! SAPT
  call gen_extend_bas_frag(frags(1), frags(3), tmp_frag1)
  call gen_extend_bas_frag(frags(2), frags(3), tmp_frag2)
  call copy_frag(tmp_frag1, frags(1)) ! overwrite frags(1)
  call copy_frag(tmp_frag2, frags(2)) ! overwrite frags(2)
 end if

 ! generate these SCF .gjf files
 do i = 1, nfrag, 1
  write(frags(i)%fname,'(I3.3,A)') i, '-'//TRIM(gjfname)
  guess_read = .false.
!  if(eda_type==3 .and. i>nfrag0 .and. i<nfrag) guess_read = .true.
  call gen_gjf_from_type_frag(frags(i), guess_read, basname)
 end do ! for i

 j = index(frags(1)%fname, '.gjf', back=.true.)

 ! do SCF computations one by one
 do i = 1, nfrag, 1
  call do_scf_and_read_e(gau_path, gau_path, frags(i)%fname, frags(i)%noiter, &
                         frags(i)%e, frags(i)%ssquare)
  if(i < nfrag) then
   write(iout,'(A,I3,A,F18.9,A,F6.2)') 'i=', i, ', frags(i)%e = ', frags(i)%e,&
                                       ', frags(i)%ssquare=', frags(i)%ssquare
  end if
  chkname = frags(i)%fname(1:j-1)//'.chk'
  fchname = frags(i)%fname(1:j-1)//'.fch'
  buf     = frags(i)%fname(1:j-1)//'.gjf'
  logname = frags(i)%fname(1:j-1)//'.log'
  call delete_file(logname)
  if(i < nfrag) call delete_file(buf)
  if(i==nfrag .and. eda_type==1) then ! MOROKUMA-EDA
   call formchk(chkname, fchname)
   call delete_files(2, [chkname, buf])
  end if
 end do ! for i

 ! For GKS-EDA, construct supermolecule MOs from direct sum of fragment MOs
 if(eda_type == 3) then
  i = nfrag
  call direct_sum_frag_mo2super_mo(nfrag0, frags(1:nfrag0)%fname, &
        frags(1:nfrag0)%wfn_type, frags(1:nfrag0)%pos, frags(i)%fname, &
        frags(i)%wfn_type)
  frags(i)%noiter = .false.
  call do_scf_and_read_e(gau_path, gau_path, frags(i)%fname, frags(i)%noiter, &
                         frags(i)%e, frags(i)%ssquare)
  write(iout,'(A,I3,A,F18.9,A,F6.2)') 'i=', i, ', frags(i)%e = ', frags(i)%e,&
                                      ', frags(i)%ssquare=', frags(i)%ssquare

 else if(eda_type == 4) then ! For SAPT, add SCF density of two fragments to
  ! obtain approximate total SCF density
  call sum_frag_density_and_prt_into_fch(2, frags(1:2)%fname, frags(3)%fname)
  k = 1
  if(frags(3)%wfn_type == 3) k = 0
  call gen_no_using_density_in_fch(fchname, k)
  call unfchk(fchname, chkname)
  frags(3)%noiter = .false.
  call modify_guess_only_in_gjf(frags(3)%fname)
  call do_scf_and_read_e(gau_path, gau_path, frags(3)%fname, frags(3)%noiter, &
                         frags(3)%e, frags(3)%ssquare)
  write(iout,'(A,F18.9,A,F6.2)') 'i=  3, frags(i)%e = ', frags(3)%e, &
                                 ', frags(i)%ssquare=', frags(3)%ssquare
 end if
 call delete_files(2, [frags(nfrag)%fname, logname])

 if(eda_type == 4) then ! SAPT in PSI4
  write(iout,'(A)') 'All SCF computations finished. Generating PSI4 input file...'
 else
  write(iout,'(A)') 'All SCF computations finished. Generating GAMESS input file...'
 end if
 return
end subroutine frag_guess_wfn

! read the parameter eda_type from a given .gjf file
subroutine read_eda_type_from_gjf(gjfname, eda_type)
 use print_id, only: iout
 implicit none
 integer :: i, j, k, fid
 integer, intent(out) :: eda_type
 character(len=8) :: str
 character(len=240) :: buf
 character(len=240), intent(in) :: gjfname

 eda_type = 0
 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
 end do ! for while

 read(fid,'(A)') buf ! Title Card line
 close(fid)

 i = index(buf, '{')
 j = index(buf, '}')
 k = index(buf,',')
 if(k > 0) j = k
 if(j < i) then
  write(iout,'(A)') 'ERROR in subroutine read_eda_type_from_gjf: wrong syntax&
                    & in file '//TRIM(gjfname)
  write(iout,'(A)') 'buf='//TRIM(buf)
  stop
 end if

 str = buf(i+1:j-1)
 call lower(str)

 select case(str)
 case('frag')
 case('morokuma')
  eda_type = 1
 case('lmo')
  eda_type = 2
 case('gks')
  eda_type = 3
 case('sapt')
  eda_type = 4
 case default
  write(iout,'(/,A)') "ERROR in subroutine read_eda_type_from_gjf: invalid&
                     & keyword '"//buf(i+1:j-1)//"' in Title Card"
  write(iout,'(A)') 'line of file '//TRIM(gjfname)//'.'
  write(iout,'(A)') 'You are supposed to write one of {morokuma}, {gks} or {&
                    &sapt,bronze} in the'
  write(iout,'(A)') 'Title Card line.'
  stop
 end select

 return
end subroutine read_eda_type_from_gjf

! read memory and nprocshared from a given .gjf file
subroutine read_mem_and_nproc_from_gjf(gjfname, mem, np)
 use print_id, only: iout
 implicit none
 integer :: i, j, k, fid
 integer, intent(out) :: mem, np
 character(len=240) :: buf
 character(len=240), intent(in) :: gjfname

 ! default settings
 mem = 1000 ! 1000 MB
 np = 1     ! 1 core

 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  call lower(buf)
  if(buf(1:1) == '#') exit

  j = LEN_TRIM(buf)
  k = index(buf,'=')
  if(buf(1:4) == '%mem') then
   read(buf(k+1:j-2),*) mem
   select case(buf(j-1:j))
   case('gb')
    mem = 1000*mem ! You like 1024? I prefer 1000
   case('mb')
   case('gw')
    mem = 1000*8*mem
   case('mw')
    mem = 8*mem
   case default
    write(iout,'(/,A)') 'ERROR in subroutine read_mem_and_nproc_from_gjf: memory&
                       & unit cannot be recognized.'
    write(iout,'(A)') "Only 'GB', 'MB', 'GW', and 'MW' are accepted."
    write(iout,'(A)') 'unit = '//TRIM(buf(j-1:j))
    stop
   end select
  else if(buf(1:6) == '%nproc') then
   read(buf(k+1:),*) np
  end if
 end do ! for while

 close(fid)

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine read_mem_and_nproc_from_gjf: incomplete&
                   & file '//TRIM(gjfname)
  stop
 end if
 return
end subroutine read_mem_and_nproc_from_gjf

! read the method and basis set from a string
! Note: please call subroutine lower before calling this subroutine,
!       in order to transform all letters to lower case
subroutine read_method_and_basis_from_buf(buf, method, basis, wfn_type)
 use print_id, only: iout
 implicit none
 integer :: i, j
 integer, intent(out) :: wfn_type ! 0/1/2/3 for undetermined/RHF/ROHF/UHF
 character(len=1200), intent(in) :: buf
 character(len=9), intent(out) :: method
 character(len=21), intent(out) :: basis

 j = index(buf, '/')
 if(j == 0) then
  write(iout,'(/,A)') "ERROR in subroutine read_method_and_basis_from_buf: no&
                     & '/' symbol found in '#' line."
  write(iout,'(A)') "Failed to identify the method name. This utility does not&
                   & support syntax like 'M062X cc-pVDZ'."
  write(iout,'(A)') "You must use the '/' symbol like 'M062X/cc-pVDZ'."
  stop
 end if

 i = index(buf(1:j-1), ' ', back=.true.)
 method = ' '
 method = buf(i+1:j-1)

 if(method(1:1) == 'u') then
  wfn_type = 3 ! UHF
  method = method(2:)
 else if(method(1:2) == 'ro') then
  wfn_type = 2 ! ROHF
  method = method(3:)
 else if(method(1:1) == 'r') then
  wfn_type = 1 ! RHF
  method = method(2:)
 else
  wfn_type = 0 ! undetermined
 end if

 i = index(buf(j+1:), ' ')
 basis = ' '
 basis = buf(j+1:j+i-1)

 return
end subroutine read_method_and_basis_from_buf

! read nfrag (the number of fragments) from a string
! Note: please call subroutine lower before calling this subroutine,
!       in order to transform all letters to lower case
subroutine read_nfrag_from_buf(buf, nfrag)
 use print_id, only: iout
 implicit none
 integer :: i, j, k
 integer, intent(out) :: nfrag
 character(len=1200), intent(in) :: buf

 i = index(buf, 'guess')
 if(i == 0) then
  write(iout,'(/,A)') "ERROR in subroutine read_nfrag_from_buf: keyword 'guess'&
                     & not found in buf."
  write(iout,'(A)') 'buf = '//TRIM(buf)
  write(iout,'(/,A)') "You should write 'guess(fragment=N)' in route section&
                     & to specify the number of fragments."
  stop
 end if

 i = i + 4
 nfrag = 0
 j = i+1 + index(buf(i+2:), '=')

 select case(buf(i+1:i+1))
 case('=')
  k = i+1 + index(buf(i+2:), ' ')
 case('(')
  k = i+1 + index(buf(i+2:), ')')
 case default
  write(iout,'(A)') 'ERROR in subroutine read_nfrag_from_buf: syntax error.'
  write(iout,'(A)') 'buf = '//TRIM(buf)
  stop
 end select
 read(buf(j+1:k-1),*) nfrag

 if(nfrag == 1) then
  write(iout,'(A)') 'ERROR in subroutine read_nfrag_from_buf: nfrag = 1.'
  write(iout,'(A)') 'Only one fragment. No need to use fragment-guess.'
  stop
 else if(nfrag < 1) then
  write(iout,'(A,I0)') 'ERROR in subroutine read_nfrag_from_buf: invalid&
                      & nfrag = ', nfrag
  stop
 end if

 return
end subroutine read_nfrag_from_buf

! generate a .gjf file from a type frag
subroutine gen_gjf_from_type_frag(frag0, guess_read, basname)
 use frag_info, only: frag
 use print_id, only: iout
 use theory_level, only: mem, nproc, method, basis, scrf, sph, eda_type, disp_type
 implicit none
 integer :: i, k(3),fid
 character(len=9) :: method0
 character(len=240), intent(in) :: basname
 type(frag), intent(in) :: frag0
 logical, intent(in) :: guess_read

 k = [index(method,'o3lyp'),index(method,'ohse2pbe'),index(method,'ohse1pbe')]

 if( ANY(k > 0) ) then
  write(iout,'(/,A)') 'ERROR in subroutine gen_gjf_from_type_frag: O3LYP,&
                     & OHSE2PBE, and OHSE1PBE functionals'
  write(iout,'(A)') 'are not supported. Because these functional names are confusing.'
  stop
 end if

 method0 = ' '
 select case(frag0%wfn_type)
 case(1)
  method0 = 'R'//TRIM(method)
 case(2)
  method0 = 'RO'//TRIM(method)
 case(3)
  method0 = 'U'//TRIM(method)
 case default
  write(iout,'(A,I0)') 'ERROR in subroutine gen_gjf_from_type_frag: invalid&
                      & wfn_type=', frag0%wfn_type
  stop
 end select

 if(eda_type==1 .and. TRIM(method0)/='Rhf') then
  write(iout,'(/,A)') 'ERROR in subroutine gen_gjf_from_type_frag: Morokuma-EDA&
                     & can only be applied to RHF.'
  write(iout,'(A)') 'But found method='//TRIM(method0)
  stop
 end if

 open(newunit=fid,file=TRIM(frag0%fname),status='replace')
 i = index(frag0%fname, '.gjf', back=.true.)
 write(fid,'(A)') '%chk='//frag0%fname(1:i-1)//'.chk'
 write(fid,'(A,I0,A)') '%mem=', mem, 'MB'
 write(fid,'(A,I0)') '%nprocshared=', nproc
 write(fid,'(A)',advance='no') '#p '//TRIM(method0)//'/'//TRIM(basis)//' nosymm'

 call lower(method0)
 select case(TRIM(method0))
 case('rhf','rohf','uhf')
  write(fid,'(A)',advance='no') ' int=nobasistransform'
 case default
  write(fid,'(A)',advance='no') ' int(nobasistransform,ultrafine)'
 end select

 select case(frag0%wfn_type)
 case(1) ! RHF
  write(fid,'(A)',advance='no') ' scf(xqc,maxcycle=300)'
 case(2) ! ROHF
  write(fid,'(A)',advance='no') ' scf(maxcycle=300,NoIncFock,NoVarAcc)'
 case(3) ! UHF
  write(fid,'(A)',advance='no') ' scf(xqc,maxcycle=300,NoIncFock,NoVarAcc,conver=7)'
 case default
  write(iout,'(A)') 'ERROR in subroutine gen_gjf_from_type_frag: wfn_type&
                   & out of range.'
  write(iout,'(A,I0)') 'frag0%wfn_type=', frag0%wfn_type
  stop
 end select

 if(LEN_TRIM(scrf) > 0) write(fid,'(A)',advance='no') ' '//TRIM(scrf)

 if(sph) then
  write(fid,'(A)',advance='no') ' 5D 7F'
 else
  write(fid,'(A)',advance='no') ' 6D 10F'
 end if

 if(disp_type == 1) then
  write(fid,'(A)',advance='no') ' em=GD3'
 else if(disp_type == 2) then
  write(fid,'(A)',advance='no') ' em=GD3BJ'
 end if

 if(guess_read) then
  write(fid,'(A)') ' guess=read'
 else
  if(frag0%noiter) then
   write(fid,'(A)') ' guess(only,save)'
  else
   if(frag0%wfn_type==3 .and. frag0%mult==1) then
    write(fid,'(A)') ' guess=mix'
   else
    write(fid,'(/)',advance='no')
   end if
  end if
 end if

 write(fid,'(/,A,/)') 'auto-generated file by frag_guess_wfn of MOKIT'
 write(fid,'(I0,1X,I0)') frag0%charge, frag0%mult

 do i = 1, frag0%natom, 1
  write(fid,'(A)',advance='no') TRIM(frag0%elem(i))
  if(frag0%ghost(i)) write(fid,'(A3)',advance='no') '-Bq'
  write(fid,'(3(1X,F18.10))') frag0%coor(1:3,i)
 end do ! for i

 write(fid,'(/)',advance='no')
 close(fid)

 call copy_gen_basis_bas2gjf(basname, frag0%fname)

 if(ANY(frag0%ghost .eqv. .true.) .and. index(basis,'genecp')>0) &
  call del_ecp_of_ghost_in_gjf(frag0%fname)

 if((.not.guess_read) .and. (.not.frag0%noiter) .and. frag0%wfn_type==3) then
  open(newunit=fid,file=TRIM(frag0%fname),status='old',position='append')
  write(fid,'(/,A)') '--Link1--'
  i = index(frag0%fname, '.gjf', back=.true.)
  write(fid,'(A)') '%chk='//frag0%fname(1:i-1)//'.chk'
  write(fid,'(A,I0,A)') '%mem=', mem, 'MB'
  write(fid,'(A,I0)') '%nprocshared=', nproc
  write(fid,'(A)',advance='no') '#p '//TRIM(method0)//' chkbasis guess=read&
   & geom=allcheck nosymm stable=opt scf(xqc,maxcycle=300,NoIncFock,NoVarAcc)'

  if(LEN_TRIM(scrf) > 0) write(fid,'(A)',advance='no') ' '//TRIM(scrf)
  if(disp_type == 1) then
   write(fid,'(A)',advance='no') ' em=GD3'
  else if(disp_type == 2) then
   write(fid,'(A)',advance='no') ' em=GD3BJ'
  end if

  select case(TRIM(method0))
  case('rhf','rohf','uhf')
   write(fid,'(A,/)') ' int=nobasistransform'
  case default
   write(fid,'(A,/)') ' int(nobasistransform,ultrafine)'
  end select
  close(fid)
 end if

 return
end subroutine gen_gjf_from_type_frag

! generate extended basis set gjf files for a fragment
subroutine gen_extend_bas_frag(frag_i, frag_n, frag_k)
 use frag_info, only: frag
 implicit none
 integer :: i, natom
 type(frag), intent(in) :: frag_i, frag_n
 type(frag), intent(out) :: frag_k

 frag_k%charge   = frag_i%charge
 frag_k%mult     = frag_i%mult
 frag_k%wfn_type = frag_i%wfn_type
 frag_k%coor     = frag_n%coor
 frag_k%elem     = frag_n%elem

 natom = frag_n%natom
 frag_k%natom = natom
 allocate(frag_k%ghost(natom))
 frag_k%ghost = .true.

 do i = 1, natom, 1
  if( ANY(frag_i%atm_map == i) ) frag_k%ghost(i) = .false.
 end do ! for i
 return
end subroutine gen_extend_bas_frag

! generate .fch files from .chk files, and create GAMESS .inp file (if EDA requested)
subroutine gen_inp_of_frags()
 use print_id, only: iout
 use frag_info, only: nfrag0, nfrag, frags
 use util_wrapper, only: fch2psi_wrap, fch2inp_wrap
 use theory_level, only: mem, wfn_type, eda_type, scrf
 use periodic_table, only: read_elem_from_gjf, elem2vdw_radii
 implicit none
 integer :: i, k
 integer :: natom  ! No. of atoms
 integer :: na, nb ! No. of alpha/beta electrons
 real(kind=8), allocatable :: radii(:)
 character(len=2), allocatable :: elem(:)
 character(len=240) :: buf, fchname, fileA, fileB
 character(len=240) :: inpname1, inpname2
 logical :: r2u
 logical, allocatable :: ghost(:)

 if(eda_type == 0) return
 k = index(frags(1)%fname, '.gjf', back=.true.)
 ! do not change k below in this subroutine

 do i = 1, nfrag, 1
  fchname = frags(i)%fname(1:k-1)//'.fch'
  if(eda_type == 4) then ! SAPT in PSI4
   call fch2psi_wrap(fchname)
   inpname1 = frags(i)%fname(1:k-1)//'_psi.inp'
   if(i < nfrag) call delete_files(2, [fchname, inpname1])
  else
   call fch2inp_wrap(fchname, .false., 0, 0)
   if(eda_type==1 .and. i==nfrag) call delete_file(fchname)
  end if
 end do ! for i

 i = index(frags(nfrag)%fname, '-')
 if(eda_type == 4) then ! SAPT in PSI4
  inpname1 = frags(nfrag)%fname(1:k-1)//'_psi.inp'
 else
  inpname1 = frags(nfrag)%fname(1:k-1)//'.inp'
 end if
 inpname2 = frags(nfrag)%fname(i+1:k-1)//'.inp'

 allocate(radii(natom), source=0d0)
 if(LEN_TRIM(scrf) > 0) then
  if(index(scrf,'smd') == 0) then
   buf = frags(nfrag)%fname(i+1:k+3)
   call read_natom_from_gjf(buf, natom)
   allocate(elem(natom), ghost(natom))
   call read_elem_from_gjf(buf, natom, elem, ghost)
   deallocate(ghost)
   forall(i = 1:natom) radii(i) = elem2vdw_radii(elem(i))
   deallocate(elem)
  end if
 end if

 if(eda_type == 4) then ! SAPT in PSI4
  do i = 1, nfrag0, 1
   if(wfn_type==3 .and. frags(i)%wfn_type/=3) then ! RHF copied to UHF
    fileA = frags(i)%fname(1:k-1)//'.A'
    fileB = frags(i)%fname(1:k-1)//'.B'
    call copy_file(fileA, fileB, .false.)
   end if
  end do ! for i

  inpname1 = frags(nfrag)%fname(1:k-1)//'_psi.inp'
  call copy_and_modify_psi4_sapt_file(inpname1, inpname2)
 else ! not SAPT

  call copy_and_modify_gamess_eda_file(natom, radii, inpname1, inpname2)
  if(eda_type == 1) call delete_file(inpname1)
  if(eda_type == 3) then ! GKS-EDA
   call copy_vec_to_append_another_inp(inpname1,inpname2,0,.false.,.true.,.true.,.false.)
  end if

  do i = 1, nfrag0, 1
   r2u = .false.
   if(wfn_type==3 .and. frags(i)%wfn_type/=3) r2u = .true.
   inpname1 = frags(i)%fname(1:k-1)//'.inp'
   call copy_vec_to_append_another_inp(inpname1,inpname2,i,.false.,.true.,.true.,r2u)

   if(eda_type == 3) then ! GKS-EDA
    inpname1 = frags(i+nfrag0)%fname(1:k-1)//'.inp'
    call copy_vec_to_append_another_inp(inpname1,inpname2,i,.true.,.true.,.true.,r2u)
   end if
  end do ! for i
 end if

 deallocate(radii, frags)
 write(iout,'(/,A)',advance='no') 'Done. Now you can submit file '//TRIM(inpname2)//' to '
 if(eda_type == 4) then ! SAPT in PSI4
  write(iout,'(A)') "PSI4 (Don't delete *.A and *.B files)."
 else
  write(iout,'(A)') 'GAMESS.'
 end if

 return
end subroutine gen_inp_of_frags

! copy content of a provided .inp file and modify it to be SAPT job
subroutine copy_and_modify_psi4_sapt_file(inpname1, inpname2)
 use print_id, only: iout
 use frag_info, only: nfrag0, nfrag, frags
 use theory_level, only: mem, basis, wfn_type
 implicit none
 integer :: i, j, natom, fid1, fid2
 character(len=240) :: buf, fileA, fileB
 character(len=240), intent(in) :: inpname1, inpname2

 open(newunit=fid1,file=TRIM(inpname1),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname2),status='replace')
 write(fid2,'(A)') '# auto-generated file by utility frag_guess_wfn of MOKIT'
 write(fid2,'(A,I0,A)') 'memory ',mem,' MB'
 write(fid2,'(/,A)') 'molecule mymol {'
 write(fid2,'(A)') 'symmetry c1'
 write(fid2,'(A)') 'no_reorient'

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(1:11) == 'no_reorient') exit
 end do ! for while
 read(fid1,'(A)') buf

 do i = 1, nfrag0, 1
  if(i > 1) write(fid2,'(A)') '--'
  write(fid2,'(I0,1X,I0)') frags(i)%charge, frags(i)%mult
  natom = COUNT(frags(i)%ghost .eqv. .false.)
  do j = 1, natom, 1
   read(fid1,'(A)') buf
   write(fid2,'(A)') TRIM(buf)
  end do ! for j
 end do ! for i
 read(fid1,'(A)') buf
 write(fid2,'(A)') '}'

 do while(.true.)
  read(fid1,'(A)') buf
  write(fid2,'(A)') TRIM(buf)
  if(buf(1:1) == '}') exit
 end do ! for while
 close(fid1, status='delete')

 write(fid2,'(/,A)') 'dimer = psi4.get_active_molecule()'
 write(fid2,'(/,A)') 'set {'
 write(fid2,'(A)') ' scf_type pk'
 write(fid2,'(A)') ' s_tolerance 1e-6'
 write(fid2,'(A)') ' e_convergence 1e5'
 write(fid2,'(A)') ' d_convergence 1e5'
 i = frags(3)%wfn_type
 select case(i)
 case(1)
  write(fid2,'(A)') ' reference rhf'
 case(2)
  write(fid2,'(A)') ' reference rohf'
 case(3)
  write(fid2,'(A)') ' reference uhf'
 case default
  write(iout,'(A)') 'ERROR in subroutine copy_and_modify_psi4_sapt_file: wfn_ty&
                    &pe out of range.'
  write(iout,'(A,I0)') 'wfn_type=', i
  stop
 end select
 write(fid2,'(A)') ' df_basis_sapt '//TRIM(basis)//'-ri'
 write(fid2,'(A)') '}'

 write(fid2,'(/,A)') 'set df_ints_io save'
 write(fid2,'(A)') "psi4.IO.set_default_namespace('dimer')"
 write(fid2,'(A)') "Edim, wfn_dimer = energy('scf', molecule=dimer, return_wfn=True)"
 write(fid2,'(A)') 'set df_ints_io load'
 write(fid2,'(A)') '# the above scf makes every array allocated'
 i = index(inpname2, '.inp', back=.true.)
 write(fileA,'(I3.3,A)') nfrag, '-'//inpname2(1:i-1)//'.A'
 write(fileB,'(I3.3,A)') nfrag, '-'//inpname2(1:i-1)//'.B'
 write(fid2,'(A)') "wfn_dimer.Ca().load('"//TRIM(fileA)//"')"
 if(wfn_type == 3) write(fid2,'(A)') "wfn_dimer.Cb().load('"//TRIM(fileB)//"')"
 write(fid2,'(A)') 'wfn_dimer.to_file(wfn_dimer.get_scratch_filename(180))'
 write(fid2,'(A)') 'set {'
 write(fid2,'(A)') ' guess read'
 write(fid2,'(A)') ' e_convergence 1e-8'
 write(fid2,'(A)') ' d_convergence 1e-6'
 write(fid2,'(A)') '}'
 write(fid2,'(A)') "Edim, wfn_dimer = energy('scf', molecule=dimer, return_wfn=True)"

 write(fid2,'(/,A)') 'monomerA = dimer.extract_subsets(1,2)'
 write(fid2,'(A)') "psi4.IO.change_file_namespace(97, 'dimer', 'monomerA')"
 write(fid2,'(A)') "psi4.IO.set_default_namespace('monomerA')"
 write(fid2,'(A)') 'set {'
 write(fid2,'(A)') ' e_convergence 1e5'
 write(fid2,'(A)') ' d_convergence 1e5'
 write(fid2,'(A)') '}'
 write(fid2,'(A)') "EmonA, wfn_monA = energy('scf', molecule=monomerA, return_wfn=True)"
 write(fid2,'(A)') '# the above scf makes every array allocated'
 i = index(inpname2, '.inp', back=.true.)
 write(fileA,'(A)') '001-'//inpname2(1:i-1)//'.A'
 write(fileB,'(A)') '001-'//inpname2(1:i-1)//'.B'
 write(fid2,'(A)') "wfn_monA.Ca().load('"//TRIM(fileA)//"')"
 if(wfn_type == 3) write(fid2,'(A)') "wfn_monA.Cb().load('"//TRIM(fileB)//"')"
 write(fid2,'(A)') 'wfn_monA.to_file(wfn_monA.get_scratch_filename(180))'
 write(fid2,'(A)') 'set {'
 write(fid2,'(A)') ' guess read'
 write(fid2,'(A)') ' e_convergence 1e-8'
 write(fid2,'(A)') ' d_convergence 1e-6'
 write(fid2,'(A)') '}'
 write(fid2,'(A)') "EmonA, wfn_monA = energy('scf', molecule=monomerA, return_wfn=True)"

 write(fid2,'(/,A)') 'monomerB = dimer.extract_subsets(2,1)'
 write(fid2,'(A)') "psi4.IO.change_file_namespace(97, 'monomerA', 'monomerB')"
 write(fid2,'(A)') "psi4.IO.set_default_namespace('monomerB')"
 write(fid2,'(A)') 'set {'
 write(fid2,'(A)') ' e_convergence 1e5'
 write(fid2,'(A)') ' d_convergence 1e5'
 write(fid2,'(A)') '}'
 write(fid2,'(A)') "EmonB, wfn_monB = energy('scf', molecule=monomerB, return_wfn=True)"
 write(fid2,'(A)') '# the above scf makes every array allocated'
 i = index(inpname2, '.inp', back=.true.)
 write(fileA,'(A)') '002-'//inpname2(1:i-1)//'.A'
 write(fileB,'(A)') '002-'//inpname2(1:i-1)//'.B'
 write(fid2,'(A)') "wfn_monB.Ca().load('"//TRIM(fileA)//"')"
 if(wfn_type == 3) write(fid2,'(A)') "wfn_monB.Cb().load('"//TRIM(fileB)//"')"
 write(fid2,'(A)') 'wfn_monB.to_file(wfn_monB.get_scratch_filename(180))'
 write(fid2,'(A)') 'set {'
 write(fid2,'(A)') ' guess read'
 write(fid2,'(A)') ' e_convergence 1e-8'
 write(fid2,'(A)') ' d_convergence 1e-6'
 write(fid2,'(A)') '}'
 write(fid2,'(A)') "EmonB, wfn_monB = energy('scf', molecule=monomerB, return_wfn=True)"

 write(fid2,'(/,A)') "psi4.IO.change_file_namespace(97, 'monomerB', 'dimer')"
 write(fid2,'(A)') "psi4.IO.set_default_namespace('dimer')"

 write(fid2,'(/,A)') "aux_basis = psi4.core.BasisSet.build(wfn_dimer.molecule(),&
                   & ""DF_BASIS_SAPT"","
 write(fid2,'(A)') " psi4.core.get_global_option(""DF_BASIS_SAPT""),"
 write(fid2,'(A)') " ""RIFIT"", psi4.core.get_global_option(""BASIS""))"
 write(fid2,'(A)') "wfn_dimer.set_basisset(""DF_BASIS_SAPT"", aux_basis)"
 write(fid2,'(A)') "wfn_dimer.set_basisset(""DF_BASIS_ELST"", aux_basis)"

 write(fid2,'(/,A)') 'psi4.sapt(wfn_dimer, wfn_monA, wfn_monB)'
 close(fid2)
 return
end subroutine copy_and_modify_psi4_sapt_file

! copy content of a provided .inp file and modify it to be EDA job
subroutine copy_and_modify_gamess_eda_file(natom, radii, inpname1, inpname2)
 use print_id, only: iout
 use frag_info, only: nfrag0, frags
 use theory_level, only: method, eda_type, scrf, solvent, disp_type
 implicit none
 integer :: i, j, m, fid1, fid2
 integer, intent(in) :: natom
 real(kind=8), intent(in) :: radii(natom)
 character(len=9) :: dft_in_gms
 character(len=240) :: buf1, buf2
 character(len=240), intent(in) :: inpname1, inpname2

 open(newunit=fid1,file=TRIM(inpname1),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname2),status='replace')

 read(fid1,'(A)') buf1
 i = index(buf1, 'RUNTYP=ENERGY')
 buf2 = buf1(1:i+6)//'EDA'//TRIM(buf1(i+13:))
 write(fid2,'(A)') TRIM(buf2)

 call convert_dft_name_gau2gms(method, dft_in_gms)
 read(fid1,'(A)') buf1

 i = index(buf1, '$END')
 if(i == 0) then
  write(iout,'(A)') "ERROR in subroutine gen_inp_of_frags: no '$END' found&
                   & in the 2nd line of file "//TRIM(inpname1)
  close(fid1)
  close(fid2,status='delete')
  stop
 end if

 write(fid2,'(A)') buf1(1:i-1)//'DFTTYP='//TRIM(dft_in_gms)//' $END'
 if(TRIM(dft_in_gms) /= 'NONE') then
  write(fid2,'(A)',advance='no') ' $DFT NRAD0=99 NLEB0=590 NRAD=99 NLEB=590'
  if(disp_type == 2) then
   write(fid2,'(A)',advance='no') ' IDCVER=3'
  else
   write(fid2,'(A)',advance='no') ' IDCVER=4'
  end if
  write(fid2,'(A)') ' $END'
 end if

 do while(.true.)
  read(fid1,'(A)') buf1
  if(buf1(2:7) == '$GUESS') exit
  if(buf1(2:5) == '$SCF') then
   if(eda_type == 1) then
    buf1 = ' $SCF DIRSCF=.F. $END'
   else
    if(TRIM(dft_in_gms) /= 'NONE') buf1 = ' $SCF DIRSCF=.T. DIIS=.F. SOSCF=.T. $END'
   end if
  end if
  write(fid2,'(A)') TRIM(buf1)
 end do ! for while

 select case(eda_type)
 case(1) ! MOROKUMA-EDA
  write(fid2,'(A)',advance='no') ' $MOROKM IATM(1)='
  do i = 1, nfrag0, 1
   if(i > 1) write(fid2,'(A)',advance='no') ','
   write(fid2,'(I0)',advance='no') frags(i)%natom
  end do ! for i

  write(fid2,'(A)',advance='no') ' ICHM(1)='
  do i = 1, nfrag0, 1
   if(i > 1) write(fid2,'(A)',advance='no') ','
   write(fid2,'(I0)',advance='no') frags(i)%charge
  end do ! for i

  write(fid2,'(A)') ' $END'

 case(2) ! LMO-EDA
  write(iout,'(A)') 'ERROR in subroutine gen_inp_of_frags: LMO-EDA not&
                   & supported surrently.'
  stop

 case(3) ! GKS-EDA
  write(fid2,'(A)',advance='no') ' $LMOEDA MATOM(1)='
  do i = 1, nfrag0, 1
   if(i > 1) write(fid2,'(A)',advance='no') ','
   write(fid2,'(I0)',advance='no') frags(i)%natom
  end do ! for i

  write(fid2,'(A)',advance='no') ' MCHARG(1)='
  do i = 1, nfrag0, 1
   if(i > 1) write(fid2,'(A)',advance='no') ','
   write(fid2,'(I0)',advance='no') frags(i)%charge
  end do ! for i

  write(fid2,'(A)',advance='no') ' MMULT(1)='
  do i = 1, nfrag0, 1
   if(i > 1) write(fid2,'(A)',advance='no') ','
   if(frags(i)%pos) then
    write(fid2,'(I0)',advance='no') frags(i)%mult
   else
    write(fid2,'(I0)',advance='no') -frags(i)%mult
   end if
  end do ! for i

  write(fid2,'(/,A)') '  EDATYP=GKS RDVECM=.T. $END'
  write(fid2,'(A)') ' $GUESS GUESS=HCORE $END'

 case default
  write(iout,'(A,I0)') 'ERROR in subroutine gen_inp_of_frags: invalid eda_type=',&
                       eda_type
  close(fid1)
  close(fid2)
  stop
 end select

 if(LEN_TRIM(scrf) > 0) then
  write(fid2,'(A)',advance='no') ' $PCM IEF=-3 SOLVNT='//TRIM(solvent)
  if(index(scrf,'smd') > 0) then
   write(fid2,'(A)',advance='no') ' SMD=.T. $END'
   write(iout,'(/,A)') 'Warning from subroutine gen_inp_of_frags: scrf=SMD detected.'
   write(iout,'(A)') 'If you encounter SCF convergence problems later in GAMESS,&
                    & you can change'
   write(iout,'(A)') 'to scrf=PCM or remove scrf=SMD (i.e. gas phase computation).'
  else
   write(fid2,'(A,/,A)') ' $END',' $PCMCAV ALPHA(1)=1.1'
   j = natom/8
   do i = 1, j, 1
    write(fid2,'(A,I0,A,F6.4,7(A1,F6.4))') '  RIN(',8*i-7,')=',radii(8*i-7),(',',radii(m),m=8*i-6,8*i)
   end do ! for i
   if(natom-8*j > 0) then
    write(fid2,'(A,I0,A,F6.4,7(A1,F6.4))') '  RIN(',8*j+1,')=',radii(8*j+1),(',',radii(m),m=8*j+2,natom)
   end if
   write(fid2,'(A)') ' $END'
  end if
 end if

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf1
  if(i /= 0) exit
  if(buf1(2:5) == '$VEC') exit
  write(fid2,'(A)') TRIM(buf1)
 end do ! for while

 close(fid1)
 close(fid2)
 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine gen_inp_of_frags: no '$VEC' found&
                   & in file "//TRIM(inpname1)
  stop
 end if

 return
end subroutine copy_and_modify_gamess_eda_file

! copy $VEC from file inpname1 into inpname2, by appending
subroutine copy_vec_to_append_another_inp(inpname1, inpname2, ivec, extended, &
                                          deleted, occ, r2u)
 use print_id, only: iout
 implicit none
 integer :: i, nline, fid1, fid2
 integer :: na, nb, nif, nbf, nbf1
 ! Note: nbf1>=nbf always holds. nbf1 is always the number of basis functions
 !  of Cartesian functions, but nbf depends on the keywords(ISPHER=1) used
 integer, intent(in) :: ivec ! ivec>=-1
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname1, inpname2
 logical :: uhf
 logical, intent(in) :: extended ! $VEC or $VEA
 logical, intent(in) :: deleted  ! whether to delete inpname1
 logical, intent(in) :: occ      ! whether to copy only occupied orbitals or all orbitals
 logical, intent(in) :: r2u      ! whether to copy nb orbitals from Alpha MOs

 call check_uhf_in_gms_inp(inpname1, uhf)
 if(uhf .and. r2u) then
  write(iout,'(A)') "ERROR in subroutine copy_vec_to_append_another_inp: both&
                   & logical variables 'uhf' and"
  write(iout,'(A)') "'r2u' are .True. This is not allowed."
  write(iout,'(A)') 'inpname1='//TRIM(inpname1)//', inpname2='//TRIM(inpname2)
  stop
 end if

 if(occ) then
  call read_na_nb_nif_nbf_from_gms_inp(inpname1, na, nb, nif, nbf)
  call read_nbf_from_dat(inpname1, nbf1)
!  write(*,'(5(A,I0))') 'na=',na,',nb=',nb,',nif=',nif,',nbf=',nbf,',nbf1=',nbf1
 end if

 open(newunit=fid1,file=TRIM(inpname1),status='old',position='rewind')
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:5) == '$VEC') exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine copy_vec_to_append_another_inp:'
  write(iout,'(A)') "No '$VEC' found in file "//TRIM(inpname1)
  close(fid1)
  stop
 end if

 open(newunit=fid2,file=TRIM(inpname2),status='old',position='append')
 if(extended) then
  if(ivec > 0) then
   write(fid2,'(A,I0)') ' $VEA',ivec
  else
   write(iout,'(A)') 'ERROR in subroutine copy_vec_to_append_another_inp:'
   write(iout,'(A)') 'ivec<=0 when extended=.True. This is not allowed.'
   close(fid1)
   close(fid2)
   stop
  end if
 else ! not extended basis set
  if(ivec > 0) then
   write(fid2,'(A,I0)') ' $VEC',ivec
  else if(ivec == 0) then
   write(fid2,'(A)') ' $VEC0'
  else
   write(fid2,'(A)') ' $VEC'
  end if
 end if

 if(occ) then ! only copy occupied MOs of a fragment
  nline = nbf1/5
  if(nbf1-5*nline > 0) nline = nline + 1

  do i = 1, nline*na, 1
   read(fid1,'(A)') buf
   write(fid2,'(A)') TRIM(buf)
  end do ! for i

  if(uhf) then
   do i = 1, nline*(nif-na), 1
    read(fid1,'(A)') buf
   end do ! for i
   do i = 1, nline*nb, 1
    read(fid1,'(A)') buf
    write(fid2,'(A)') TRIM(buf)
   end do ! for i
  end if

  if(r2u) then ! copy beta MOs from alpha MOs
   do while(.true.)
    BACKSPACE(fid1)
    BACKSPACE(fid1)
    read(fid1,'(A)') buf
    if(buf(2:5) == '$VEC') exit
   end do ! for while
   do i = 1, nline*nb, 1
    read(fid1,'(A)') buf
    write(fid2,'(A)') TRIM(buf)
   end do ! for i
  end if

 else ! copy all MOs of a fragment

  do while(.true.)
   read(fid1,'(A)') buf
   if(buf(2:5) == '$END') exit
   write(fid2,'(A)') TRIM(buf)
  end do ! for while

  if(r2u) then ! copy beta MOs from alpha MOs
   do while(.true.)
    BACKSPACE(fid1)
    BACKSPACE(fid1)
    read(fid1,'(A)') buf
    if(buf(2:5) == '$VEC') exit
   end do ! for while
   do while(.true.)
    read(fid1,'(A)') buf
    if(buf(2:5) == '$END') exit
    write(fid2,'(A)') TRIM(buf)
   end do ! for while
  end if
 end if

 write(fid2,'(A)') ' $END'
 if(deleted) then
  close(fid1,status='delete')
 else
  close(fid1)
 end if

 close(fid2)
 return
end subroutine copy_vec_to_append_another_inp

! read the scrf string from a longbuf
subroutine read_scrf_string_from_buf(longbuf, scrf)
 implicit none
 integer :: i, j
 character(len=1200), intent(in) :: longbuf
 character(len=40), intent(out) :: scrf

 scrf = ' '
 i = index(longbuf, 'scrf')
 if(i == 0) return

 j = index(longbuf(i:), ' ')
 if(j == 0) then
  j = i + 3
 else
  j = j + i - 2
 end if

 scrf = longbuf(i:j)
 return
end subroutine read_scrf_string_from_buf

! determine implicit solvent model from the scrf string
! Note: the input scrf must be in lower case
subroutine determine_solvent_from_gau2gms(scrf, solvent)
 use print_id, only: iout
 implicit none
 integer :: i, j, k(3)
 character(len=20) :: solvent_gau
 character(len=40), intent(in) :: scrf
 character(len=20), intent(out) :: solvent

 solvent = ' '
 if(LEN_TRIM(scrf) == 0) return

 i = index(scrf,'solvent')
 if(i == 0) then
  solvent = 'H2O'
  return
 end if

 k = [index(scrf(i+8:),')'),index(scrf(i+8:),','),index(scrf(i+8:),' ')]
 forall(j=1:3, k(j)==0) k(j) = LEN_TRIM(scrf)
 j = MINVAL(k)
 solvent_gau = ' '
 solvent_gau = scrf(i+8:i+6+j)

 select case(TRIM(solvent_gau))
 case('water','h2o')
  solvent = 'H2O'
 case('methanol','ch3oh')
  solvent = 'CH3OH'
 case('ethanol','ch3ch2oh')
  solvent = 'C2H5OH'
 case('chloroform')
  solvent = 'CHCl3'
 case('carbontetrachloride','ccl4')
  solvent = 'CCl4'
 case('dichloromethane','methylenechloride','ch2cl2')
  solvent = 'CH2Cl2'
 case('dichloroethane','1,2-dichloroethane','ch2clch2cl')
  solvent = 'C2H4Cl2'
 case('benzene','c6h6')
  solvent = 'C6H6'
 case('toluene','c6h5ch3')
  solvent = 'C6H5CH3'
 case('chlorobenzene','c6h5cl')
  solvent = 'C6H5Cl'
 case('nitromethane','ch3no2')
  solvent = 'CH3NO2'
 case('heptane','n-heptane')
  solvent = 'C7H16'
 case('cyclohexane','c6h12')
  solvent = 'C6H12'
 case('aniline','c6h5nh2')
  solvent = 'C6H5NH2'
 case('acetone','ch3coch3')
  solvent = 'CH3COCH3'
 case('tetrahydrofuran','thf')
  solvent = 'THF'
 case('dimethylsulfoxide','dmso','ch3soch3')
  solvent = 'DMSO'
 case default
  write(iout,'(A)') 'Warning in subroutine determine_solvent_from_gau2gms:&
                   & implicit solvent name in Gaussian'
  write(iout,'(A)') '.gjf file cannot be identified. You should open GAMESS&
                   & .inp file and add it by yourself.'
  write(iout,'(A)') 'Solvent in .gjf file: '//TRIM(solvent_gau)
  solvent = 'INPUT'
 end select

 return
end subroutine determine_solvent_from_gau2gms

! convert DFT name of Gaussian to that of GAMESS
subroutine convert_dft_name_gau2gms(method, dft_in_gms)
 use print_id, only: iout
 implicit none
 integer :: i
 character(len=9), intent(in) :: method
 character(len=9), intent(out) :: dft_in_gms

 dft_in_gms = ' '

 select case(TRIM(method))
 case('b3lyp','blyp','x3lyp','b3p86','b3pw91','wB97','wB97X','b98','m05','m06',&
      'm11','bmk','bp86','tpssh','revtpss')
  dft_in_gms = method
 case('m052x')
  dft_in_gms = 'M05-2X'
 case('m06hf')
  dft_in_gms = 'M06-HF'
 case('m062x')
  dft_in_gms = 'M06-2X'
 case('m06l')
  dft_in_gms = 'M06-L'
 case('m08hx')
  dft_in_gms = 'M08-HX'
 case('m11l')
  dft_in_gms = 'M11-L'
 case('wb97xd')
  dft_in_gms = 'wB97X-D'
 case('b971')
  dft_in_gms = 'B97-1'
 case('b972')
  dft_in_gms = 'B97-2'
 case('b97d')
  dft_in_gms = 'B97-D'
 case('pbepbe')
  dft_in_gms = 'PBE'
 case('pbe1pbe')
  dft_in_gms = 'PBE0'
 case('cam-b3lyp')
  dft_in_gms = 'CAMB3LYP'
 case('tpsstpss')
  dft_in_gms = 'TPSS'
 case('hf','rhf','rohf','uhf')
  dft_in_gms = 'NONE'
 case default
  write(iout,'(A)') 'Warning in subroutine convert_dft_name_gau2gms: functional name&
                   & cannot be recognized.'
  write(iout,'(A)') 'DFTTYP=NONE will be set in .inp file. You can modify it by yourself.'
  dft_in_gms = 'NONE'
 end select

 return
end subroutine convert_dft_name_gau2gms

! delete ECP/PP of ghost atoms in a given .gjf file
subroutine del_ecp_of_ghost_in_gjf(gjfname)
 use periodic_table, only: read_elem_from_gjf
 implicit none
 integer :: i, j, nbat1, nbat2, fid, fid1, RENAME
 integer :: nblank, natom
 character(len=2), allocatable :: elem(:)
 character(len=3) :: str
 character(len=6) :: str1
 character(len=240) :: buf, gjfname1
 character(len=240), intent(in) :: gjfname
 logical :: skipped
 logical, allocatable :: ghost(:) ! .True./.False. for ghost atom or not

 call read_natom_from_gjf(gjfname, natom)
 allocate(elem(natom), ghost(natom))
 call read_elem_from_gjf(gjfname, natom, elem, ghost)

 gjfname1 = TRIM(gjfname)//'.t'
 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(gjfname1),status='replace')
 nblank = 0

 do while(.true.)
  read(fid,'(A)') buf
  write(fid1,'(A)') TRIM(buf)
  if(LEN_TRIM(buf) == 0) nblank = nblank + 1
  if(nblank == 4) exit
 end do ! for while

 read(fid,'(A)') buf

 if(LEN_TRIM(buf) > 0) then
  do while(.true.)
   read(buf,*) str
   str = ADJUSTL(str)
   if(str(1:1) == '-') str = str(2:)

   j = IACHAR(str(1:1))
   if(j>96 .and. j<123) str(1:1) = ACHAR(j-32)
   j = IACHAR(str(2:2))
   if(j>64 .and. j<91) str(2:2) = ACHAR(j+32)

   skipped = .false.
   do i = 1, natom, 1
    if(elem(i) == TRIM(str)) then
     if(ghost(i)) then
      skipped = .true.
      exit
     end if
    end if
   end do ! for i

   if(.not. skipped) write(fid1,'(A)') TRIM(buf)
   read(fid,'(A)') buf
   read(buf,*,iostat=i) str1, nbat1

   if(skipped) then
    if(i == 0) then ! ECP/PP data, not name
     nbat1 = nbat1 + 1
     do i = 1, nbat1, 1
      read(fid,'(A)') buf
      read(fid,*) nbat2
      do j = 1, nbat2, 1
       read(fid,'(A)') buf
      end do ! for j
     end do ! for i
    end if
   else ! not skipped
    write(fid1,'(A)') TRIM(buf)

    if(i == 0) then ! ECP/PP data, not name
     nbat1 = nbat1 + 1
     do i = 1, nbat1, 1
      read(fid,'(A)') buf
      write(fid1,'(A)') TRIM(buf)
      read(fid,'(A)') buf
      write(fid1,'(A)') TRIM(buf)
      read(buf,*) nbat2
      do j = 1, nbat2, 1
       read(fid,'(A)') buf
       write(fid1,'(A)') TRIM(buf)
      end do ! for j
     end do ! for i
    end if
   end if

   read(fid,'(A)') buf
   if(LEN_TRIM(buf) == 0) exit
  end do ! for while
 end if

 write(fid1,'(/)',advance='no')
 close(fid,status='delete')
 close(fid1)
 deallocate(elem, ghost)
 i = RENAME(TRIM(gjfname1), TRIM(gjfname))
 return
end subroutine del_ecp_of_ghost_in_gjf

! construct supermolecule occ MOs from direct sum of fragment occ MOs
! construct supermolecule vir MOs from direct sum of fragment vir MOs
subroutine direct_sum_frag_mo2super_mo(n, gjfname0, wfn_type0, pos, gjfname, wfn_type)
 use print_id, only: iout
 use util_wrapper, only: formchk, unfchk
 implicit none
 integer :: i, j, fid0, fid, RENAME
 integer :: k1, k2, k3, k4, k5, k6, k7
 integer :: nbf0, nif0, nbf, nif
 integer :: na0, nb0, na, nb   ! No. of alpha/beta electrons
 integer, intent(in) :: n
 integer, intent(in) :: wfn_type0(n), wfn_type ! 1/2/3 for RHF/ROHF/UHF
 character(len=240) :: buf, fchname0, fchname, chkname
 character(len=240), intent(in) :: gjfname0(n), gjfname
 real(kind=8), allocatable :: mo_a0(:,:), mo_b0(:,:) ! fragment MOs
 real(kind=8), allocatable :: mo_a(:,:), mo_b(:,:)   ! supermolecule MOs
 logical :: alive
 logical, intent(in) :: pos(n)

 i = index(gjfname, '.gjf', back=.true.)
 if(i == 0) then
  write(iout,'(A)') "ERROR in subroutine direct_sum_frag_mo2super_mo: no '.gjf'&
                  & string found in file "//TRIM(gjfname)
  stop
 end if

 if(ANY(wfn_type0 == 0) .or. wfn_type==0) then
  write(iout,'(A)') 'ERROR in subroutine direct_sum_frag_mo2super_mo: there exists&
                   & at least one file'
  write(iout,'(A)') 'whose wfn_type is 0. Not allowed.'

  do i = 1, n, 1
   write(iout,'(2(A,I0))') 'i=',i,','//TRIM(gjfname0(i))//',wfn_type0(i)=',wfn_type0(i)
  end do ! for i
  write(iout,'(A,I0)') 'wfn_type=', wfn_type
  stop
 end if

 ! modify 'guess(only,save)' in gjf file into 'guess=read geom=allcheck'
 fchname = TRIM(gjfname)//'.t'
 open(newunit=fid0,file=TRIM(gjfname),status='old',position='rewind')
 open(newunit=fid,file=TRIM(fchname),status='replace')
 do while(.true.)
  read(fid0,'(A)') buf
  i = index(buf,'guess')
  if(i > 0) then
   buf = buf(1:i+4)//'=read geom=allcheck'
   i = index(buf,'/')
   j = index(buf(i+1:),' ')
   buf = buf(1:i)//'chkbasis'//buf(i+j:)

   i = index(buf,',conver=')
   if(i > 0) then
    j = index(buf(i+1:),')')
    buf = buf(1:i-1)//')'//buf(i+j+1:)
   end if

   if(wfn_type == 3) buf = TRIM(buf)//' stable=opt'
  end if
  write(fid,'(A)') TRIM(buf)
  if(LEN_TRIM(buf) == 0) exit
 end do ! for while
 close(fid0,status='delete')
 close(fid)
 i = RENAME(TRIM(fchname), TRIM(gjfname))
 ! modify done

 i = index(gjfname, '.gjf', back=.true.)
 fchname = gjfname(1:i-1)//'.fch'
 chkname = gjfname(1:i-1)//'.chk'
 inquire(file=TRIM(fchname), exist=alive)
 if(.not. alive) call formchk(chkname)
 call read_na_and_nb_from_fch(fchname, na, nb)
 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(mo_a(nbf,nif), source=0d0)
 if(wfn_type == 3) allocate(mo_b(nbf,nif), source=0d0)
 k1 = 0; k2 = 0; k3 = na; k4 = 0; k5 = 0; k6 = nb; k7 = 0

 do i = 1, n, 1
  j = index(gjfname0(i), '.gjf', back=.true.)
  fchname0 = gjfname0(i)(1:j-1)//'.fch'
  call read_nbf_and_nif_from_fch(fchname0, nbf0, nif0)
  allocate(mo_a0(nbf0,nif0))
  if(pos(i)) then
   call read_na_and_nb_from_fch(fchname0, na0, nb0)
   call read_mo_from_fch(fchname0, nbf0, nif0, 'a', mo_a0)
  else ! negative spin
   call read_na_and_nb_from_fch(fchname0, nb0, na0)
   if(wfn_type0(i) == 3) then
    call read_mo_from_fch(fchname0, nbf0, nif0, 'b', mo_a0)
   else ! ROHF
    call read_mo_from_fch(fchname0, nbf0, nif0, 'a', mo_a0)
   end if
  end if

  if(k1+nbf0 > nbf) then
   write(iout,'(A)') 'ERROR in subroutine direct_sum_frag_mo2super_mo: the 1st&
                    & dimension of array mo_a out of range!'
   write(iout,'(5(A,I0))') 'k1=',k1,',nbf0=',nbf0,',nbf=',nbf,',i=',i,',n=',n
   stop
  end if
  k4 = k3 + nif0 - na0
  if(k4 > nif) then
   write(iout,'(A)') 'ERROR in subroutine direct_sum_frag_mo2super_mo: the 2nd&
                    & dimension of array mo_a out of range!'
   write(iout,'(3(A,I0))') 'k3=',k3,',nif0=',nif0,',na0=',na0
   write(iout,'(4(A,I0))') 'k4=',k4,',nif=',nif,',i=',i,',n=',n
   stop
  end if
  mo_a(k1+1:k1+nbf0, k2+1:k2+na0) = mo_a0(:,1:na0)
  mo_a(k1+1:k1+nbf0, k3+1:k4) = mo_a0(:,na0+1:nif0)

  if(wfn_type == 3) then ! UHF
   k7 = k6 + nif0 - nb0
   if(k7 > nif) then
    write(iout,'(A)') 'ERROR in subroutine direct_sum_frag_mo2super_mo: the 2nd&
                     & dimension of array mo_b out of range!'
    write(iout,'(3(A,I0))') 'k6=',k6,',nif0=',nif0,',nb0=',nb0
    write(iout,'(4(A,I0))') 'k7=',k7,',nif=',nif,',i=',i,',n=',n
    stop
   end if
   if(wfn_type0(i) == 3) then ! UHF fragment
    allocate(mo_b0(nbf0,nif0))
    if(pos(i)) then
     call read_mo_from_fch(fchname0, nbf0, nif0, 'b', mo_b0)
    else ! negative spin
     call read_mo_from_fch(fchname0, nbf0, nif0, 'a', mo_b0)
    end if
    mo_b(k1+1:k1+nbf0, k5+1:k5+nb0) = mo_b0(:,1:nb0)
    mo_b(k1+1:k1+nbf0, k6+1:k7) = mo_b0(:,nb0+1:nif0)
    deallocate(mo_b0)
   else ! RHF/ROHF fragment
    mo_b(k1+1:k1+nbf0, k5+1:k5+nb0) = mo_a0(:,1:nb0)
    mo_b(k1+1:k1+nbf0, k6+1:k7) = mo_a0(:,nb0+1:nif0)
   end if
   k5 = k5 + nb0
   k6 = k7
  end if

  deallocate(mo_a0)
  k1 = k1 + nbf0
  k2 = k2 + na0
  k3 = k4
 end do ! for i

 call write_mo_into_fch(fchname, nbf, nif, 'a', mo_a)
 deallocate(mo_a)

 if(wfn_type == 3) then
  call write_mo_into_fch(fchname, nbf, nif, 'b', mo_b)
  deallocate(mo_b)
 end if

 call unfchk(fchname, chkname)
 return
end subroutine direct_sum_frag_mo2super_mo

! sum fragment Total SCF Densities and print the 
subroutine sum_frag_density_and_prt_into_fch(n, fname0, fname)
 use print_id, only: iout
 use util_wrapper, only: formchk
 implicit none
 integer :: i, j, nif, nbf, nbf1
 integer, intent(in) :: n
 character(len=240) :: chkname
 character(len=240), intent(in) :: fname0(n), fname
 character(len=240), allocatable :: fchname(:)
 real(kind=8), allocatable :: dm(:,:), dm1(:,:)
 logical :: alive
 logical, allocatable :: has_spin_density(:)

 allocate(fchname(n+1))

 ! if fname0(n) and fname are .gjf files, convert filenames into .fch
 do i = 1, n, 1
  j = index(fname0(i), '.fch', back=.true.)
  if(j == 0) then
   j = index(fname0(i), '.gjf', back=.true.)
   fchname(i) = fname0(i)(1:j-1)//'.fch'
  else ! j > 0
   fchname(i) = fname0(i)
  end if
 end do ! for i

 j = index(fname, '.fch', back=.true.)
 if(j == 0) then
  j = index(fname, '.gjf', back=.true.)
  fchname(n+1) = fname(1:j-1)//'.fch'
 else ! j > 0
  fchname(n+1) = fname
  j = index(fname, '.fch', back=.true.)
 end if
 chkname = fname(1:j-1)//'.chk'
 inquire(file=TRIM(fchname(n+1)), exist=alive)
 if(.not. alive) call formchk(chkname, fchname(n+1))

 call read_nbf_and_nif_from_fch(fchname(n+1), nbf, nif)
 allocate(dm(nbf,nbf), dm1(nbf,nbf))
 dm = 0d0

 do i = 1, n, 1
  call read_nbf_and_nif_from_fch(fchname(i), nbf1, nif)
  if(nbf1 /= nbf) then
   write(iout,'(A)') 'ERROR in subroutine sum_frag_density_and_prt_into_fch:&
                    & nbf1 /= nbf.'
   write(iout,'(A,I0,A)') 'Inconsistent nbf between fragment ',i,' and the&
                         & total system.'
   write(iout,'(2(A,I0))') 'nbf1=', nbf1, ', nbf=', nbf
   stop
  end if
  call read_density_from_fch(fchname(i), 1, nbf, dm1)
  dm = dm + dm1
 end do ! for i
 call write_density_into_fch(fchname(n+1), nbf, .true., dm)

 allocate(has_spin_density(n))
 has_spin_density = .false.
 dm = 0d0

 do i = 1, n, 1
  call detect_spin_scf_density_in_fch(fchname(i), has_spin_density(i))
  if(.not. has_spin_density(i)) cycle
  call read_density_from_fch(fchname(i), 2, nbf, dm1)
  dm = dm + dm1
 end do ! for i

 call detect_spin_scf_density_in_fch(fchname(n+1), alive)

 if(ANY(has_spin_density .eqv. .true.)) then
  if(alive) then
   call write_density_into_fch(fchname(n+1), nbf, .false., dm)
  else ! no Spin SCF Density in the total system .fch file
   write(iout,'(A)') 'ERROR in subroutine sum_frag_density_and_prt_into_fch: so&
                     &me fragment has UHF-type wave'
   write(iout,'(A)') 'function, but the total system has RHF-type wave function.&
                    & Inconsistency detected.'
   stop
  end if
 else ! all fragments has RHF-type wave function
  if(alive) call write_density_into_fch(fchname(n+1), nbf, .false., dm)
 end if

 deallocate(fchname, dm, dm1, has_spin_density)
 return
end subroutine sum_frag_density_and_prt_into_fch

! modify 'guess(only,save)' in a given .gjf file
subroutine modify_guess_only_in_gjf(gjfname)
 implicit none
 integer :: i, j, fid1, fid2, RENAME
 character(len=240) :: buf, fname
 character(len=240), intent(in) :: gjfname

 fname = TRIM(gjfname)//'.t'
 open(newunit=fid1,file=TRIM(gjfname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(fname),status='replace')

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(1:1) == '#') exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 close(fid1,status='delete')

 i = index(buf, '/')
 j = index(buf(i+1:),' ') + i - 1
 buf = buf(1:i-1)//' chkbasis'//buf(j+1:)

 i = index(buf, 'guess')
 buf(i:) = 'guess=read geom=allcheck'

 ! assuming '#p R...' or '#p U...'
 if(buf(4:4)=='u' .or. buf(4:4)=='U') then
  buf = TRIM(buf)//' stable=opt'
  i = index(buf, ',conver=6')
  if(i == 0) i = index(buf, ',conver=7')
  if(i > 0) buf = buf(1:i-1)//buf(i+9:)
 end if

 write(fid2,'(A,/)') TRIM(buf)
 close(fid2)
 i = RENAME(TRIM(fname), TRIM(gjfname))
 return
end subroutine modify_guess_only_in_gjf

