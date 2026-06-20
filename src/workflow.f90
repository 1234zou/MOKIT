! firstly written by jxzou at 20260608
! TODO:
! 1) support specification of DFT in two steps
! 2) support specification of temperature in pNMR
! 3) check whether Gaussian g-tensor is correct and can be used alternatively

module mol_and_calc_info
 implicit none
 integer :: mem = 4      ! GB
 integer :: nproc = 1    ! CPU core(s)
 integer :: charge = 0   ! the net charge of the molecule
 integer :: mult = 1     ! spin multiplicity
 integer :: natom = 0    ! the number of atoms
 integer :: dis_type = 2 ! dispersion correction, 0/1/2 for none/D3(0)/D3(BJ)
 integer :: nmr_bas_choice = 0 ! pcSseg-1 by default
 ! 0/1/2 for pcSseg-1/pcSseg-2/def2TZVP for all atoms
 integer, allocatable :: nuc(:), def2_bas_typ(:)
 real(kind=8) :: pnmr_temperature = 298d0 ! 298 K
 real(kind=8), allocatable :: coor(:,:), sigma_orb(:), sigma_para(:)
 character(len=2), allocatable :: elem(:)

 ! method for the initial single-point, opt, and freq
 character(len=5), parameter :: method_def1 = 'B3LYP'
 character(len=8), parameter :: method_def2 = 'TPSSTPSS'
 character(len=30) :: method

 ! method for the NMR chemical shielding tensor
 character(len=4), parameter :: method_nmr_def = 'B972'
 character(len=30) :: method_nmr

 character(len=6) :: imp_sol_model = 'iefpcm' ! implicit solvent model
 character(len=30) :: solvent = 'chloroform' ! default solvent
 character(len=240) :: mokit_root = ' '
 character(len=240) :: gau_path = ' '
 character(len=240) :: orca_path = ' '

 logical :: uhf = .false. ! .True.: UHF/UKS; .False.: RHF/RKS
 logical :: skip_opt_freq = .false.
 ! whether to skip opt+freq and directly calculate the desired property

 logical, allocatable :: calc_pnmr(:) ! size natom
 ! calc_pnmr(i) = True/False: whether the sigma_para for the i-th atom has been
 !                            calculated
end module mol_and_calc_info

subroutine find_def2_bas_typ(natom, nuc, def2_bas_typ)
 implicit none
 integer :: i
 integer, intent(in) :: natom
 integer, intent(in) :: nuc(natom)
 integer, intent(out) :: def2_bas_typ(natom)
 ! 0: all-electron def2-SVP
 ! 1: def2-SVP with ECP
 ! 2: all-electron def2-TZVP
 ! 3: def2-TZVP with ECP

 def2_bas_typ = 0
 do i = 1, natom, 1
  select case(nuc(i))
  case(37,38,49:56,81:86)
   def2_bas_typ(i) = 1
  case(21:30)
   def2_bas_typ(i) = 2
  case(39:48,57:80)
   def2_bas_typ(i) = 3
  end select
 end do ! for i
end subroutine find_def2_bas_typ

! print gen/genecp def2SVP and def2TZVP basis sets into a specified .gjf
subroutine prt_gen_def2_bas2gjf(fid, natom, nuc, def2_bas_typ)
 use fch_content, only: nuc2elem
 use mol_and_calc_info, only: elem
 implicit none
 integer :: i, k
 integer, intent(in) :: fid, natom
 integer, intent(in) :: nuc(natom), def2_bas_typ(natom)
 logical :: prt, alive

 if(.not. allocated(elem)) then
  allocate(elem(natom))
  do i = 1, natom, 1
   elem(i) = nuc2elem(nuc(i))
  end do ! for i
 end if

 ! assuming that we just printed the Cartesian coordinates, so we need one
 ! blank line
 write(fid,'(/)',advance='no')

 if(ANY(def2_bas_typ < 2)) then
  k = 0
  do i = 1, natom, 1
   prt = .false.
   if(i == 1) then
    if(def2_bas_typ(1) < 2) prt = .true.
   else
    if(ALL(nuc(1:i-1)/=nuc(i)) .and. def2_bas_typ(i)<2) prt = .true.
   end if
   if(prt) then
    k = k + 1
    if(k == 1) then
     write(fid,'(A)',advance='no') TRIM(elem(i))
    else
     write(fid,'(A)',advance='no') ' '//TRIM(elem(i))
    end if
   end if
  end do ! for i
  write(fid,'(A,/,A,/,A)') ' 0', 'def2SVP', '****'
 end if

 if(ANY(def2_bas_typ > 1)) then
  k = 0
  do i = 1, natom, 1
   prt = .false.
   if(i == 1) then
    if(def2_bas_typ(1) > 1) prt = .true.
   else
    if(ALL(nuc(1:i-1)/=nuc(i)) .and. def2_bas_typ(i)>1) prt = .true.
   end if
   if(prt) then
    k = k + 1
    if(k == 1) then
     write(fid,'(A)',advance='no') TRIM(elem(i))
    else
     write(fid,'(A)',advance='no') ' '//TRIM(elem(i))
    end if
   end if
  end do ! for i
  write(fid,'(A,/,A,/,A)') ' 0', 'def2TZVP', '****'
 end if

 if(ANY(nuc > 36)) then
  write(fid,'(/)',advance='no')
  k = 0
  do i = 1, natom, 1
   prt = .false.
   alive = (def2_bas_typ(i)==1 .or. def2_bas_typ(i)==3)
   if(i == 1) then
    if(alive) prt = .true.
   else
    if(ALL(nuc(1:i-1)/=nuc(i)) .and. alive) prt = .true.
   end if
   if(prt) then
    k = k + 1
    if(k == 1) then
     write(fid,'(A)',advance='no') TRIM(elem(i))
    else
     write(fid,'(A)',advance='no') ' '//TRIM(elem(i))
    end if
   end if
  end do ! for i
  write(fid,'(/,A)') 'def2'
 end if

 write(fid,'(/)')
end subroutine prt_gen_def2_bas2gjf

! write/create a .gjf file for a single-point calculation
subroutine write_sp_gjf(gjfname)
 use fch_content, only: elem2nuc, nuc2elem
 use mol_and_calc_info, only: mem, nproc, method, charge, mult, natom, elem, &
  nuc, coor, def2_bas_typ, dis_type, imp_sol_model, solvent
 implicit none
 integer :: i, fid
 character(len=6) :: str6
 character(len=34), parameter :: error_warn='ERROR in subroutine write_sp_gjf: '
 character(len=240) :: chkname
 character(len=240), intent(in) :: gjfname

 if(.not. (allocated(nuc) .or. allocated(elem))) then
  write(6,'(/,A)') error_warn//'neither of `nuc` and `elem` array is allocated.'
  write(6,'(A)') 'gjfname='//TRIM(gjfname)
  stop
 else if(allocated(elem) .and. (.not.allocated(nuc))) then
  allocate(nuc(natom))
  do i = 1, natom, 1
   nuc(i) = elem2nuc(elem(i))
  end do ! for i
 else if(allocated(nuc) .and. (.not.allocated(elem))) then
  allocate(elem(natom))
  do i = 1, natom, 1
   elem(i) = nuc2elem(nuc(i))
  end do ! for i
 end if

 str6 = 'gen   '
 if(ANY(nuc > 86)) then
  write(6,'(/,A)') error_warn//'there exists some element which is beyond'
  write(6,'(A)') 'the range of def2- basis sets and def2-ECP.'
  write(6,'(A)') 'gjfname='//TRIM(gjfname)
  stop
 else if(ANY(nuc > 36)) then ! > Kr
  str6 = 'genecp'
 end if

 call find_specified_suffix(gjfname, '.gjf', i)
 chkname = gjfname(1:i-1)//'.chk'

 open(newunit=fid,file=TRIM(gjfname),status='replace')
 write(fid,'(A)') '%chk='//TRIM(chkname)
 write(fid,'(A,I0,A)') '%mem=', mem, 'GB'
 write(fid,'(A,I0)') '%nprocshared=', nproc
 write(fid,'(A)',advance='no') '#p '
 if(mult > 1) write(fid,'(A)',advance='no') 'U'
 write(fid,'(A)',advance='no') TRIM(method)//'/'//TRIM(str6)

 select case(dis_type)
 case(0) ! no dispersion correction
 case(1) ! D3zero, D3(0)
  write(fid,'(A)',advance='no') ' em=GD3'
 case(2) ! D3BJ, D3(BJ)
  write(fid,'(A)',advance='no') ' em=GD3BJ'
 case default
  write(fid,'(/,A,I0)') error_warn//'invalid dis_type=', dis_type
  close(fid)
  stop
 end select

 if(LEN_TRIM(solvent) > 0) then
  i = LEN_TRIM(imp_sol_model)
  if(i==0 .or. (i==6 .and. imp_sol_model(1:6)=='iefpcm')) then
   write(fid,'(A)',advance='no') ' scrf(solvent='//TRIM(solvent)//')'
  else
   write(fid,'(A)',advance='no') ' scrf('//TRIM(imp_sol_model)//',solvent='//&
                                 TRIM(solvent)//')'
  end if
 end if
 write(fid,'(A)') ' nosymm int=nobasistransform scf(xqc,maxcycle=500) stable=opt'
 write(fid,'(/,A,/)') 'generated by MOKIT workflow'
 write(fid,'(I0,1X,I0)') charge, mult

 do i = 1, natom, 1
  write(fid,'(A2,3(1X,F18.8))') elem(i), coor(:,i)
 end do ! for i

 allocate(def2_bas_typ(natom))
 call find_def2_bas_typ(natom, nuc, def2_bas_typ)
 ! TODO: also check the connectivity of each atom. For hypervalent atoms, change
 ! to def2-TZVP automatically.
 call prt_gen_def2_bas2gjf(fid, natom, nuc, def2_bas_typ)

 close(fid)
end subroutine write_sp_gjf

! write/create a .gjf file for the single-point calculation using method_nmr
subroutine write_nmr_gjf1(gjfname, old_chk)
 use mol_and_calc_info, only: mem, nproc, uhf, method_nmr, imp_sol_model, &
  solvent, mokit_root, natom, elem, nuc, nmr_bas_choice, def2_bas_typ
 implicit none
 integer :: i, k, fid
 character(len=6) :: str6
 character(len=36), parameter :: error_warn = 'ERROR in subroutine write_nmr_gj&
                                              &f1: '
 character(len=240) :: chkname
 character(len=240), intent(in) :: gjfname, old_chk
 logical :: prt, alive

 if(.not. allocated(def2_bas_typ)) then
  write(6,'(/,A)') error_warn//'the integer array def2_bas_typ'
  write(6,'(A)') 'is not allocated. Please initialize it before using this subr&
                 &outine.'
  write(6,'(A)') 'gjfname='//TRIM(gjfname)
  stop
 end if
 str6 = 'gen   '
 if(ANY(def2_bas_typ==1) .or. ANY(def2_bas_typ==3)) str6 = 'genecp'

 call find_specified_suffix(gjfname, '.gjf', i)
 chkname = gjfname(1:i-1)//'.chk'
 if(TRIM(old_chk) /= TRIM(chkname)) then
  call sys_copy_file(TRIM(old_chk), TRIM(chkname), .false.)
 end if

 call find_specified_suffix(gjfname, '.gjf', i)
 chkname = gjfname(1:i-1)//'.chk'

 open(newunit=fid,file=TRIM(gjfname),status='replace')
 write(fid,'(A)') '%chk='//TRIM(chkname)
 write(fid,'(A,I0,A)') '%mem=', mem, 'GB'
 write(fid,'(A,I0)') '%nprocshared=', nproc
 write(fid,'(A)',advance='no') '#p '
 if(uhf) write(fid,'(A)',advance='no') 'U'
 write(fid,'(A)',advance='no') TRIM(method_nmr)//'/'//TRIM(str6)
 ! DFT-D3 does not affect the NMR chemical shielding tensor, so here we will not
 ! print em=GD3/GD3BJ

 if(LEN_TRIM(solvent) > 0) then
  i = LEN_TRIM(imp_sol_model)
  if(i==0 .or. (i==6 .and. imp_sol_model(1:6)=='iefpcm')) then
   write(fid,'(A)',advance='no') ' scrf(solvent='//TRIM(solvent)//')'
  else
   write(fid,'(A)',advance='no') ' scrf('//TRIM(imp_sol_model)//',solvent='//&
                                 TRIM(solvent)//')'
  end if
 end if

 write(fid,'(A,/)') ' nosymm int=nobasistransform scf(xqc,maxcycle=500,NoVarAcc&
                     &) guess=read geom=allcheck stable=opt'

 select case(nmr_bas_choice)
 case(0) ! pcSseg-1 for all atoms
  write(fid,'(A)') '@'//TRIM(mokit_root)//'/mokit/basis/pcSseg-1/N'
 case(1) ! pcSseg-2 for all atoms
  write(fid,'(A)') '@'//TRIM(mokit_root)//'/mokit/basis/pcSseg-2/N'
 case(2) ! def2-TZVP for all atoms
  k = 0
  do i = 1, natom, 1
   prt = .false.
   if(i == 1) then
    prt = .true.
   else ! i > 1
    if(ALL(nuc(1:i-1) /= nuc(i))) prt = .true.
   end if
   if(prt) then
    k = k + 1
    if(k == 1) then
     write(fid,'(A)',advance='no') TRIM(elem(i))
    else
     write(fid,'(A)',advance='no') ' '//TRIM(elem(i))
    end if
   end if
  end do ! for i
  write(fid,'(A,/,A,/,A)') ' 0', 'def2TZVP', '****'

  if(ANY(nuc > 36)) then
   write(fid,'(/)',advance='no')
   k = 0
   do i = 1, natom, 1
    prt = .false.
    alive = (def2_bas_typ(i)==1 .or. def2_bas_typ(i)==3)
    if(i == 1) then
     if(alive) prt = .true.
    else
     if(ALL(nuc(1:i-1)/=nuc(i)) .and. alive) prt = .true.
    end if
    if(prt) then
     k = k + 1
     if(k == 1) then
      write(fid,'(A)',advance='no') TRIM(elem(i))
     else
      write(fid,'(A)',advance='no') ' '//TRIM(elem(i))
     end if
    end if
   end do ! for i
   write(fid,'(/,A)') 'def2'
  end if
 case default
  write(6,'(/,A,I0)') error_warn//'invalid nmr_bas_choice=', nmr_bas_choice
  write(6,'(A)') 'Currently only 0/1/2 is allowed.'
  write(6,'(A)') 'gjfname='//TRIM(gjfname)
  stop
 end select

 write(fid,'(/)')
 close(fid)
end subroutine write_nmr_gjf1

! write/create a .gjf file for NMR shielding tensor calculation using method_nmr
subroutine write_nmr_gjf2(gjfname, old_chk)
 use mol_and_calc_info, only: mem, nproc, uhf, method_nmr, imp_sol_model, solvent
 implicit none
 integer :: i, fid
 character(len=36), parameter :: error_warn = 'ERROR in subroutine write_nmr_gj&
                                              &f2: '
 character(len=240) :: chkname
 character(len=240), intent(in) :: gjfname, old_chk

 call find_specified_suffix(gjfname, '.gjf', i)
 chkname = gjfname(1:i-1)//'.chk'
 if(TRIM(old_chk) /= TRIM(chkname)) then
  call sys_copy_file(TRIM(old_chk), TRIM(chkname), .false.)
 end if

 open(newunit=fid,file=TRIM(gjfname),status='replace')
 write(fid,'(A)') '%chk='//TRIM(chkname)
 write(fid,'(A,I0,A)') '%mem=', mem, 'GB'
 write(fid,'(A,I0)') '%nprocshared=', nproc
 write(fid,'(A)',advance='no') '#p '
 if(uhf) write(fid,'(A)',advance='no') 'U'
 write(fid,'(A)',advance='no') TRIM(method_nmr)//' chkbasis'
 ! The NMR basis set was defined in subroutine write_nmr_gjf1. Here we simply
 ! read the basis set using `chkbasis`.
 ! DFT-D3 does not affect the NMR chemical shielding tensor, so here we will not
 ! print em=GD3/GD3BJ

 if(LEN_TRIM(solvent) > 0) then
  i = LEN_TRIM(imp_sol_model)
  if(i==0 .or. (i==6 .and. imp_sol_model(1:6)=='iefpcm')) then
   write(fid,'(A)',advance='no') ' scrf(solvent='//TRIM(solvent)//')'
  else
   write(fid,'(A)',advance='no') ' scrf('//TRIM(imp_sol_model)//',solvent='//&
                                 TRIM(solvent)//')'
  end if
 end if

 write(fid,'(A,/)') ' nosymm int=nobasistransform scf=NoVarAcc guess=read geom=&
                    &allcheck NMR'
 close(fid)
end subroutine write_nmr_gjf2

! write/create a .gjf file for routine calculations
subroutine write_routine_gjf(gjfname, old_chk, calc_type, hessian_type)
 use mol_and_calc_info, only: mem, nproc, method, uhf, dis_type, imp_sol_model,&
  solvent
 implicit none
 integer :: i, fid
 integer, intent(in) :: calc_type, hessian_type
 ! calc_type: stable(rrhf,opt)/stable=opt/opt/freq/opt+freq
 ! hessian_type: none/calcfc/rcfc/recalc=5
 character(len=39), parameter :: error_warn = 'ERROR in subroutine write_routin&
                                              &e_gjf: '
 character(len=240) :: chkname
 character(len=240), intent(in) :: gjfname, old_chk

 call find_specified_suffix(gjfname, '.gjf', i)
 chkname = gjfname(1:i-1)//'.chk'
 if(TRIM(old_chk) /= TRIM(chkname)) then
  call sys_copy_file(TRIM(old_chk), TRIM(chkname), .false.)
 end if

 open(newunit=fid,file=TRIM(gjfname),status='replace')
 write(fid,'(A)') '%chk='//TRIM(chkname)
 write(fid,'(A,I0,A)') '%mem=', mem, 'GB'
 write(fid,'(A,I0)') '%nprocshared=', nproc

 select case(calc_type)
 case(1,2) ! stable(rrhf,opt), stable=opt
  write(fid,'(A)',advance='no') '#p '
 case(3,5) ! opt/opt+freq
  select case(hessian_type)
  case(1)
   write(fid,'(A)',advance='no') '#p opt(maxcycles=250) '
  case(2)
   write(fid,'(A)',advance='no') '#p opt(calcfc,maxcycles=250) '
  case(3)
   write(fid,'(A)',advance='no') '#p opt(rcfc,maxcycles=250) '
  case(4)
   write(fid,'(A)',advance='no') '#p opt(recalc=5,maxcycles=250) '
  case default
   write(6,'(/,A,I0)') error_warn//'invalid hessian_type=', hessian_type
   write(6,'(A)') 'gjfname='//TRIM(gjfname)
   stop
  end select
  if(calc_type == 5) write(fid,'(A)',advance='no') 'freq '
 case(4) ! freq
  write(fid,'(A)',advance='no') '#p freq '
 case default
  write(6,'(/,A,I0)') error_warn//'invalid calc_type=', calc_type
  write(6,'(A)') 'gjfname='//TRIM(gjfname)
  stop
 end select

 if(uhf) write(fid,'(A)',advance='no') 'U'
 write(fid,'(A)',advance='no') TRIM(method)//' chkbasis'

 select case(dis_type)
 case(0) ! no dispersion correction
 case(1) ! D3zero, D3(0)
  write(fid,'(A)',advance='no') ' em=GD3'
 case(2) ! D3BJ, D3(BJ)
  write(fid,'(A)',advance='no') ' em=GD3BJ'
 case default
  write(fid,'(/,A,I0)') error_warn//'invalid dis_type=', dis_type
  close(fid)
  stop
 end select

 if(LEN_TRIM(solvent) > 0) then
  i = LEN_TRIM(imp_sol_model)
  if(i==0 .or. (i==6 .and. imp_sol_model(1:6)=='iefpcm')) then
   write(fid,'(A)',advance='no') ' scrf(solvent='//TRIM(solvent)//')'
  else
   write(fid,'(A)',advance='no') ' scrf('//TRIM(imp_sol_model)//',solvent='//&
                                 TRIM(solvent)//')'
  end if
 end if
 write(fid,'(A)',advance='no') ' nosymm int=nobasistransform scf(xqc,maxcycle=5&
                               &00,NoVarAcc) guess=read geom=allchec'

 select case(calc_type)
 case(1) ! stable(rrhf,opt)
  write(fid,'(A,/)') 'k stable(rrhf,opt)'
 case(2) ! stable=opt
  write(fid,'(A,/)') 'k stable=opt'
 case(3,4,5) ! opt/freq/opt+freq
  write(fid,'(A,/)') 'k'
 end select

 close(fid)
end subroutine write_routine_gjf

! increase the integer after the `s` character in filename
! e.g. CuO_s1.gjf -> CuO_s2.gjf
subroutine increase_int_after_s_in_filename(fname)
 implicit none
 integer :: i, j, k, m, n
 character(len=54), parameter :: error_warn = 'ERROR in subroutine increase_int&
                                              &_after_s_in_filename: '
 character(len=240) :: buf
 character(len=240), intent(inout) :: fname

 m = LEN_TRIM(fname)
 j = INDEX(fname(1:m), 's', back=.true.)
 k = INDEX(fname(1:m), '.', back=.true.)
 if(j==0 .or. k==0 .or. (k-1<j+1)) then
  write(6,'(/,A)') error_warn//'illegal filename for this subroutine.'
  write(6,'(A)') 'fname='//TRIM(fname)
  stop
 end if

 i = 0
 read(fname(j+1:k-1),fmt=*,iostat=n) i
 if(n /= 0) then
  write(6,'(/,A)') error_warn//'no integer is found in filename'
  write(6,'(A)') TRIM(fname)
  stop
 end if

 write(buf,'(A,I0,A)') fname(1:j),i+1,fname(k:m)
 fname = buf
end subroutine increase_int_after_s_in_filename

! perform `opt` and `stable=opt` using a loop
subroutine do_routine_opt_and_stable_opt(proname, max_nstep, istep)
 use mol_and_calc_info, only: gau_path
 use util_wrapper, only: formchk
 implicit none
 integer :: k
 integer, intent(in) :: max_nstep
!f2py intent(in) :: max_nstep
 integer, intent(inout) :: istep
!f2py intent(inout) :: istep
 real(kind=8) :: scf_e(2), ssquare(2)
 real(kind=8), parameter :: e_diff = 1d-4
 character(len=240) :: inpname, old_chk, chkname, outname
 character(len=240), intent(in) :: proname
!f2py intent(in) :: proname

 k = LEN_TRIM(proname)
 inpname = proname(1:k)//'.gjf'
 old_chk = proname(1:k)//'.chk'
 chkname = old_chk
#ifdef _WIN32
  outname = proname(1:k)//'.out'
#else
  outname = proname(1:k)//'.log'
#endif

 k = 1
 do while(k <= max_nstep)
  ! read converged MOs and perform the geometry optimization for the input
  ! geometry
  istep = istep + 1
  write(6,'(/,A,I0,A)') 'Step ',istep,':'
  call increase_int_after_s_in_filename(inpname)
  call increase_int_after_s_in_filename(chkname)
  call increase_int_after_s_in_filename(outname)
  call write_routine_gjf(inpname, old_chk, 3, 1)
  call submit_gau_job(gau_path, inpname, .true.)
  call formchk(chkname)
  call read_hf_e_and_ss_from_gau_log(outname, .false., scf_e(1), ssquare(1))
  write(6,'(A,F20.8,A,F10.3)') 'The electronic energy of the optimized geometry &
                               &is: ', scf_e(1), ', <S^2>=', ssquare(1)

  ! Do `stable=opt`, i.e. check the wave function stability of the optimized
  ! geometry. If not stable, `stable=opt` will get stable wave function.
  istep = istep + 1
  write(6,'(/,A,I0,A)') 'Step ',istep,':'
  call increase_int_after_s_in_filename(inpname)
  call increase_int_after_s_in_filename(old_chk)
  call increase_int_after_s_in_filename(chkname)
  call increase_int_after_s_in_filename(outname)
  call write_routine_gjf(inpname, old_chk, 2, 1)
  call submit_gau_job(gau_path, inpname, .true.)
  call formchk(chkname)

  ! Now we check the SCF energies in the output file. If there exists one lower
  ! energy, we need to read this new energy as well as new MOs, and continue to
  ! optimize the geometry.
  call read_hf_e_and_ss_from_gau_log(outname, .true., scf_e(1), ssquare(1))
  call read_hf_e_and_ss_from_gau_log(outname, .false., scf_e(2), ssquare(2))
  if(scf_e(1)-scf_e(2) > e_diff) then
   write(6,'(A)') 'Remark: `stable=opt` finds a lower energy for the current opt&
                  &imized geometry.'
   write(6,'(A)') 'Therefore, we need to continue the geometry optimization.'
  else
   write(6,'(A)') 'The wave function of the optimized geometry is stable.'
   exit
  end if

  k = k + 1
 end do ! for while

 if(k > max_nstep) then
  write(6,'(/,A)') 'ERROR in subroutine do_routine_opt_and_stable_opt: the time&
                   &s of `opt` and'
  write(6,'(A,I0)') '`stable=opt` exceed max_nstep=',max_nstep
  write(6,'(A)') 'proname='//TRIM(proname)
  stop
 end if
end subroutine do_routine_opt_and_stable_opt

! check whether there is any imaginary frequency in a specified Gaussian output file
subroutine check_imag_freq_in_log(logname, has_imag_freq)
 implicit none
 integer :: i, fid
 real(kind=8) :: r
 character(len=240) :: buf
 character(len=240), intent(in) :: logname
 logical, intent(out) :: has_imag_freq

 has_imag_freq = .false.; r = 0d0
 open(newunit=fid,file=TRIM(logname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:12) == 'Frequencies') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine check_imag_freq_in_log: failed to locat&
                   &e "Frequencies"'
  write(6,'(A)') 'in file '//TRIM(logname)
  stop
 end if

 i = INDEX(buf, '--')
 read(buf(i+2:),*) r
 if(r < 0d0) has_imag_freq = .true.
end subroutine check_imag_freq_in_log

! modify a given ORCA .inp file to make it suitable for pNMR calculation
! Note: the input `solvent` must be in ORCA convention
subroutine modify_orca_pnmr_inp(inpname, mult, natom, nuc, solvent)
 use mol_and_calc_info, only: mem, nproc, calc_pnmr
 use fch_content, only: nuc2elem
 implicit none
 integer :: i, nproc1, mem1, fid, fid1, RENAME
 integer, intent(in) :: mult, natom
 integer, intent(in) :: nuc(natom)
 character(len=30), intent(in) :: solvent
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname
 logical :: has_h, has_c

 has_h = .false.; has_c = .false.
 if(ANY(nuc == 1)) has_h = .true.
 if(ANY(nuc == 6)) has_c = .true.

 call find_specified_suffix(inpname, '.inp', i)
 inpname1 = inpname(1:i-1)//'.t'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 ! Here we use only half of CPU cores for ORCA, since it utilizes the MPI
 ! parallelism and requires larger memory when the number of CPU cores increases.
 call reduce_nproc_and_enlarge_mem(nproc, mem, 2, nproc1, mem1)
 write(fid1,'(A,I0,A)') '%pal nprocs ', nproc1, ' end'
 write(fid1,'(A,I0)') '%maxcore ', mem1

 do i = 1, 3
  read(fid,'(A)') buf ! skip 3 lines
 end do ! for i

 if(buf(1:1) /= '!') then
  close(fid)
  close(fid1,status='delete')
  write(6,'(/,A)') 'ERROR in subroutine modify_orca_pnmr_inp: unexpected string.'
  write(6,'(A)') 'buf='//TRIM(buf)
  write(6,'(A)') 'inpname='//TRIM(inpname)
  stop
 end if

 if(LEN_TRIM(solvent) > 0) then
  write(fid1,'(A)') TRIM(buf)//' CPCM('//TRIM(solvent)//')'
 else
  write(fid1,'(A)') TRIM(buf)
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 write(fid1,'(A)') '%eprnmr'
 write(fid1,'(A)') ' gtensor 1'
 if(mult > 2) write(fid1,'(A)') ' dtensor ssandso'
 if(has_h) write(fid1,'(A)') ' Nuclei = all H { aiso,adip,aorb }'
 if(has_c) write(fid1,'(A)') ' Nuclei = all C { aiso,adip,aorb }'

 if(has_h .or. has_c) then ! calculate sigma_para for H and C atoms
  allocate(calc_pnmr(natom), source=.false.)
  where(nuc==1 .or. nuc==6) calc_pnmr = .true.
 else                      ! calculate sigma_para for all atoms
  write(fid1,'(A)') ' Nuclei = all '//TRIM(nuc2elem(nuc(1)))//' { aiso,adip,aor&
                    &b }'
  do i = 2, natom, 1
   if(ANY(nuc(1:i-1) == nuc(i))) cycle
   write(fid1,'(A)') ' Nuclei = all '//TRIM(nuc2elem(nuc(i)))//' { aiso,adip,ao&
                     &rb }'
  end do ! for i
  allocate(calc_pnmr(natom), source=.true.)
 end if
 write(fid1,'(A)') 'end'

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine modify_orca_pnmr_inp

! read NMR workflow parameters from the Title Card line of a specified .gjf
subroutine read_nmr_wf_para_from_gjf(gjfname, prt)
 use mol_and_calc_info, only: mult, nmr_bas_choice, method_def1, method_def2, &
  method_nmr_def, method, method_nmr
 implicit none
 integer :: i, j, fid
 character(len=47), parameter :: error_warn = 'ERROR in subroutine read_nmr_wf_&
                                              &para_from_gjf: '
 character(len=240) :: buf
 character(len=240), intent(in) :: gjfname
 logical, intent(in) :: prt

 nmr_bas_choice = 0
 if(mult == 1) then
  method = method_def1
 else
  method = method_def2
 end if
 method_nmr = method_nmr_def

 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'wrong syntax in file'
  write(6,'(A)') 'gjfname='//TRIM(gjfname)
  close(fid)
  stop
 end if

 ! the Title Card line is supposed to be after the 1st blank line
 read(fid,'(A)') buf
 close(fid)

 i = INDEX(buf, '{'); j = INDEX(buf, '}', back=.true.)
 if(i==0 .and. j==0) then
  if(prt) then
   write(6,'(/,A,I0)') 'method='//TRIM(method)//', method_nmr='//&
             TRIM(method_nmr)//', nmr_bas_choice=', nmr_bas_choice
  end if
  return
 end if

 if((i==0 .and. j>0) .or. (i>j)) then
  write(6,'(/,A)') error_warn//'wrong {} in file'
  write(6,'(A)') TRIM(gjfname)
  stop
 end if

 if(buf(i+1:i+15) == 'nmr_bas_choice=') then
  read(buf(i+16:j-1),*) nmr_bas_choice
  select case(nmr_bas_choice)
  case(0,1,2)
  case default
   write(6,'(/,A,I0)') error_warn//'invalid nmr_bas_choice=', nmr_bas_choice
   write(6,'(A)') 'gjfname='//TRIM(gjfname)
   stop
  end select
 else
  write(6,'(/,A)') error_warn//'currently only `nmr_bas_choice`'
  write(6,'(A)') 'is supported in {}.'
  write(6,'(A)') 'gjfname='//TRIM(gjfname)
  stop
 end if

 if(prt) then
  write(6,'(/,A,I0)') 'method='//TRIM(method)//', method_nmr='//&
            TRIM(method_nmr)//', nmr_bas_choice=', nmr_bas_choice
 end if
end subroutine read_nmr_wf_para_from_gjf

! This is an NMR calculation workflow, which includes three steps: opt, freq and
!  NMR. The conformational search is not considered in this subroutine.
! The density functional in opt and freq will be read from `gjfname`.
! The basis set used in opt and freq will be fixed as def2-SVP for non-metal
!  atoms, and def2-TZVP for metal atoms.
! The theory level used in NMR will be fixed as B972/pcSseg-2.
! For singlet molecule, the wave function stability will be checked after the
!  RKS geometry optimization. If instability is found, broken-symmetry UKS will
!  be invoked. For non-singlet molecule, the wave function stability will be
!  checked both for the initial geometry and optimized geometry.
! The wave function stability will also be checked in the B972/pcSseg-2 single-
!  point calculation, followed by the NMR calculation at the same theory level.
subroutine nmr_workflow(gjfname)
 use mol_and_calc_info, only: mem, nproc, natom, charge, mult, nuc, elem, coor,&
  dis_type, def2_bas_typ, uhf, imp_sol_model, solvent, mokit_root, gau_path, &
  orca_path, method_nmr, calc_pnmr, sigma_orb, sigma_para, pnmr_temperature
 use util_wrapper, only: formchk, fch2mkl_wrap, mkl2gbw
 implicit none
 integer :: i, k, istep
 integer, parameter :: max_nstep = 12 ! 12 steps
 integer, parameter :: max_nstep2 = 4 ! 4 steps
 integer(kind=4) :: hostnm
 real(kind=8) :: scf_e, ssquare, r(3)
 character(len=8) :: hostname
 character(len=24) :: data_string
 character(len=34), parameter :: error_warn='ERROR in subroutine nmr_workflow: '
 character(len=240), external :: get_mokit_root 
 character(len=240) :: proname, fchname, mklname, pnmr_out, qro_name, uno_name,&
  unso_name, ao_xyz1, ao_xyz2
 character(len=240), allocatable :: inpname(:), chkname(:), outname(:)
 character(len=240), intent(in) :: gjfname
!f2py intent(in) :: gjfname
 character(len=500) :: longbuf
 logical :: has_imag_freq

 hostname = ' '
 data_string = ' '
 i = hostnm(hostname)
 call fdate(data_string)
 write(6,'(A)') 'HOST '//TRIM(hostname)//', '//TRIM(data_string)
 write(6,'(A)') 'NMR/pNMR calculation workflow of MOKIT.'

 write(6,'(/,A)') 'Read program paths from environment variables:'
 mokit_root = get_mokit_root()
 call get_gau_path(gau_path)
 call get_exe_path('orca', orca_path)
 write(6,'(A)') 'MOKIT_ROOT= '//TRIM(mokit_root)
 write(6,'(A)') 'gau_path  = '//TRIM(gau_path)
 write(6,'(A)') 'orca_path = '//TRIM(orca_path)
 call check_exe_exist(gau_path)
 call check_exe_exist(orca_path)
 allocate(inpname(max_nstep), chkname(max_nstep), outname(max_nstep))

 call find_specified_suffix(gjfname, '.gjf', k)
 do i = 1, max_nstep, 1
  write(inpname(i),'(A,I0,A)') gjfname(1:k-1)//'_s',i,'.gjf'
  write(chkname(i),'(A,I0,A)') gjfname(1:k-1)//'_s',i,'.chk'
#ifdef _WIN32
  write(outname(i),'(A,I0,A)') gjfname(1:k-1)//'_s',i,'.out'
#else
  write(outname(i),'(A,I0,A)') gjfname(1:k-1)//'_s',i,'.log'
#endif
 end do ! for i
 ! here `_s1` means step 1

 call read_mem_and_nproc_from_gjf(gjfname, mem, nproc)
 ! Note: the return `mem` is in unit MB, here we convert it to unit GB
 mem = NINT(DBLE(mem)/1d3)
 write(6,'(/,2(A,I0))') 'mem=', mem, 'GB, nproc=', nproc

 call read_natom_from_gjf(gjfname, natom)
 allocate(elem(natom), nuc(natom), coor(3,natom))
 call read_elem_and_coor_from_gjf(gjfname, natom, elem, nuc, coor, charge, mult)
 write(6,'(2(A,I0))') 'charge=', charge, ', mult=', mult

 call read_nmr_wf_para_from_gjf(gjfname, .true.)
 write(6,'(A)',advance='no') 'Dispersion correction: '
 select case(dis_type)
 case(0) ! no dispersion correction
  write(6,'(A)') 'none'
 case(1) ! D3(0)
  write(6,'(A)') 'em=GD3, i.e. D3(0)'
 case(2) ! D3BJ
  write(6,'(A)') 'em=GD3BJ, i.e. D3(BJ)'
 case default
  write(6,'(/,A,I0)') error_warn//'invalid dis_type=', dis_type
  stop
 end select

 call read_scrf_model_and_solvent_from_gjf(gjfname, imp_sol_model, solvent)
 if(LEN_TRIM(imp_sol_model) == 0) then
  write(6,'(A)') 'Implicit solvent model: none'
 else
  write(6,'(A)') 'Implicit solvent model: '//TRIM(imp_sol_model)
 end if
 if(LEN_TRIM(solvent) == 0) then
  write(6,'(A)') 'Implicit solvent: none'
 else
  write(6,'(A)') 'Implicit solvent: '//TRIM(solvent)
 end if

 ! perform a single-point calculation for the input geometry
 write(6,'(/,A)') 'Step 1:'
 write(6,'(A)') 'Single-point calculation for the initial geometry.'
 call write_sp_gjf(inpname(1))
 call submit_gau_job(gau_path, inpname(1), .true.)
 call find_specified_suffix(chkname(1), '.chk', i)
 fchname = chkname(1)(1:i-1)//'.fch'
 call formchk(chkname(1), fchname)
 call check_uhf_in_fch(fchname, uhf)
 call read_hf_e_and_ss_from_gau_log(outname(1), .false., scf_e, ssquare)
 write(6,'(A,F20.8,A,F10.3)') 'The electronic energy of the initial geometry is&
                              &: ', scf_e, ', <S^2>=', ssquare
 istep = 1

 ! perform `opt` and `stable=opt` in using a loop
 call find_specified_suffix(inpname(1), '.gjf', i)
 proname = inpname(1)(1:i-1)
 call do_routine_opt_and_stable_opt(proname, max_nstep2, istep)
 ! Note: do not use inpname(1)(1:i-1) directly since its length is less than
 ! 240. Use proname.

 ! perform the harmonic frequency calculation, i.e. freq
 istep = istep + 1
 write(6,'(/,A,I0,A)') 'Step ',istep,':'
 call write_routine_gjf(inpname(istep), chkname(istep-1), 4, 1)
 call submit_gau_job(gau_path, inpname(istep), .true.)
 call formchk(chkname(istep))
 ! Check whether there is any imaginary frequency. If any, read Hessian and
 ! continue the geometry optimization.
 call check_imag_freq_in_log(outname(istep), has_imag_freq)
 if(has_imag_freq) then
  write(6,'(A)') 'Remark: the harmonic frequency calculation is done. But there&
                 & exists one or'
  write(6,'(A)') 'more imaginary frequency. This workflow will read Hessian fro&
                 &m .chk file and'
  write(6,'(A)') 'continue the geometry optimization.'
  stop
 else
  write(6,'(A)') 'The harmonic frequency calculation is accomplished and there &
                 &is no imaginary'
  write(6,'(A)') 'frequency.'
 end if

 istep = istep + 1
 write(6,'(/,A,I0,A)') 'Step ',istep,':'
 write(6,'(A)') 'Use the NMR theory level to perform a single-point calculation&
                & and check wave'
 write(6,'(A)') 'function stability.'
 call write_nmr_gjf1(inpname(istep), chkname(istep-1))
 call submit_gau_job(gau_path, inpname(istep), .true.)
 call formchk(chkname(istep))

 istep = istep + 1
 write(6,'(/,A,I0,A)') 'Step ',istep,':'
 write(6,'(A)') 'Calculate the NMR shielding tensor sigma_orb, i.e. orbital con&
                &tribution.'
 call write_nmr_gjf2(inpname(istep), chkname(istep-1))
 call submit_gau_job(gau_path, inpname(istep), .true.)
 call formchk(chkname(istep))
 allocate(sigma_orb(natom), source=0d0)
 call read_nmr_iso_shield_from_gau_log(outname(istep), natom, sigma_orb)

 if(mult == 1) then ! singlet
  write(6,'(/,A)') 'NMR isotropic shielding of all atoms (ppm):'
  write(6,'(A)') 'Label Element   sigma_orb'
  do i = 1, natom, 1
   write(6,'(I5,6X,A2,1X,F11.3)') i, elem(i), sigma_orb(i)
  end do ! for i
 else               ! For non-singlet, we need to calculate sigma_para
  istep = istep + 1
  write(6,'(/,A,I0,A)') 'Step ',istep,':'
  if(LEN_TRIM(imp_sol_model) > 0) then
   write(6,'(A)') 'Remark: `scrf` is specified in Gaussian. No matter which imp&
                  &licit solvent'
   write(6,'(A)') 'model used in Gaussian, here we will use CPCM in ORCA.'
  end if
  call find_specified_suffix(chkname(istep-1), '.chk', i)
  fchname = chkname(istep-1)(1:i-1)//'.fch'
  mklname = chkname(istep-1)(1:i-1)//'.mkl'
  call increase_int_after_s_in_filename(mklname)
  call fch2mkl_wrap(fchname, mklname, method_nmr, .false.)
  call mkl2gbw(mklname)
  call change_suffix_in_fname(inpname(istep), '.gjf', '.inp')
#ifndef _WIN32
  call change_suffix_in_fname(outname(istep), '.log', '.out')
#endif
  call modify_orca_pnmr_inp(inpname(istep), mult, natom, nuc, solvent)
  call submit_orca_job(orca_path, inpname(istep), .true., .false., .false.)
  call find_specified_suffix(inpname(istep), '.inp', i)
  proname = inpname(istep)(1:i-1)
  pnmr_out = inpname(istep)(1:i-1)//'_pnmr.out'
  qro_name = inpname(istep)(1:i-1)//'.qro'
  uno_name = inpname(istep)(1:i-1)//'.uno'
  unso_name = inpname(istep)(1:i-1)//'.unso'
  write(ao_xyz1,'(A,I0,A)') inpname(istep)(1:i-1)//'.Dori_SCF.mult', mult, &
                            '.dsoc0.re.el..AO.xyz'
  write(ao_xyz2,'(A,I0,A)') inpname(istep)(1:i-1)//'.gori_SCF.mult', mult, &
                            '.deribgiaoX.re.spin..AO.xyz'
  longbuf = 'orca_pnmr '//TRIM(proname)//' >'//TRIM(pnmr_out)//' 2>&1'
  call run_command(TRIM(longbuf), .false., .false.)
  allocate(sigma_para(natom))
  call read_para_shield_from_orca_pnmr_out2(pnmr_out, natom, sigma_para)
  write(6,'(/,A,F7.2,A)') 'Temperature used in calculating sigma_para: ', &
                           pnmr_temperature, ' K'
  write(6,'(A)') 'pNMR isotropic orbital shielding and paramagnetic shielding o&
                 &f target atoms (ppm):'
  write(6,'(A)') 'Label Element   sigma_orb  sigma_para   sigma_tot'
  do i = 1, natom, 1
   if(calc_pnmr(i)) then
    r(1) = sigma_orb(i); r(2) = sigma_para(i); r(3) = r(1)+r(2)
    write(6,'(I5,6X,A2,3(1X,F11.3))') i, elem(i), r
   end if
  end do ! for i

  call delete_files(5, [qro_name, uno_name, unso_name, ao_xyz1, ao_xyz2])
  deallocate(sigma_para, calc_pnmr)
 end if

 ! remove all .chk files since they are useless
 do i = 1, max_nstep, 1
  call delete_file(TRIM(chkname(i)))
 end do ! for i
 deallocate(elem, nuc, coor, def2_bas_typ, inpname, chkname, outname, sigma_orb)

 call fdate(data_string)
 write(6,'(/,A)') 'Please remember to cite MOKIT in your future publication reg&
                  &arding the NMR/pNMR'
 write(6,'(A)') 'calculation workflow.'
 write(6,'(A)') 'Normal termination of NMR workflow at '//TRIM(data_string)
end subroutine nmr_workflow

