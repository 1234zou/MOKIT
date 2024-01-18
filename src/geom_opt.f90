! written by jxzou at 2031201: perform geometry optimization using Gaussian as
!  the optimizer

module mol_comp_info
 implicit none
 integer :: charge = 0
 integer :: mult = 1
 integer :: natom = 0
 integer :: mem = 4000 ! unit in MB
 integer :: nproc = 4
 integer :: job_type = 1 ! 1/2/3/4/5 for opt/freq/opt+freq/IRC/rigid_scan
 ! Single point calculation is not allowed in this program.
 integer :: wfn_type = 0
 real(kind=8), allocatable :: coor(:,:)   ! coor(3,natom)
 character(len=2), allocatable :: elem(:) ! elem(natom)
 character(len=11) :: method = 'wB97M-V'
 character(len=21) :: basis = 'def2TZVP'
 character(len=240) :: gau_path
 logical :: numfreq = .false. ! numerical frequency
 logical :: init_scf = .false.
 ! True: perform an SCF computation using Gaussian for the initial geometry, and
 !       transfer converged MOs to ORCA
 ! False: generate only SCF initial guess using Gaussian for the initial geometry
 !       and transfer initial guess MOs to ORCA
contains

subroutine get_job_type_from_gjf(gjfname)
 implicit none
 integer :: i, j, fid
 character(len=1) :: str
 character(len=240) :: buf
 character(len=240), intent(in) :: gjfname

 str = ' '
 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:1) == '#') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine get_job_type_from_gjf: problematic file.'
  write(6,'(A)') 'Filename='//TRIM(gjfname)
  close(fid)
  stop
 end if

 read(fid,'(A)') str
 close(fid)
 if(str /= ' ') then
  write(6,'(/,A)') 'ERROR in subroutine get_job_type_from_gjf: please write all&
                   & keywords in one'
  write(6,'(A)') 'line. Multiple lines of keywords are not supported in geom_opt.'
  write(6,'(A)') 'str='//str
  stop
 end if

 call lower(buf)
 call read_method_and_basis_from_buf(buf, method, basis, wfn_type)

 i = INDEX(buf, 'opt'); j = INDEX(buf, 'freq')
 if(i>0 .and. j>0) then
  job_type = 3
 else if(j > 0) then
  job_type = 2
 end if

 i = INDEX(buf, 'irc')
 if(i > 0) job_type = 4
 i = INDEX(buf, 'scan')
 if(i > 0) job_type = 5

 i = INDEX(buf, 'freq=numer'); j = INDEX(buf, 'freq(numer')
 if(i>0 .or. j>0) numfreq = .true.

 if(TRIM(method)=='wb97m-v' .and. (job_type==2 .or. job_type==3)) then
  numfreq = .true.
 end if
end subroutine get_job_type_from_gjf

subroutine init_mol_comp_info(gjfname)
 implicit none
 integer, allocatable :: nuc(:)
 character(len=240), intent(in) :: gjfname

 call read_mem_and_nproc_from_gjf(gjfname, mem, nproc)
 call read_natom_from_gjf(gjfname, natom)
 allocate(elem(natom), nuc(natom), coor(3,natom))
 call read_elem_and_coor_from_gjf(gjfname, natom, elem, nuc, coor, charge, mult)
 deallocate(nuc)
 call get_job_type_from_gjf(gjfname)
end subroutine init_mol_comp_info
end module mol_comp_info

program geom_opt
 use mol_comp_info, only: init_mol_comp_info, gau_path
 use util_wrapper, only: gbw2mkl, mkl2fch_wrap
 implicit none
 integer :: i, fid
 character(len=240) :: gjfname, gjfname1, gjfname2, hfile, bakfile, numfile, &
  gbwname, gesname, engrad, mklname, dm_name, fchname1, fchname2, propfile

 i = iargc()
 if(i /= 1) then
  write(6,'(/,A)') ' ERROR in program geom_opt: wrong command line argument!'
  write(6,'(A,/)') ' Example: geom_opt h2o_wB97MV.gjf'
  stop
 end if

 gjfname = ' '
 call getarg(1, gjfname)
 call require_file_exist(gjfname)
 call find_specified_suffix(gjfname, '.gjf', i)
 gjfname2 = gjfname(1:i-1)//'_g.gjf'
 gjfname1 = gjfname(1:i-1)//'_o.gjf'
 gbwname = gjfname(1:i-1)//'_o.gbw'
 gesname = gjfname(1:i-1)//'_o.ges'
 engrad = gjfname(1:i-1)//'_o.engrad'
 mklname = gjfname(1:i-1)//'_o.mkl'
 dm_name = gjfname(1:i-1)//'_o.densities'
 fchname1 = gjfname(1:i-1)//'_o.fch'
 fchname2 = gjfname(1:i-1)//'_g.fch'
 hfile = gjfname(1:i-1)//'_o.nc'
 bakfile = gjfname(1:i-1)//'_o.bak'
 numfile = gjfname(1:i-1)//'_o.num'
 propfile = gjfname(1:i-1)//'_o_property.txt'

 call init_mol_comp_info(gjfname)
 call gen_guess_only_gjf(gjfname, gjfname1)
 call gen_orca_inp_gbw(gjfname1)
 call gen_gau_opt_gjf(gjfname2)

 ! Create an empty file. This file would be detected by the utility gau_external
 ! and it means that the 1st geometry need not be changed.
 open(newunit=fid,file=TRIM(hfile),status='replace')
 close(fid)

 ! Create an empty file. This file would be detected by the utility gau_external.
 ! It means that the Cartesian coordinates and the MOs of each new geometry would
 ! be saved into a corresponding .fch file. This takes some extra time for IO.
 open(newunit=fid,file=TRIM(bakfile),status='replace')
 close(fid)

 ! Create a file which contains only one integer. This integer indicates the
 ! number of steps.
 open(newunit=fid,file=TRIM(numfile),status='replace')
 write(fid,'(A1)') '0'
 close(fid)

 call submit_gau_job(gau_path, gjfname2, .false.)
 call delete_files(8, [gbwname, gesname, mklname, fchname1, bakfile, dm_name, &
                       engrad, propfile])
end program geom_opt

subroutine gen_guess_only_gjf(gjfname, gjfname1)
 use mol_comp_info, only: init_scf
 implicit none
 integer :: i, fid, fid1
 character(len=240) :: buf, chkname
 character(len=240), intent(in) :: gjfname, gjfname1

 i = INDEX(gjfname1, '.gjf', back=.true.)
 chkname = gjfname1(1:i-1)//'.chk'

 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(gjfname1),status='replace')

 write(fid1,'(A)') '%chk='//TRIM(chkname)

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:1) == '#') exit
  if(buf(1:4) == '%chk') cycle
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 call lower(buf)

 i = INDEX(buf, 'wb97m-v')
 if(i > 0) buf = buf(1:i-1)//'wB97X'//TRIM(buf(i+7:))

 call del_kywrd_in_buf(buf, 'opt')
 call del_kywrd_in_buf(buf, 'freq')

 buf = TRIM(buf)//' nosymm int=nobasistransform'
 if(.not. init_scf) buf = TRIM(buf)//' guess(only,save)'
 write(fid1,'(A)') TRIM(buf)

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for i

 close(fid)
 close(fid1)
end subroutine gen_guess_only_gjf

! delete a keyword (with its subkeywords, if any) in buf
! Note: do not include any space in key
subroutine del_kywrd_in_buf(buf, key)
 implicit none
 integer :: i, j, k1, k2
 character(len=*), intent(in) :: key
 character(len=240), intent(inout) :: buf

 i = INDEX(buf, key)
 if(i == 0) return
 ! now i > 0
 k1 = LEN_TRIM(key)
 k2 = LEN_TRIM(buf)
 j = i + k1 - 2 + INDEX(buf(i+k1-1:), ' ')
 if(j == k2+1) then
  buf(i:) = ' '
 else
  buf = buf(1:i-1)//TRIM(buf(j+1:))
 end if
end subroutine del_kywrd_in_buf

subroutine gen_orca_inp_gbw(gjfname)
 use mol_comp_info, only: mem, nproc, method, gau_path
 use util_wrapper, only: formchk, fch2mkl_wrap, mkl2gbw
 implicit none
 integer :: i, k, fid, fid1
 character(len=240) :: buf, proname, chkname, fchname, mklname, old_inp, &
  inpname, logname
 character(len=240), intent(in) :: gjfname

 call get_gau_path(gau_path)
 call submit_gau_job(gau_path, gjfname, .false.)

 i = INDEX(gjfname, '.gjf', back=.true.)
 if(i == 0) i = INDEX(gjfname, '.com', back=.true.)
 proname = gjfname(1:i-1)
 chkname = gjfname(1:i-1)//'.chk'
 fchname = gjfname(1:i-1)//'.fch'
 logname = gjfname(1:i-1)//'.log'
 mklname = gjfname(1:i-1)//'.mkl'
 old_inp = gjfname(1:i-1)//'_o.inp'
 inpname = gjfname(1:i-1)//'.inp'
 call formchk(chkname)
 call simplify_fch(fchname)
 call fch2mkl_wrap(fchname, mklname)
 call mkl2gbw(mklname)
 call delete_files(4, [chkname, gjfname, mklname, logname])

 open(newunit=fid,file=TRIM(old_inp),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname),status='replace')

 write(fid1,'(A,I0,A)') '%pal nprocs ', nproc, ' end'
 write(fid1,'(A,I0)') '%maxcore ', INT(DBLE(mem)/DBLE(nproc))

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:1) == '!') exit
 end do ! for while

 i = INDEX(buf, 'HF')
 k = INDEX(buf(1:i-1), ' ', back=.true.)
 buf = buf(1:k)//TRIM(method)//TRIM(buf(i+2:))

 i = INDEX(buf, 'VeryTightSCF')
 if(i > 0) buf = buf(1:i-1)//TRIM(buf(i+4:))
 buf = TRIM(buf)//' defgrid3 RIJCOSX def2/J EnGrad'
 write(fid1,'(A)') TRIM(buf)

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
end subroutine gen_orca_inp_gbw

subroutine gen_gau_opt_gjf(gjfname)
 use mol_comp_info, only: mem, charge, mult, natom, elem, coor, numfreq
 implicit none
 integer :: i, fid
 character(len=240), intent(in) :: gjfname

 open(newunit=fid,file=TRIM(gjfname),status='replace')
 write(fid,'(A,I0,A)') '%mem=', min(mem,4000), 'MB'
 write(fid,'(A)') '%nprocshared=1'
 write(fid,'(A)',advance='no') '# opt(nomicro,maxcycles=300)'
 if(numfreq) write(fid,'(A)',advance='no') ' freq=numer'
 write(fid,'(A)') " def2SVPP nosymm external='gau_external'"
 ! TODO: automatically switch to UGBS when elements are out of range of def2SVPP
 write(fid,'(/,A)') 'Using Gaussian as the geometry optimizer'
 write(fid,'(/,I0,1X,I0)') charge, mult
 do i = 1, natom, 1
  write(fid,'(A2,3(1X,F18.8))') elem(i), coor(:,i)
 end do ! for i
 write(fid,'(/)',advance='no')
 close(fid)
end subroutine gen_gau_opt_gjf

