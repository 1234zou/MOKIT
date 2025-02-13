! written by jxzou at 20210705

! generate natural orbitals(NO) using given density in .fch file
subroutine gen_no_using_density_in_fch(fchname, itype)
 implicit none
 integer :: i, j, nbf, nif
 integer, intent(in) :: itype ! 1/2 for Total/Spin SCF Density
!f2py intent(in) :: itype
 ! itype has values [1,10] in subroutine read_dm_from_fch
 ! here I use itype=0 to stand for alpha/beta natural spin orbitals
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname
 real(kind=8) :: na, nb, ne1, ne2
 real(kind=8), allocatable :: dm(:,:), spin_dm(:,:), dm_a(:,:), dm_b(:,:)
 real(kind=8), allocatable :: coeff_a(:,:), coeff_b(:,:), S(:,:)
 real(kind=8), allocatable :: noon_a(:), noon_b(:)

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 call read_na_and_nb_from_fch(fchname, i, j)
 na = DBLE(i); nb = DBLE(j)

 if(itype == 0) then
  allocate(dm(nbf,nbf), spin_dm(nbf,nbf))
  call read_dm_from_fch(fchname, 1, nbf, dm) ! Total SCF Density
  call read_dm_from_fch(fchname, 2, nbf, spin_dm) ! Spin SCF Density
  allocate(dm_a(nbf,nbf), dm_b(nbf,nbf))
  dm_a = 0.5d0*(dm + spin_dm)
  dm_b = 0.5d0*(dm - spin_dm)
  deallocate(dm, spin_dm)
 else
  allocate(dm_a(nbf,nbf))
  call read_dm_from_fch(fchname, itype, nbf, dm_a)
 end if

 allocate(S(nbf,nbf))
 call get_ao_ovlp_using_fch(fchname, nbf, S)
 select case(itype)
 case(0)
  call get_ne_from_PS(nbf, dm_a, S, ne1)
  call get_ne_from_PS(nbf, dm_b, S, ne2)
  call check_two_real8_eq(na, ne1, 1d-5)
  call check_two_real8_eq(nb, ne2, 1d-5)
 case(1,3,5,7,9)
  call get_ne_from_PS(nbf, dm_a, S, ne1)
  call check_two_real8_eq(na+nb, ne1, 1d-5)
 case default ! 2,4,6,8,10
  call get_ne_from_PS(nbf, dm_a, S, ne1)
  call check_two_real8_eq(na-nb, ne1, 1d-5)
 end select

 allocate(coeff_a(nbf,nif), noon_a(nif))
 call gen_no_from_density_and_ao_ovlp(nbf, nif, dm_a, S, noon_a, coeff_a)
 call write_eigenvalues_to_fch(fchname, nif, 'a', noon_a, .true.)
 call write_mo_into_fch(fchname, nbf, nif, 'a', coeff_a)
 deallocate(dm_a, coeff_a, noon_a)

 if(itype == 0) then
  allocate(coeff_b(nbf,nif), noon_b(nif))
  call gen_no_from_density_and_ao_ovlp(nbf, nif, dm_b, S, noon_b, coeff_b)
  call write_eigenvalues_to_fch(fchname, nif, 'b', noon_b, .true.)
  call write_mo_into_fch(fchname, nbf, nif, 'b', coeff_b)
  deallocate(dm_b, coeff_b, noon_b)
 end if

 deallocate(S)
end subroutine gen_no_using_density_in_fch

! get electronic dipole using specified density in a given .fch file
! Note: the nuclear dipole is calculated using subroutine get_nuc_dipole in
!       rwgeom.f90
subroutine get_e_dipole_using_density_in_fch(fchname, itype, dipole)
 implicit none
 integer :: nbf, nif
 integer, intent(in) :: itype
!f2py intent(in) :: itype
 ! itype has values [1,10] in subroutine read_density_from_fch
! real(kind=8) :: n_dipole(3)
 real(kind=8), intent(out) :: dipole(3)
!f2py intent(out) :: dipole
 real(kind=8), allocatable :: dm(:,:), D(:,:,:)
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)

 ! read desired AO-basis density matrix
 allocate(dm(nbf,nbf))
 call read_dm_from_fch(fchname, itype, nbf, dm)

 ! generate and read AO-basis dipole integrals
 allocate(D(nbf,nbf,3))
 call get_ao_dipole_using_fch(fchname, nbf, D)

 ! calculate electronic dipole moment
 call get_e_dipole_from_PD(nbf, dm, D, dipole)

 deallocate(dm, D)
end subroutine get_e_dipole_using_density_in_fch

! Read orthonormal basis (also called orthonormal atomic orbitals, OAO) from a
! given .fch(k) file.
subroutine read_orthonormal_basis_from_fch(fchname, nbf, mo, success)
 implicit none
 integer :: i, fid
 integer, intent(in) :: nbf
 real(kind=8), intent(out) :: mo(nbf,nbf)
 character(len=11) :: str
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
 logical, intent(out) :: success

 success = .false.; mo = 0d0
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  str = buf(1:11)
  select case(str)
  case('Orthonormal')
   exit
  case('Total SCF D','Mulliken Ch')
   i = -1
   exit
  end select
 end do ! for while

 if(i /= 0) then
  close(fid)
  return
 end if

 read(fid,'(5(1X,ES15.8))') mo
 close(fid)
 success = .true.
end subroutine read_orthonormal_basis_from_fch

! call Gaussian to compute AO-basis overlap integrals using the given .fch file
subroutine get_ao_ovlp_using_fch(fchname, nbf, S)
 implicit none
 integer :: nif
 integer, intent(in) :: nbf
 real(kind=8), intent(out) :: S(nbf,nbf)
 real(kind=8), allocatable :: mo(:,:)
 character(len=240), intent(in) :: fchname
 character(len=240) :: file47
 logical :: success

 call read_nif_from_fch(fchname, nif)

 if(nif == nbf) then ! no linear dependence
  allocate(mo(nbf,nbf))
  call read_orthonormal_basis_from_fch(fchname, nbf, mo, success)
  if(success) call solve_ovlp_from_cct(nbf, mo, S)
  deallocate(mo)
 end if

 if(nif<nbf .or. (nif==nbf .and. (.not.success))) then
  call call_gaussian_gen47_from_fch(fchname, file47)
  call read_ao_ovlp_from_47(file47, nbf, S)
  call delete_file(TRIM(file47))
 end if
end subroutine get_ao_ovlp_using_fch

! call Gaussian to compute AO-basis dipole integrals using the given .fch file
subroutine get_ao_dipole_using_fch(fchname, nbf, D)
 use phys_cons, only: Bohr_const
 implicit none
 integer, intent(in) :: nbf
 real(kind=8), intent(out) :: D(nbf,nbf,3)
 character(len=240), intent(in) :: fchname
 character(len=240) :: file47

 call call_gaussian_gen47_from_fch(fchname, file47)
 call read_ao_dipole_from_47(file47, nbf, D)
 call delete_file(file47)
 D = D/Bohr_const   ! important
end subroutine get_ao_dipole_using_fch

! call Gaussian program to generate .47 file from a given .fch(k) file
subroutine call_gaussian_gen47_from_fch(fchname, file47)
 use util_wrapper, only: unfchk
 implicit none
 integer :: i, k, fid, system
 character(len=10) :: str
 character(len=24) :: mem
 character(len=240), intent(in) :: fchname
 character(len=240), intent(out) :: file47
 character(len=240) :: gau_path, proname, chkname, gjfname, logname

 call get_a_random_int(i)
 write(str,'(I0)') i

 call find_specified_suffix(fchname, '.fch', i)
 proname = fchname(1:i-1)//TRIM(str)
 chkname = TRIM(proname)//'.chk'
 gjfname = TRIM(proname)//'.gjf'
 file47  = TRIM(proname)//'.47'
 call upper(file47)
#ifdef _WIN32
 logname = TRIM(proname)//'.out'
#else
 logname = TRIM(proname)//'.log'
#endif

 call getenv('GAUSS_MEMDEF', mem)
 if(LEN_TRIM(mem) == 0) then
  mem = '1GB'
 else
  mem = ADJUSTL(mem)
  k = LEN_TRIM(mem)
  if(mem(k-1:k-1) == 'G') then
   read(mem(1:k-2),*) i
   if(i > 4) mem = '4GB'
  end if
 end if

 open(newunit=fid,file=TRIM(gjfname),status='replace')
 write(fid,'(A)') '%nprocshared=1'
 write(fid,'(A)') '%mem='//TRIM(mem)
 write(fid,'(A)') '%chk='//TRIM(chkname)
 write(fid,'(A)') '# chkbasis nosymm int=nobasistransform guess(read,only)&
                 & geom=allcheck pop(nboread)'
 write(fid,'(/,A)') '$NBO'
 write(fid,'(A)')   ' NOBOND'
 write(fid,'(A)')   ' SKIPBO'
 write(fid,'(A)')   ' archive file='//TRIM(proname)
 write(fid,'(A,/)') '$END'
 close(fid)
 call unfchk(fchname, chkname)

 call get_gau_path(gau_path)
 i = SYSTEM(TRIM(gau_path)//' '//TRIM(gjfname))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine call_gaussian_gen47_from_fch: Gaussian&
                  & job failed.'
  write(6,'(A)') 'You can open file '//TRIM(logname)//' and check why.'
  stop
 end if

 call delete_files(3, [chkname, gjfname, logname])
end subroutine call_gaussian_gen47_from_fch

! read AO-basis dipole integrals from a given .47 file
subroutine read_ao_dipole_from_47(file47, nbf, D)
 implicit none
 integer :: i, j, k, fid
 integer, intent(in) :: nbf
 real(kind=8), intent(out) :: D(nbf,nbf,3)
 character(len=240) :: buf
 character(len=240), intent(in) :: file47

 D = 0d0
 open(newunit=fid,file=TRIM(file47),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:8) == '$DIPOLE') exit
 end do ! for while

 if(i /= 0) then
  close(fid)
  write(6,'(/,A)') "ERROR in subroutine read_ao_dipole_from_47: no '$DIPOLE'&
                  & found in file "//TRIM(file47)
  stop
 end if

 do i = 1, 3
  read(fid,'(2X,5E15.7)') ((D(k,j,i),k=1,j),j=1,nbf)
 end do ! for i
 close(fid)

 ! remember to symmetrize the dipole integral matrix
 forall(i=1:3, j=1:nbf-1, k=1:nbf, j<k) D(k,j,i) = D(j,k,i)
end subroutine read_ao_dipole_from_47

! calculate the electronic dipole moment from AO-basis density matrix and
! dipole integrals
subroutine get_e_dipole_from_PD(nbf, P, D, e_dipole)
 implicit none
 integer :: i, j
 integer, intent(in) :: nbf
 real(kind=8), intent(in) :: P(nbf,nbf), D(nbf,nbf,3)
 real(kind=8), intent(out) :: e_dipole(3) ! x,y,z 3-components
 real(kind=8), allocatable :: r(:)

 allocate(r(nbf))

 do i = 1, 3
  forall(j = 1:nbf) r(j) = DOT_PRODUCT(P(:,j), D(:,j,i))
  e_dipole(i) = -SUM(r)
 end do ! for i

 deallocate(r)
end subroutine get_e_dipole_from_PD

! read the path of the Gaussian binary executable file 
subroutine get_gau_path(gau_path)
 implicit none
 integer :: i
 character(len=240) :: g(3)
 character(len=240), intent(out) :: gau_path
 character(len=480) :: buf
 logical :: alive

 gau_path = ' '
 call getenv('GAUSS_EXEDIR', buf)
 ! in case that the $GAUSS_EXEDIR is too long, use buf to store it

 if(LEN_TRIM(buf) == 0) then
  gau_path = 'NOT FOUND'
  return
 end if

#ifdef _WIN32
 i = INDEX(buf, '\', back=.true.)
 if(i == 0) then
  write(6,'(/,A)') "ERROR in subroutine get_gau_path: no '\' symbol found in ga&
                   &u_path="//TRIM(buf)
  stop
 end if
 buf = """"//TRIM(buf)//'\g'//buf(i+2:i+3)//".exe"""
 if(LEN_TRIM(buf) > 240) then
  write(6,'(/,A)') 'ERROR in subroutine get_gau_path: Gaussian path is too long!'
  write(6,'(A)') 'Please install your Gaussian program in a shorter path.'
  stop
 end if
 gau_path = TRIM(buf)

#else
 i = INDEX(buf, ':', back=.true.)
 if(i == 0) then
  write(6,'(/,A)') "ERROR in subroutine get_gau_path: no ':' symbol found in g&
                   &au_path="//TRIM(buf)
  write(6,'(/,A)') 'This error often occurs when your machine has no (or has i&
                   &ncorrect) Gaussian'
  write(6,'(A)') 'environment variables. Here a correct example is shown:'
  write(6,'(A)') REPEAT('-',45)
  write(6,'(A)') ' export g16root=/opt'
  write(6,'(A)') ' source $g16root/g16/bsd/g16.profile'
  write(6,'(A)') ' export GAUSS_SCRDIR=/scratch/$USER/gaussian'
  write(6,'(A)') REPEAT('-',45)
  write(6,'(A)') 'Please check your Gaussian environment variables according to&
                 & the example shown'
  write(6,'(A)') 'above.'
  write(6,'(A,/)') "Also note: DO NOT write 'export GAUSS_EXEDIR', it is useless."
  stop
 end if

 buf = buf(i+1:)
 if(LEN_TRIM(buf) > 240) then
  write(6,'(/,A)') 'ERROR in subroutine get_gau_path: Gaussian path is too long!'
  write(6,'(A)') 'Please install your Gaussian program in a shorter path.'
  stop
 end if
 gau_path = TRIM(buf)
 g(1) = TRIM(gau_path)//'/g03'
 g(2) = TRIM(gau_path)//'/g09'
 g(3) = TRIM(gau_path)//'/g16'
 ! maybe extended to g23 in 2023

 alive = .false.
 do i = 3, 1, -1
  inquire(file=TRIM(g(i)),exist=alive)
  if(alive) then
   gau_path = g(i)
   exit
  end if
 end do ! for i

 if(.not. alive) then
  write(6,'(/,A)') 'ERROR in subroutine get_gau_path: Gaussian is not found in&
                  & your node/machine.'
  write(6,'(A)') 'Please check whether you have installed Gaussian correctly.'
  stop
 end if
#endif
end subroutine get_gau_path

! find/get the absolute path of a specified binary executable
subroutine get_exe_path(exe_name, exe_path)
 implicit none
 integer :: i, fid, SYSTEM
 character(len=30) :: fname
 character(len=*), intent(in) :: exe_name
 character(len=240), intent(out) :: exe_path

 exe_path = ' '
 call get_a_random_int(i)
 write(fname,'(A,I0)') 'mokit.'//TRIM(exe_name)//'.', i

 i = SYSTEM('which '//TRIM(exe_name)//' >'//TRIM(fname)//" 2>&1")
 if(i /= 0) then
  call delete_file(TRIM(fname))
  exe_path = 'NOT FOUND'
  return
 end if

 open(newunit=fid,file=TRIM(fname),status='old',position='rewind')
 read(fid,'(A)',iostat=i) exe_path
 close(fid,status='delete')

 if(i /= 0) then
  exe_path = 'NOT FOUND'
 else
  if(LEN_TRIM(exe_path) == 0) then
   exe_path = 'NOT FOUND'
  else if(INDEX(exe_path,'no '//TRIM(exe_name)) > 0) then
   exe_path = 'NOT FOUND'
  end if
 end if
end subroutine get_exe_path

subroutine get_psi4_path(psi4_path)
 implicit none
 character(len=240), intent(out) :: psi4_path

 psi4_path = ' '
 call getenv('PSI4', psi4_path)
 if(LEN_TRIM(psi4_path) > 0) return
 ! $PSI4 has higher priority than 'which psi4'

 call get_exe_path('psi4', psi4_path)
end subroutine get_psi4_path

! compute the number of electrons by tracing the product of density matrix and
! AO-basis overlap
subroutine get_ne_from_PS(nbf, P, S, ne)
 implicit none
 integer :: i
 integer, intent(in) :: nbf
 real(kind=8), allocatable :: ne0(:)
 real(kind=8), intent(in) :: P(nbf,nbf), S(nbf,nbf)
 real(kind=8), intent(out) :: ne

 allocate(ne0(nbf))
 forall(i = 1:nbf) ne0(i) = DOT_PRODUCT(P(:,i),S(:,i))
 ne = SUM(ne0)
 deallocate(ne0)
end subroutine get_ne_from_PS

! calculate the total number of electrons using total density in a .fch file
subroutine get_ne_from_fch(fchname)
 implicit none
 integer :: nbf, nif
 real(kind=8) :: ne
 real(kind=8), allocatable :: den(:,:), S(:,:)
 character(len=240), intent(in) :: fchname

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(den(nbf,nbf), S(nbf,nbf))
 call read_dm_from_fch(fchname, 1, nbf, den)
 call get_ao_ovlp_using_fch(fchname, nbf, S)
 call get_ne_from_PS(nbf, den, S, ne)
 deallocate(den, S)
 write(6,'(A,F11.4)') 'ne = ', ne
end subroutine get_ne_from_fch

! check whether two double precision values equal to each other
subroutine check_two_real8_eq(r1, r2, thres)
 implicit none
 real(kind=8), intent(in) :: r1, r2, thres

 if(thres < 0d0) then
  write(6,'(/,A)') 'ERROR in subroutine check_two_real8_eq: input thres<0.'
  write(6,'(A,F14.8,A)') 'thres=',thres,', but thres>0 is required.'
  stop
 end if

 if(DABS(r1-r2) > thres) then
  write(6,'(/,A)') 'ERROR in subroutine check_two_real8_eq: two double precisi&
                    &on values differ larger'
  write(6,'(A,3F14.8)') 'than the threshold. r1, r2, thres=', r1, r2, thres
  stop
 end if
end subroutine check_two_real8_eq

! submit a GAMESS frozen core GVB job
subroutine submit_gms_fzgvb_job(gms_path, gms_scr_path, inpname, nproc)
 implicit none
 integer :: i, fid, fid1
 integer, intent(in) :: nproc
 character(len=240) :: buf, inpname1, datname1
 character(len=240), intent(in) :: gms_path, gms_scr_path, inpname

 call find_specified_suffix(inpname, '.inp', i)
 inpname1 = inpname(1:i-1)//'_f.inp'
 datname1 = inpname(1:i-1)//'_f.dat'
 call submit_gms_job(gms_path, gms_scr_path, inpname1, nproc)

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(datname1),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:5) == '$VEC') exit
 end do ! for while

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(2:5) == '$VEC') exit
 end do ! for while

 do while(.true.)
  read(fid1,'(A)') buf
  write(fid,'(A)') TRIM(buf)
  if(buf(2:5) == '$END') exit
 end do ! for while

 close(fid)
 close(fid1)
end subroutine submit_gms_fzgvb_job

! submit a GAMESS job
subroutine submit_gms_job(gms_path, gms_scr_path, inpname, nproc)
 implicit none
 integer :: i, system
 integer, intent(in) :: nproc
 character(len=240) :: datname, gmsname, hs1, hs2, trj
 character(len=500) :: longbuf
 character(len=240), intent(in) :: gms_path, gms_scr_path, inpname

 call find_specified_suffix(inpname, '.inp', i)

 ! delete scratch files, if any
 datname = inpname(1:i-1)//'.dat'
 gmsname = inpname(1:i-1)//'.gms'
 hs1 = inpname(1:i-1)//'.hs1'
 hs2 = inpname(1:i-1)//'.hs2'
 trj = inpname(1:i-1)//'.trj'
 call delete_files_in_path(gms_scr_path, 4, [datname, hs1, hs2, trj])

 ! now we can submit the GAMESS job
 write(longbuf,'(A,I0,A)') TRIM(inpname)//' 01 ',nproc,' >'//TRIM(gmsname)//" 2>&1"
 write(6,'(A)') '$$GMS '//TRIM(longbuf)

 i = SYSTEM(TRIM(gms_path)//' '//TRIM(longbuf))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine submit_gms_job: GAMESS job failed.'
  write(6,'(A)') 'You can open file '//TRIM(gmsname)//' and check why.'
  stop
 end if

 ! move the .dat file into current directory
 i = SYSTEM('mv '//TRIM(gms_scr_path)//'/'//TRIM(datname)//' .')
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine submit_gms_job: fail to move file. Poss&
                   &ibly wrong gms_scr_path.'
  write(6,'(A)') 'gms_scr_path='//TRIM(gms_scr_path)
  stop
 end if
end subroutine submit_gms_job

! submit a MOKIT automr job
subroutine submit_automr_job(gjfname)
 implicit none
 integer :: i, system
 character(len=240) :: buf, outname
 character(len=240), intent(in) :: gjfname

 call find_specified_suffix(gjfname, '.gjf', i)
 outname = gjfname(1:i-1)//'.out'
 buf = 'automr '//TRIM(gjfname)//' >'//TRIM(outname)//" 2>&1"
 write(6,'(A)') '$'//TRIM(buf)

 i = SYSTEM(TRIM(buf))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine submit_automr_job: automr job failed.'
  write(6,'(A)') 'You can open file '//TRIM(outname)//' and check why.'
  stop
 end if
end subroutine submit_automr_job

subroutine submit_gau_job(gau_path, gjfname, prt)
 implicit none
 integer :: i, system
 character(len=240) :: logname
 character(len=240), intent(in) :: gau_path, gjfname
 logical, intent(in) :: prt

 call find_specified_suffix(gjfname, '.gjf', i)
#ifdef _WIN32
 logname = gjfname(1:i-1)//'.out'
#else
 logname = gjfname(1:i-1)//'.log'
#endif

 if(prt) write(6,'(A)') '$'//TRIM(gau_path)//' '//TRIM(gjfname)
 i = SYSTEM(TRIM(gau_path)//' '//TRIM(gjfname))

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subrouitine submit_gau_job: Gaussian job failed.'
  write(6,'(A)') 'Please open file '//TRIM(logname)//' and check.'
  stop
 end if

 call delete_file('fort.7')
end subroutine submit_gau_job

subroutine submit_orca_job(orca_path, inpname, prt, del_den, del_prop)
 implicit none
 integer :: i, system
 character(len=240) :: outname, gesname, denf1, denf2, propf1, propf2
 character(len=480) :: buf
 character(len=240), intent(in) :: orca_path, inpname
 logical, intent(in) :: prt, del_den, del_prop

 call find_specified_suffix(inpname, '.inp', i)
 outname = inpname(1:i-1)//'.out'
 gesname = inpname(1:i-1)//'.ges'
 denf1 = inpname(1:i-1)//'.densities'
 denf2 = inpname(1:i-1)//'.densitiesinfo'
 propf1 = inpname(1:i-1)//'_property.txt'
 propf2 = inpname(1:i-1)//'.property.txt'

 write(buf,'(A)') TRIM(inpname)//' >'//TRIM(outname)//" 2>&1"
 if(prt) write(6,'(A)') '$orca '//TRIM(buf)

 i = SYSTEM(TRIM(orca_path)//' '//TRIM(buf))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subrouitine submit_orca_job: ORCA job failed.'
  write(6,'(A)') 'Please open file '//TRIM(outname)//' and check.'
  stop
 end if

 call delete_file(TRIM(gesname))
 if(del_den) call delete_files(2, [denf1,denf2])
 if(del_prop) call delete_files(2, [propf2,propf2])
end subroutine submit_orca_job

! need to distinguish between OpenMP version and MPI version of OpenMolcas
subroutine submit_molcas_job(inpname, mem, nproc, omp)
 implicit none
 integer :: i, fid, system
 integer, intent(in) :: mem, nproc ! mem in GB
 character(len=240) :: shname, statname, gssorb, gss_h5, gss_molden, outname
 character(len=500) :: buf
 character(len=240), intent(in) :: inpname
 logical, intent(in) :: omp

 call find_specified_suffix(inpname, '.inp', i)
 shname = inpname(1:i-1)//'.sh'
 statname = inpname(1:i-1)//'.status'
 gssorb = inpname(1:i-1)//'.GssOrb'
 gss_h5 = inpname(1:i-1)//'.guessorb.h5'
 gss_molden = inpname(1:i-1)//'.guessorb.molden'
 outname = inpname(1:i-1)//'.out'

 open(newunit=fid,file=TRIM(shname),status='replace')
 if(omp) then ! OpenMP version
  write(fid,'(A)') '#OpenMP'
  write(fid,'(A,I0,A)') 'export MOLCAS_MEM=',mem,'Gb'
  write(fid,'(A)') 'export MOLCAS_NPROCS=1'
  write(fid,'(A,I0)') 'export OMP_NUM_THREADS=',nproc
  write(buf,'(A,I0)') 'pymolcas -nt ', nproc
 else ! MPI version
  write(fid,'(A)') '#MPI'
  write(fid,'(A,I0)') 'export MOLCAS_MEM=',INT(DBLE(mem)*1d3/DBLE(nproc))
  write(fid,'(A,I0)') 'export MOLCAS_NPROCS=',nproc
  write(fid,'(A)') 'export OMP_NUM_THREADS=1'
  write(buf,'(A,I0)') 'pymolcas -nt 1 -np ', nproc
 end if

 buf = TRIM(buf)//' '//TRIM(inpname)//' >'//TRIM(outname)//" 2>&1"
 write(6,'(A)') '$'//TRIM(buf)
 write(fid,'(A)') TRIM(buf)
 close(fid)

 i = SYSTEM('/bin/bash '//TRIM(shname))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subrouitine submit_molcas_job: OpenMolcas job fail&
                   &ed. Please'
  write(6,'(A)') 'open file '//TRIM(outname)//' and check.'
  stop
 end if

 call delete_files(5, [shname, statname, gssorb, gss_h5, gss_molden])
end subroutine submit_molcas_job

! some default settings are different among different versions of PSI4, so we
! have to find its version
subroutine get_psi4_version(version)
 implicit none
 integer :: i, fid, SYSTEM
 character(len=10), intent(out) :: version
 character(len=30) :: fname
 character(len=240) :: psi4_path

 version = ' '; psi4_path = ' '
 call get_a_random_int(i)
 write(fname,'(A,I0)') 'psi4.ver.', i

 call getenv('PSI4', psi4_path)
 if(LEN_TRIM(psi4_path) > 0) then
  i = SYSTEM('$PSI4 --version >'//TRIM(fname)//" 2>&1")
 else
  i = SYSTEM('psi4 --version >'//TRIM(fname)//" 2>&1")
 end if

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine get_psi4_version: failed to find PSI4 v&
                   &ersion.'
  write(6,'(A)') 'Please check whether your PSI4 has been installed correctly.'
  call delete_file(TRIM(fname))
  stop
 end if

 open(newunit=fid,file=TRIM(fname),status='old',position='rewind')
 read(fid,'(A)') version
 close(fid,status='delete')
end subroutine get_psi4_version

subroutine submit_psi4_job(psi4_path, inpname, nproc)
 implicit none
 integer :: i, system
 integer, intent(in) :: nproc
 character(len=240) :: outname
 character(len=480) :: buf
 character(len=240), intent(in) :: psi4_path, inpname

 call find_specified_suffix(inpname, '.inp', i)
 outname = inpname(1:i-1)//'.out'

 write(buf,'(A,I0)') TRIM(inpname)//' '//TRIM(outname)//' -n ', nproc
 write(6,'(A)') '$'//TRIM(psi4_path)//' '//TRIM(buf)

 i = SYSTEM(TRIM(psi4_path)//' '//TRIM(buf))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subrouitine submit_psi4_job: PSI4 job failed.'
  write(6,'(A)') 'Please open file '//TRIM(outname)//' and check.'
  stop
 end if

 call delete_file('timer.dat')
 call delete_file('ijk.dat')
end subroutine submit_psi4_job

subroutine submit_qchem_job(inpname, nproc)
 implicit none
 integer :: i, system
 integer, intent(in) :: nproc
 character(len=240) :: scr_dir, outname
 character(len=500) :: buf
 character(len=240), intent(in) :: inpname

 call find_specified_suffix(inpname, '.in', i)
 scr_dir = inpname(1:i-1)
 outname = inpname(1:i-1)//'.out'

 write(buf,'(A,I0,A)') 'qchem -nt ',nproc,' -np 1 '//TRIM(inpname)//' '//&
                       TRIM(outname)//' '//TRIM(scr_dir)//' >junk 2>&1'
 write(6,'(A)') '$'//TRIM(buf)

 i = SYSTEM(TRIM(buf))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subrouitine submit_qchem_job: QChem job failed.'
  write(6,'(A)') 'Please open file '//TRIM(outname)//' and check.'
  stop
 end if

 call delete_file('junk')
 call delete_file('pathtable')
end subroutine submit_qchem_job

subroutine submit_molpro_job(inpname, mem, nproc)
 implicit none
 integer :: i, SYSTEM
 integer, intent(in) :: mem ! total memory in GB
 integer, intent(in) :: nproc ! number of processors
 character(len=270) :: buf = ' '
 character(len=240) :: outname, xmlname
 character(len=240), intent(in) :: inpname

 call find_specified_suffix(inpname, '.com', i)
 outname = inpname(1:i-1)//'.out'
 xmlname = inpname(1:i-1)//'.xml'
 call delete_files(2, [outname, xmlname]) ! clean old files (if any)

 i = CEILING(DBLE(mem*125)/DBLE(nproc))
 write(buf,'(2(A,I0),A)') 'molpro -W ./ -t 1 -n ',nproc,' -m ',i,'m '//TRIM(inpname)
 write(6,'(A)') '$'//TRIM(buf)

 i = SYSTEM(TRIM(buf))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine submit_molpro_job: Molpro job failed.'
  write(6,'(A)') 'You can open the file '//TRIM(outname)//' and check why.'
  stop
 end if
end subroutine submit_molpro_job

! submit a GVB-BCCI job
subroutine submit_gvb_bcci_job(nproc, ci_order, inpname, outname)
 implicit none
 integer :: i, fid, system
 integer, intent(in) :: nproc, ci_order
 character(len=240) :: shname
 character(len=240), intent(in) :: inpname, outname

 if(ci_order /= 2) then
  write(6,'(A)') 'ERROR in subroutine submit_gvb_bcci_job: only CI_order=2 is&
                & allowed.'
  write(6,'(A,I0)') 'Input CI_order=', ci_order
  stop
 end if

 call find_specified_suffix(inpname, '.inp', i)
 shname = inpname(1:i-1)//'.sh'
 open(newunit=fid,file=TRIM(shname),status='replace')
 write(fid,'(A)') 'export OMP_STACKSIZE=2G'
 write(fid,'(A,I0)') 'export OMP_NUM_THREADS=',nproc
 write(fid,'(A)') 'gvb_bcci2b '//TRIM(inpname)//' 2 >'//TRIM(outname)//" 2>&1"
 close(fid)

 i = SYSTEM('/bin/bash '//TRIM(shname))
 call delete_file(shname)
 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine submit_gvb_bcci_job: Linearized BCCC job&
                & failed.'
  write(6,'(A)') 'You can open file '//TRIM(outname)//' and check why.'
  stop
 end if
end subroutine submit_gvb_bcci_job

! submit a GVB-BCCC job
subroutine submit_gvb_bccc_job(mult, nproc, cc_order, inpname, outname)
 implicit none
 integer :: i, fid, system
 integer, intent(in) :: mult, nproc, cc_order
 character(len=240) :: bccc_prog, shname
 character(len=240), intent(in) :: inpname, outname

 select case(cc_order)
 case(2)
  if(mult == 1) then
   bccc_prog = 'gvb_bccc2b'
  else
   bccc_prog = 'gvb_bccc2b_T'
  end if
 case(3)
  bccc_prog = 'gvb_bccc3b'
 case default
  write(6,'(A)') 'ERROR in subroutine submit_gvb_bccc_job: only CC_order=2 or&
                & 3 is allowed.'
  write(6,'(A,I0)') 'Input CC_order=', cc_order
  stop
 end select

 call find_specified_suffix(inpname, '.inp', i)
 shname = inpname(1:i-1)//'.sh'
 open(newunit=fid,file=TRIM(shname),status='replace')
 write(fid,'(A)') 'export OMP_STACKSIZE=2G'
 write(fid,'(A,I0)') 'export OMP_NUM_THREADS=',nproc
 write(fid,'(A)') TRIM(bccc_prog)//' '//TRIM(inpname)//' >'//TRIM(outname)//&
                  " 2>&1"
 close(fid)

 i = SYSTEM('/bin/bash '//TRIM(shname))
 call delete_file(shname)
 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine submit_gvb_bccc_job: Linearized BCCC job&
                & failed.'
  write(6,'(A)') 'You can open file '//TRIM(outname)//' and check why.'
  stop
 end if
end subroutine submit_gvb_bccc_job

subroutine submit_pyscf_job(pyname, prt)
 implicit none
 integer :: i, SYSTEM
 character(len=240) :: outname
 character(len=240), intent(in) :: pyname
 character(len=500) :: buf
 logical, intent(in) :: prt

 call find_specified_suffix(pyname, '.py', i)
 outname = pyname(1:i-1)//'.out'

 write(buf,'(A)') 'python '//TRIM(pyname)//' >'//TRIM(outname)//" 2>&1"
 if(prt) write(6,'(A)') '$'//TRIM(buf)
 i = SYSTEM(TRIM(buf))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subrouitine submit_pyscf_job: PySCF job failed.'
  write(6,'(A)') 'Please open file '//TRIM(outname)//' and check.'
  stop
 end if
end subroutine submit_pyscf_job

! Submit a CFOUR job. Note: the CFOUR input file is always ZMAT, so there is no
!  need to specify an input name, but specifying an outname is required here
subroutine submit_cfour_job(nproc, outname, omp)
 implicit none
 integer :: i, fid, SYSTEM
 integer, intent(in) :: nproc
 character(len=20) :: shname
 character(len=240), intent(in) :: outname
 logical, intent(in) :: omp

 call get_a_random_int(i)
 write(shname,'(I0,A)') i, '.sh'
 open(newunit=fid,file=TRIM(shname),status='replace')

 if(omp) then ! OpenMP
  write(fid,'(A)') 'export CFOUR_NUM_CORES=1'
  write(fid,'(A,I0)') 'export OMP_NUM_THREADS=', nproc
 else         ! MPI
  write(fid,'(A,I0)') 'export CFOUR_NUM_CORES=', nproc
  write(fid,'(A)') 'export OMP_NUM_THREADS=1'
 end if
 write(fid,'(A)') "xcfour >"//TRIM(outname)//" 2>&1"
 close(fid)

 i = SYSTEM('/bin/bash '//TRIM(shname))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine submit_cfour_job: CFOUR job failed.'
  write(6,'(A)') 'Please open file '//TRIM(outname)//' and check.'
  stop
 else
  call delete_file(TRIM(shname))
  call delete_cfour_tmp_files(nproc)
 end if
end subroutine submit_cfour_job

subroutine submit_dalton_job(proname, mem, nproc, mpi, sirius, noarch, del_sout)
 implicit none
 integer :: i, fid, SYSTEM
 integer, intent(in) :: mem, nproc
 character(len=20) :: shname
 character(len=240) :: sout
 character(len=240), intent(in) :: proname
 character(len=500) :: buf
 logical, intent(in) :: mpi, sirius, noarch, del_sout

 if(mpi) then ! MPI
  write(buf,'(2(A,I0))') 'dalton -gb ',MIN(mem,16),' -N ',nproc
 else
  write(buf,'(2(A,I0))') 'dalton -gb ',MIN(mem,16),' -omp ',nproc
 end if
 if(sirius) buf = TRIM(buf)//" -put ""SIRIUS.RST"""
 if(noarch) buf = TRIM(buf)//' -noarch'
 sout = TRIM(proname)//'.sout'
 buf = TRIM(buf)//' -ow '//TRIM(proname)//' >'//TRIM(sout)//" 2>&1"
 write(6,'(A)') '$'//TRIM(buf)

 if(mpi) then
  call get_a_random_int(i)
  write(shname,'(I0,A)') i, '.sh'
  open(newunit=fid,file=TRIM(shname),status='replace')
  write(fid,'(A)') 'unset DALTON_LAUNCHER'
  write(fid,'(A)') TRIM(buf)
  close(fid)
  buf = '/bin/bash '//TRIM(shname)
 end if

 i = SYSTEM(TRIM(buf))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine do_cas: Dalton job failed.'
  write(6,'(A)') 'Please open files '//TRIM(proname)//'.out as well as '//&
                 TRIM(sout)//' and check why.'
  stop
 end if

 if(mpi) call delete_file(TRIM(shname))
 if(del_sout) call delete_file(TRIM(sout))
 if(.not. noarch) i = SYSTEM('tar -xpf '//TRIM(proname)//'.tar.gz SIRIUS.RST')
end subroutine submit_dalton_job

! check/detect OpenMolcas is OpenMP version or MPI version
subroutine check_molcas_is_omp(molcas_omp)
 implicit none
 integer :: i, fid, SYSTEM
 character(len=3) :: str
 character(len=30) :: fname
 character(len=240) :: buf
 logical, intent(out) :: molcas_omp

 str = ' '; molcas_omp = .true.
 call get_a_random_int(i)
 write(fname,'(A,I0)') 'molcas.ver.', i
 i = SYSTEM('pymolcas --banner >'//TRIM(fname)//" 2>&1")

 ! maybe OpenMolcas not installed, assume OpenMP version
 if(i /= 0) then
  call delete_file(TRIM(fname))
  return
 end if

 ! OK, now it is installed
 open(newunit=fid,file=TRIM(fname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:8) == 'Parallel') exit
 end do ! for while

 ! maybe 'bash: pymolcas: command not found', assume OpenMP version
 if(i /= 0) then
  close(fid,status='delete')
  return
 end if

 i = INDEX(buf,':')
 read(buf(i+1:),*) str
 close(fid,status='delete')
 str = ADJUSTL(str)

 select case(TRIM(str))
 case('ON')  ! MPI version
  molcas_omp = .false.
 case('OFF') ! OpenMP version
 case default
  write(6,'(/,A)') 'ERROR in subroutine check_molcas_is_omp: unrecognized paral&
                   &lel type='//str
  stop
 end select
end subroutine check_molcas_is_omp

! check/detect Dalton is MKL version or MPI version
subroutine check_dalton_is_mpi(mpi)
 implicit none
 integer :: i, fid, SYSTEM
 character(len=30) :: fname
 character(len=240) :: buf
 logical, intent(out) :: mpi

 mpi = .false.
 call get_a_random_int(i)
 write(fname,'(A,I0)') 'dalton.ver.', i
 i = SYSTEM('ldd `(which dalton.x)` >'//TRIM(fname)//" 2>&1")

 ! maybe Dalton not installed, assume MKL version
 if(i /= 0) then
  call delete_file(TRIM(fname))
  return
 end if

 ! OK, now it is installed
 open(newunit=fid,file=TRIM(fname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(INDEX(buf,'libmpi') > 0) then
   mpi = .true.
   exit
  end if
 end do ! for while

 close(fid,status='delete')
end subroutine check_dalton_is_mpi

