! written by jxzou at 20210705

! generate natural orbitals(NO) using given density in .fch file
subroutine gen_no_using_density_in_fch(fchname, itype)
 implicit none
 integer :: i, j, nbf, nif
 integer, intent(in) :: itype ! 1/2 for Total/Spin SCF Density
 ! itype has values [1,10] in subroutine read_density_from_fch
 ! here I use itype=0 to stand for alpha/beta natural spin orbitals
 character(len=240), intent(in) :: fchname
 real(kind=8) :: na, nb, ne1, ne2
 real(kind=8), allocatable :: dm(:,:), spin_dm(:,:), dm_a(:,:), dm_b(:,:)
 real(kind=8), allocatable :: coeff_a(:,:), coeff_b(:,:), S(:,:)
 real(kind=8), allocatable :: noon_a(:), noon_b(:)

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 call read_na_and_nb_from_fch(fchname, i, j)
 na = DBLE(i); nb = DBLE(j)

 if(itype == 0) then
  allocate(dm(nbf,nbf), spin_dm(nbf,nbf))
  call read_density_from_fch(fchname, 1, nbf, dm) ! Total SCF Density
  call read_density_from_fch(fchname, 2, nbf, spin_dm) ! Spin SCF Density
  dm_a = 0.5d0*(dm + spin_dm) ! auto-allocation
  dm_b = 0.5d0*(dm - spin_dm) ! auto-allocation
  deallocate(dm, spin_dm)
 else
  allocate(dm_a(nbf,nbf))
  call read_density_from_fch(fchname, itype, nbf, dm_a)
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
 call get_no_from_density_and_ao_ovlp(nbf, nif, dm_a, S, noon_a, coeff_a)
 call write_eigenvalues_to_fch(fchname, nif, 'a', noon_a, .true.)
 call write_mo_into_fch(fchname, nbf, nif, 'a', coeff_a)
 deallocate(dm_a, coeff_a, noon_a)

 if(itype == 0) then
  allocate(coeff_b(nbf,nif), noon_b(nif))
  call get_no_from_density_and_ao_ovlp(nbf, nif, dm_b, S, noon_b, coeff_b)
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
 integer :: i, nbf, nif
 integer, intent(in) :: itype
!f2py intent(in) :: itype
 ! itype has values [1,10] in subroutine read_density_from_fch
 real(kind=8) :: n_dipole(3)
 real(kind=8), intent(out) :: dipole(3)
!f2py intent(out) :: dipole
 real(kind=8), allocatable :: dm(:,:), D(:,:,:)
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)

 ! read desired AO-basis density matrix
 allocate(dm(nbf,nbf))
 call read_density_from_fch(fchname, itype, nbf, dm)

 ! generate and read AO-basis dipole integrals
 allocate(D(nbf,nbf,3))
 call get_ao_dipole_using_fch(fchname, nbf, D)

 ! calculate electronic dipole moment
 call get_e_dipole_from_PD(nbf, dm, D, dipole)

 deallocate(dm, D)
end subroutine get_e_dipole_using_density_in_fch

! call Gaussian to compute AO-basis overlap integrals using the given .fch file
subroutine get_ao_ovlp_using_fch(fchname, nbf, S)
 implicit none
 integer, intent(in) :: nbf
 real(kind=8), intent(out) :: S(nbf,nbf)
 character(len=240), intent(in) :: fchname
 character(len=240) :: file47

 call call_gaussian_gen47_from_fch(fchname, file47)
 call read_ao_ovlp_from_47(file47, nbf, S)
 call delete_file(file47)
end subroutine get_ao_ovlp_using_fch

! call Gaussian to compute AO-basis dipole integrals using the given .fch file
subroutine get_ao_dipole_using_fch(fchname, nbf, D)
 implicit none
 integer, intent(in) :: nbf
 real(kind=8), intent(out) :: D(nbf,nbf,3)
 character(len=240), intent(in) :: fchname
 character(len=240) :: file47

 call call_gaussian_gen47_from_fch(fchname, file47)
 call read_ao_dipole_from_47(file47, nbf, D)
 call delete_file(file47)
end subroutine get_ao_dipole_using_fch

! call Gaussian program to generate .47 file from a given .fch(k) file
subroutine call_gaussian_gen47_from_fch(fchname, file47)
 use util_wrapper, only: unfchk
 implicit none
 integer :: i, fid, system
 character(len=10) :: str
 character(len=24) :: mem
 character(len=240), intent(in) :: fchname
 character(len=240), intent(out) :: file47
 character(len=240) :: gau_path, proname, chkname, gjfname, logname

 str = REPEAT(' ',10)
 call get_a_random_int(i)
 write(str,'(I0)') i

 i = index(fchname, '.fch', back=.true.)
 if(i == 0) then
  write(6,'(A)') "ERROR in subroutine call_gaussian_gen47_from_fch: no '.fch'&
                & suffix found in"
  write(6,'(A)') 'filename '//TRIM(fchname)
  stop
 end if

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
 if(LEN_TRIM(mem) == 0) mem = '1GB'

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
 i = system(TRIM(gau_path)//' '//TRIM(gjfname))
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
 forall(i=1:3,j=1:nbf,k=1:nbf,j<k) D(k,j,i) = D(j,k,i)
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

! get a random integer
subroutine get_a_random_int(i)
 implicit none
 integer :: n, clock
 integer, intent(out) :: i
 integer, allocatable :: seed(:)
 real(kind=4) :: r

 call random_seed(size=n)
 allocate(seed(n))
 call system_clock(count=clock)
 seed = clock
 call random_seed(put=seed)
 call random_number(r)
 deallocate(seed)

 i = CEILING(r*1e6)
 return
end subroutine get_a_random_int

! read the path of the Gaussian binary executable file 
subroutine get_gau_path(gau_path)
 implicit none
 integer :: i
 character(len=240) :: g(3)
 character(len=240), intent(out) :: gau_path
 logical :: alive

 gau_path = ' '
 call getenv('GAUSS_EXEDIR', gau_path)

#ifdef _WIN32
 i = index(gau_path, '\', back=.true.)
 if(i == 0) then
  write(6,'(A)') "ERROR in subroutine get_gau_path: no '\' symbol found in&
                  & gau_path="//TRIM(gau_path)
  stop
 end if
 gau_path = """"//TRIM(gau_path)//'\g'//gau_path(i+2:i+3)//".exe"""

#else
 i = index(gau_path, ':', back=.true.)
 if(i == 0) then
  write(6,'(/,A)') "ERROR in subroutine get_gau_path: no ':' symbol found&
                  & in gau_path="//TRIM(gau_path)
  write(6,'(/,A)') 'This error often occurs when your machine has no (or has&
                  & incorrect) Gaussian'
  write(6,'(A)') 'environment variables. Here I offer a correct example:'
  write(6,'(A)') REPEAT('-',45)
  write(6,'(A)') ' export g16root=/opt'
  write(6,'(A)') ' source $g16root/g16/bsd/g16.profile'
  write(6,'(A)') ' export GAUSS_SCRDIR=/scratch/$USER/gaussian'
  write(6,'(A)') REPEAT('-',45)
  write(6,'(A)') 'Please check your Gaussian environment variables according&
                   & to the example shown above.'
  write(6,'(A,/)') "Also note: DO NOT write 'export GAUSS_EXEDIR', it is useless."
  stop
 end if

 gau_path = gau_path(i+1:)
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

! check whether two double precision values equal to each other
subroutine check_two_real8_eq(r1, r2, diff)
 implicit none
 real(kind=8), intent(in) :: r1, r2, diff

 if(diff < 0d0) then
  write(6,'(A)') 'ERROR in subroutine check_two_real8_eq: input diff<0.'
  write(6,'(A,F14.8,A)') 'diff=',diff,', but diff>0 is required.'
  stop
 end if

 if(DABS(r1-r2) > diff) then
  write(6,'(A)') 'ERROR in subroutine check_two_real8_eq: two double precisi&
                    &on values differ larger than the'
  write(6,'(3(A,F14.8))') 'tolerance diff. r1=',r1,', r2=',r2,', diff=',diff
  stop
 end if
end subroutine check_two_real8_eq

! submit a GAMESS job
subroutine submit_gms_job(gms_path, gms_scr_path, inpname, nproc)
 implicit none
 integer :: i, system
 integer, intent(in) :: nproc
 character(len=240) :: datname, gmsname, hs1, hs2, trj
 character(len=500) :: longbuf
 character(len=240), intent(in) :: gms_path, gms_scr_path, inpname

 i = index(inpname, '.inp', back=.true.)
 if(i == 0) then
  write(6,'(/,A)') "ERROR in subroutine submit_gms_job: '.inp' suffix not &
                      &found in filename "//TRIM(inpname)
  stop
 end if

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

 i = system(TRIM(gms_path)//' '//TRIM(longbuf))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine submit_gms_job: GAMESS job failed.'
  write(6,'(A)') 'You can open file '//TRIM(gmsname)//' and check why.'
  stop
 end if

 ! move the .dat file into current directory
 i = system('mv '//TRIM(gms_scr_path)//'/'//TRIM(datname)//' .')
 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine submit_gms_job: fail to move file. Possi&
                 &bly wrong gms_scr_path.'
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

 i = index(gjfname,'.gjf',back=.true.)
 outname = gjfname(1:i-1)//'.out'
 buf = 'automr '//TRIM(gjfname)//' >'//TRIM(outname)//" 2>&1"
 write(6,'(A)') '$'//TRIM(buf)

 i = system(TRIM(buf))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine submit_automr_job: automr job failed.'
  write(6,'(A)') 'You can open file '//TRIM(outname)//' and check why.'
  stop
 end if
end subroutine submit_automr_job

subroutine submit_gau_job(gau_path, gjfname)
 implicit none
 integer :: i, system
 character(len=240) :: logname
 character(len=240), intent(in) :: gau_path, gjfname

 i = index(gjfname, '.gjf', back=.true.)
#ifdef _WIN32
 logname = gjfname(1:i-1)//'.out'
#else
 logname = gjfname(1:i-1)//'.log'
#endif

 write(6,'(A)') '$'//TRIM(gau_path)//' '//TRIM(gjfname)
 i = system(TRIM(gau_path)//' '//TRIM(gjfname))

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subrouitine submit_gau_job: Gaussian job failed.'
  write(6,'(A)') 'Please open file '//TRIM(logname)//' and check.'
  stop
 end if
end subroutine submit_gau_job

subroutine submit_orca_job(orca_path, inpname)
 implicit none
 integer :: i, system
 character(len=240) :: outname
 character(len=480) :: buf
 character(len=240), intent(in) :: orca_path, inpname

 i = index(inpname, '.inp', back=.true.)
 outname = inpname(1:i-1)//'.out'

 write(buf,'(A)') TRIM(inpname)//' >'//TRIM(outname)//" 2>&1"
 write(6,'(A)') '$orca '//TRIM(buf)
 i = system(TRIM(orca_path)//' '//TRIM(buf))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subrouitine submit_orca_job: ORCA job failed.'
  write(6,'(A)') 'Please open file '//TRIM(outname)//' and check.'
  stop
 end if
end subroutine submit_orca_job

subroutine submit_molcas_job(inpname)
 implicit none
 integer :: i, system
 character(len=240) :: outname
 character(len=480) :: buf
 character(len=240), intent(in) :: inpname

 i = index(inpname, '.inp', back=.true.)
 outname = inpname(1:i-1)//'.out'

 write(buf,'(A)') 'pymolcas '//TRIM(inpname)//' >'//TRIM(outname)//" 2>&1"
 write(6,'(A)') '$'//TRIM(buf)

 i = system(TRIM(buf))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subrouitine submit_molcas_job: OpenMolcas job failed.'
  write(6,'(A)') 'Please open file '//TRIM(outname)//' and check.'
  stop
 end if
end subroutine submit_molcas_job

subroutine submit_psi4_job(psi4_path, inpname, nproc)
 implicit none
 integer :: i, system
 integer, intent(in) :: nproc
 character(len=240) :: outname
 character(len=480) :: buf
 character(len=240), intent(in) :: psi4_path, inpname

 i = index(inpname, '.inp', back=.true.)
 outname = inpname(1:i-1)//'.out'

 write(buf,'(A,I0)') TRIM(inpname)//' '//TRIM(outname)//' -n ', nproc
 write(6,'(A)') '$psi4 '//TRIM(buf)

 i = system(TRIM(psi4_path)//' '//TRIM(buf))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subrouitine submit_psi4_job: PSI4 job failed.'
  write(6,'(A)') 'Please open file '//TRIM(outname)//' and check.'
  stop
 end if
end subroutine submit_psi4_job

