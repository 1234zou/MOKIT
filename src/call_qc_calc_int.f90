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
 return
end subroutine gen_no_using_density_in_fch

! call Gaussian to compute AO-basis overlap according to a given .fch file
subroutine get_ao_ovlp_using_fch(fchname, nbf, S)
 use print_id, only: iout
 use util_wrapper, only: unfchk
 implicit none
 integer :: i, fid, system
 integer, intent(in) :: nbf
 real(kind=8), intent(out) :: S(nbf,nbf)
 character(len=10) :: str, nproc
 character(len=24) :: mem
 character(len=240), intent(in) :: fchname
 character(len=240) :: gau_path, proname, chkname, gjfname, logname, file47

 str = ' '
 call get_a_random_int(i)
 write(str,'(I0)') i

 i = index(fchname, '.fch', back=.true.)
 if(i == 0) then
  write(iout,'(A)') "ERROR in subroutine get_ao_ovlp_using_fch: no '.fch'&
                   & suffix found in filename "//TRIM(fchname)
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

 call getenv('OMP_NUM_THREADS', nproc)
 if(LEN_TRIM(nproc) == 0) call getenv('MKL_NUM_THREADS', nproc)
 if(LEN_TRIM(nproc) == 0) nproc = '1'

 open(newunit=fid,file=TRIM(gjfname),status='replace')
 write(fid,'(A)') '%mem='//TRIM(mem)
 write(fid,'(A)') '%nprocshared='//TRIM(nproc)
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
  write(iout,'(/,A)') 'ERROR in subroutine get_ao_ovlp_using_fch: Gaussian job failed.'
  write(iout,'(A)') 'You can open file '//TRIM(logname)//' and check why.'
  stop
 end if
 call delete_files(3, [chkname, gjfname, logname])
 call read_ao_ovlp_from_47(file47, nbf, S)
 call delete_file(file47)
 return
end subroutine get_ao_ovlp_using_fch

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
 use print_id, only: iout
 implicit none
 integer :: i
 character(len=240), intent(out) :: gau_path

 gau_path = ' '
 call getenv('GAUSS_EXEDIR', gau_path)

#ifdef _WIN32
 i = index(gau_path, '\', back=.true.)
 if(i == 0) then
  write(iout,'(A)') "ERROR in subroutine get_gau_path: no '\' symbol found in&
                  & gau_path="//TRIM(gau_path)
  stop
 end if
 gau_path = """"//TRIM(gau_path)//'\g'//gau_path(i+2:i+3)//".exe"""

#else
 i = index(gau_path, ':', back=.true.)
 if(i == 0) then
  write(iout,'(/,A)') "ERROR in subroutine get_gau_path: no ':' symbol found&
                     & in gau_path="//TRIM(gau_path)
  write(iout,'(/,A)') 'This error often occurs when your machine has no (or has&
                     & incorrect) Gaussian'
  write(iout,'(A)') 'environment variables. Here I offer a correct example:'
  write(iout,'(A)') REPEAT('-',45)
  write(iout,'(A)') ' export g16root=/opt'
  write(iout,'(A)') ' source $g16root/g16/bsd/g16.profile'
  write(iout,'(A)') ' export GAUSS_SCRDIR=/scratch/$USER/gaussian'
  write(iout,'(A)') REPEAT('-',45)
  write(iout,'(A)') 'Please check your Gaussian environment variables according&
                   & to the example shown above.'
  write(iout,'(A)') "Also note: DO NOT write 'export GAUSS_EXEDIR', it is useless."
  stop
 end if

 gau_path = gau_path(i+1:)
 i = index(gau_path, '/', back=.true.)
 gau_path = TRIM(gau_path)//'/'//TRIM(gau_path(i+1:))
#endif

 return
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
 return
end subroutine get_ne_from_PS

! check whether two double precision values equal to each other
subroutine check_two_real8_eq(r1, r2, diff)
 use print_id, only: iout
 implicit none
 real(kind=8), intent(in) :: r1, r2, diff

 if(diff < 0d0) then
  write(iout,'(A)') 'ERROR in subroutine check_two_real8_eq: input diff<0.'
  write(iout,'(A,F14.8,A)') 'diff=',diff,', but diff>0 is required.'
  stop
 end if

 if(DABS(r1-r2) > diff) then
  write(iout,'(A)') 'ERROR in subroutine check_two_real8_eq: two double precisi&
                    &on values differ larger than the'
  write(iout,'(3(A,F14.8))') 'tolerance diff. r1=',r1,', r2=',r2,', diff=',diff
  stop
 end if

 return
end subroutine check_two_real8_eq

