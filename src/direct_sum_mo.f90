! Perform the direct sum of fragment MOs in a series of .fch files (fchname0),
!  and write merged MOs into fchname.
! Note:
!  1) the wave function type (wfn_type) of each .fch file must be provided since
!     this subroutine allows extending R(O)HF->UHF automatically
!  2) the positive/negative spin information must be provided since this subroutine
!     allows negative spins like anti-ferromagnetic coupling
subroutine direct_sum_frag_mo_in_fch(n, fchname0, wfn_type0, pos, fchname, wfn_type)
 implicit none
 integer :: i, k1, k2, k3, nbf0, nif0, nbf, nif
 integer :: na0, nb0, na, nb   ! No. of alpha/beta electrons
 integer, intent(in) :: n
 integer, intent(in) :: wfn_type0(n), wfn_type ! 1/2/3 for RHF/ROHF/UHF
 character(len=240), intent(in) :: fchname0(n), fchname
 real(kind=8), allocatable :: mo_a0(:,:), mo_b0(:,:) ! fragment MOs
 real(kind=8), allocatable :: mo_a(:,:), mo_b(:,:)   ! supermolecule MOs
 real(kind=8), allocatable :: ovlp(:,:) ! AO overlap matrix
 logical, intent(in) :: pos(n)

 if(ANY(wfn_type0==0) .or. wfn_type==0) then
  write(6,'(/,A)') 'ERROR in subroutine direct_sum_frag_mo_in_fch: there exists&
                   & at least one'
  write(6,'(A)') 'file whose wfn_type is 0 (i.e. undetermined).'

  do i = 1, n, 1
   write(6,'(2(A,I0))') 'i=', i, ','//TRIM(fchname0(i))//',wfn_type0(i)=', &
                        wfn_type0(i)
  end do ! for i
  write(6,'(A,I0)') 'wfn_type=', wfn_type
  stop
 end if

 call read_na_and_nb_from_fch(fchname, na, nb)
 call read_nbf_and_nif_from_fch(fchname, nbf, nif)

 allocate(mo_a(nbf,nif), source=0d0)
 if(wfn_type == 3) allocate(mo_b(nbf,nif), source=0d0)
 k1 = 0; k2 = 0; k3 = 0

 do i = 1, n, 1
  call read_nbf_and_nif_from_fch(fchname0(i), nbf0, nif0)
  allocate(mo_a0(nbf0,nif0))
  if(pos(i)) then ! positive spin
   call read_na_and_nb_from_fch(fchname0(i), na0, nb0)
   call read_mo_from_fch(fchname0(i), nbf0, nif0, 'a', mo_a0)
  else            ! negative spin
   call read_na_and_nb_from_fch(fchname0(i), nb0, na0)
   if(wfn_type0(i) == 3) then
    call read_mo_from_fch(fchname0(i), nbf0, nif0, 'b', mo_a0)
   else ! ROHF
    call read_mo_from_fch(fchname0(i), nbf0, nif0, 'a', mo_a0)
   end if
  end if

  if(k1+nbf0 > nbf) then
   write(6,'(/,A)') 'ERROR in subroutine direct_sum_frag_mo_in_fch: the 1st dim&
                    &ension of'
   write(6,'(A)') 'array mo_a out of range!'
   write(6,'(5(A,I0))') 'k1=',k1,',nbf0=',nbf0,',nbf=',nbf,',i=',i,',n=',n
   stop
  end if
  mo_a(k1+1:k1+nbf0, k2+1:k2+na0) = mo_a0(:,1:na0)

  if(wfn_type == 3) then ! UHF
   if(wfn_type0(i) == 3) then ! UHF fragment
    allocate(mo_b0(nbf0,nif0))
    if(pos(i)) then ! positive spin
     call read_mo_from_fch(fchname0(i), nbf0, nif0, 'b', mo_b0)
    else            ! negative spin
     call read_mo_from_fch(fchname0(i), nbf0, nif0, 'a', mo_b0)
    end if
    mo_b(k1+1:k1+nbf0, k3+1:k3+nb0) = mo_b0(:,1:nb0)
    deallocate(mo_b0)
   else ! RHF/ROHF fragment
    mo_b(k1+1:k1+nbf0, k3+1:k3+nb0) = mo_a0(:,1:nb0)
   end if
   k3 = k3 + nb0
  end if

  deallocate(mo_a0)
  k1 = k1 + nbf0
  k2 = k2 + na0
 end do ! for i

 ! In previous version, virtual MOs of the supermolecule were constructed from
 ! virtual MOs of fragments. But it did not work for linear dependence case,
 ! so I change it to construction of virtual MOs using PAO

 ! TODO: construct virtual MOs from fragments, do canonical orthonormalization
 !  if linear dependence occurs, do symmetric orthonormalization if there is no
 !  linear dependence

 allocate(ovlp(nbf,nbf), mo_a0(nbf,nif))
 call get_ao_ovlp_using_fch(fchname, nbf, ovlp)
 ! orthonormalize occupied MOs
 call orthonormalize_orb(.true., nbf, k2, ovlp, mo_a(:,1:k2), mo_a0(:,1:k2))
 call construct_vir(nbf, nif, k2+1, mo_a0, ovlp, mo_a)
 deallocate(mo_a0)
 call write_mo_into_fch(fchname, nbf, nif, 'a', mo_a)
 deallocate(mo_a)

 if(wfn_type == 3) then
  allocate(mo_b0(nbf,nif))
  call orthonormalize_orb(.true., nbf, k3, ovlp, mo_b(:,1:k3), mo_b0(:,1:k3))
  call construct_vir(nbf, nif, k3+1, mo_b0, ovlp, mo_b)
  deallocate(mo_b0)
  call write_mo_into_fch(fchname, nbf, nif, 'b', mo_b)
  deallocate(mo_b)
 end if

 deallocate(ovlp)
end subroutine direct_sum_frag_mo_in_fch

subroutine direct_sum_frag_fock_in_fch(n, fchname0, wfn_type0, pos, fchname, wfn_type)
 implicit none
 integer :: i, j, k1, k2, k3, nbf0, nif0, nbf, nif
 integer :: na0, nb0, na, nb   ! No. of alpha/beta electrons
 integer, intent(in) :: n
 integer, intent(in) :: wfn_type0(n), wfn_type ! 1/2/3 for RHF/ROHF/UHF
 character(len=240) :: logname
 character(len=240), intent(in) :: fchname0(n), fchname
 real(kind=8) :: ne
 real(kind=8), allocatable :: mo_a0(:,:), mo_b0(:,:) ! fragment MOs
 real(kind=8), allocatable :: mo_a(:,:), mo_b(:,:)   ! supermolecule MOs
 real(kind=8), allocatable :: e_a0(:), e_b0(:)       ! fragment MO energies
 real(kind=8), allocatable :: fock_a0(:,:), fock_b0(:,:) ! fragment Fock
 real(kind=8), allocatable :: fock_a(:,:), fock_b(:,:)   ! supermolecule Fock
 real(kind=8), allocatable :: ovlp(:,:) ! AO overlap matrix
 real(kind=8), allocatable :: U(:,:), X(:,:)
 logical, intent(in) :: pos(n)

 if(ANY(wfn_type0==0) .or. wfn_type==0) then
  write(6,'(/,A)') 'ERROR in subroutine direct_sum_frag_fock_in_fch: there exis&
                   &ts at least one'
  write(6,'(A)') 'file whose wfn_type is 0 (i.e. undetermined).'

  do i = 1, n, 1
   write(6,'(2(A,I0))') 'i=', i, ','//TRIM(fchname0(i))//',wfn_type0(i)=', &
                        wfn_type0(i)
  end do ! for i
  write(6,'(A,I0)') 'wfn_type=', wfn_type
  stop
 end if

 call read_na_and_nb_from_fch(fchname, na, nb)
 call read_nbf_and_nif_from_fch(fchname, nbf, nif)

 allocate(fock_a(nbf,nbf), source=0d0)
 if(wfn_type == 3) allocate(fock_b(nbf,nbf), source=0d0)

 ! initialize F = Hcore

 logname = 'DMC_trimer_rhf_hcore.log'
! call read_hcore_from_gaulog(logname, nbf, fock_a)

 k1 = 0; k2 = 0; k3 = 0

 do i = 1, n, 1
  call read_nbf_and_nif_from_fch(fchname0(i), nbf0, nif0)
  allocate(mo_a0(nbf0,nif0), e_a0(nif0), fock_a0(nbf0,nbf0))
  if(pos(i)) then ! positive spin
   call read_na_and_nb_from_fch(fchname0(i), na0, nb0)
   call read_mo_from_fch(fchname0(i), nbf0, nif0, 'a', mo_a0)
   call read_eigenvalues_from_fch(fchname0(i), nif0, 'a', e_a0)
  else            ! negative spin
   call read_na_and_nb_from_fch(fchname0(i), nb0, na0)
   if(wfn_type0(i) == 3) then
    call read_mo_from_fch(fchname0(i), nbf0, nif0, 'b', mo_a0)
    call read_eigenvalues_from_fch(fchname0(i), nif0, 'b', e_a0)
   else ! ROHF
    call read_mo_from_fch(fchname0(i), nbf0, nif0, 'a', mo_a0)
    call read_eigenvalues_from_fch(fchname0(i), nif0, 'a', e_a0)
   end if
  end if
  call solve_fock_from_ctfc(nbf0, nif0, mo_a0, e_a0, fock_a0)
  deallocate(mo_a0, e_a0)

  if(k1+nbf0 > nbf) then
   write(6,'(/,A)') 'ERROR in subroutine direct_sum_frag_fock_in_fch: the 1st d&
                    &imension of'
   write(6,'(A)') 'array mo_a out of range!'
   write(6,'(5(A,I0))') 'k1=',k1,',nbf0=',nbf0,',nbf=',nbf,',i=',i,',n=',n
   stop
  end if
  fock_a(k1+1:k1+nbf0, k1+1:k1+nbf0) = fock_a0

  if(wfn_type == 3) then ! UHF
   if(wfn_type0(i) == 3) then ! UHF fragment
    allocate(mo_b0(nbf0,nif0), e_b0(nif0), fock_b0(nbf0,nbf0))
    if(pos(i)) then ! positive spin
     call read_mo_from_fch(fchname0(i), nbf0, nif0, 'b', mo_b0)
     call read_eigenvalues_from_fch(fchname0(i), nif0, 'b', e_b0)
    else            ! negative spin
     call read_mo_from_fch(fchname0(i), nbf0, nif0, 'a', mo_b0)
     call read_eigenvalues_from_fch(fchname0(i), nif0, 'a', e_b0)
    end if
    call solve_fock_from_ctfc(nbf0, nif0, mo_b0, e_b0, fock_b0)
    deallocate(mo_b0, e_b0)
    fock_b(k1+1:k1+nbf0, k1+1:k1+nbf0) = fock_b0
    deallocate(fock_b0)
   else                       ! R(O)HF fragment
    fock_b(k1+1:k1+nbf0, k1+1:k1+nbf0) = fock_a0
   end if
   k3 = k3 + nb0
  end if

  deallocate(fock_a0)
  k1 = k1 + nbf0
  k2 = k2 + na0
 end do ! for i

 ! In previous version, virtual MOs of the supermolecule were constructed from
 ! virtual MOs of fragments. But it did not work for linear dependence case,
 ! so I change it to construction of virtual MOs using PAO

 ! TODO: construct virtual MOs from fragments, do canonical orthonormalization
 !  if linear dependence occurs, do symmetric orthonormalization if there is no
 !  linear dependence

 allocate(ovlp(nbf,nbf))
 call get_ao_ovlp_using_fch(fchname, nbf, ovlp)

 ! do 1 cycle FC=SCE using ovlp and the current fock_a
 allocate(U(nbf,nbf), source=ovlp)
 allocate(e_a0(nbf))
 call diag_get_e_and_vec(nbf, U, e_a0)
 forall(i = 1:nbf) e_a0(i) = 1d0/DSQRT(e_a0(i))
 allocate(X(nbf,nbf))   ! X = U((s')^(-1/2))
 forall(i=1:nbf, j=1:nbf) X(i,j) = U(i,j)*e_a0(j)
 call calc_CTSC(nbf, nbf, X, fock_a, U)
 deallocate(fock_a)
 call diag_get_e_and_vec(nbf, U, e_a0)
 deallocate(e_a0)
 allocate(mo_a(nbf,nif), source=0d0)
 mo_a = MATMUL(X, U)
 deallocate(X, U)

 call write_mo_into_fch(fchname, nbf, nif, 'a', mo_a)
 ovlp = MATMUL(2d0*MATMUL(mo_a(:,1:na), TRANSPOSE(mo_a(:,1:na))), ovlp)
 ne = 0d0
 do i = 1, nbf, 1
  ne = ne + ovlp(i,i)
 end do ! for i
 write(6,'(A,F18.8)') 'trace(PS) = ', ne
 deallocate(mo_a)

! if(wfn_type == 3) then
!  allocate(mo_b0(nbf,nif))
!  call orthonormalize_orb(.true., nbf, k3, ovlp, mo_b(:,1:k3), mo_b0(:,1:k3))
!  call construct_vir(nbf, nif, k3+1, mo_b0, ovlp, mo_b)
!  deallocate(mo_b0)
!  call write_mo_into_fch(fchname, nbf, nif, 'b', mo_b)
!  deallocate(mo_b)
! end if

 deallocate(ovlp)
 if(wfn_type == 3) deallocate(fock_b)
end subroutine direct_sum_frag_fock_in_fch

! Direct sum the total SCF densities of each fragment and write the sum into
! Total SCF Density section of fname. This subroutine requires each .fch file
! has fewer basis functions, and thus it is "direct sum".
! For a simple sum, please use subroutine sum_frag_dm_in_fch().
subroutine direct_sum_frag_dm_in_fch(n, fname0, fname)
 use util_wrapper, only: formchk
 implicit none
 integer :: i, j, nbf0, nif0, nbf, nif
 integer, intent(in) :: n
 real(kind=8), allocatable :: dm0(:,:), dm(:,:), ovlp(:,:), mo(:,:), occ(:)
 character(len=240) :: chkname
 character(len=240), intent(in) :: fname0(n), fname
 character(len=240), allocatable :: fchname(:)
 logical :: alive

 allocate(fchname(n+1))
 fchname(1:n) = fname0
 fchname(n+1) = fname
 ! if fname0(n) and fname are .gjf filenames, convert them into .fch
 call fnames_gjf2fch(n+1, fchname)

 call find_specified_suffix(fchname(n+1), '.fch', i)
 chkname = fchname(n+1)(1:i-1)//'.chk'
 inquire(file=TRIM(fchname(n+1)), exist=alive)
 if(.not. alive) call formchk(chkname, fchname(n+1))

 call read_nbf_and_nif_from_fch(fchname(n+1), nbf0, nif0)
 allocate(ovlp(nbf0,nbf0))
 call get_ao_ovlp_using_fch(fchname(n+1), nbf0, ovlp)
 allocate(dm0(nbf0,nbf0), source=0d0)

 j = 0
 do i = 1, n, 1
  call read_nbf_and_nif_from_fch(fchname(i), nbf, nif)
  write(6,'(A,I0)') 'nbf=', nbf
  allocate(dm(nbf,nbf))
  call read_dm_from_fch(fchname(i), 1, nbf, dm)
  dm0(j+1:j+nbf,j+1:j+nbf) = dm
  deallocate(dm)
  j = j + nbf
 end do ! for i

 if(j /= nbf0) then
  write(6,'(/,A)') 'ERROR in subroutine direct_sum_frag_dm_in_fch: the sum of n&
                   &umber of basis'
  write(6,'(A)') 'func. from each fragment is not equal to the number of basis &
                 &func. of the complex.'
  stop
 end if

 allocate(occ(nif0), mo(nbf0,nif0))
 call gen_no_from_density_and_ao_ovlp(nbf0, nif0, dm0, ovlp, occ, mo)
 call write_eigenvalues_to_fch(fchname(n+1), nif0, 'a', occ, .true.)
 deallocate(dm0, ovlp, occ)
 ! There is no need to write the occupation numbers and total density into .fch,
 ! since they do not affect the SCF computation.
 call write_mo_into_fch(fchname(n+1), nbf0, nif0, 'a', mo)
 deallocate(fchname, mo)
end subroutine direct_sum_frag_dm_in_fch

! Sum the total SCF densities of each fragment and write the sum into Total SCF
! Density section of fname. This subroutine requires each .fch file has the same
! number of basis functions, and thus it is not "direct sum", but simply sum.
! For "direct sum", please use subroutine direct_sum_frag_dm_in_fch().
subroutine sum_frag_dm_in_fch(n, fname0, pos, fname)
 use util_wrapper, only: formchk
 implicit none
 integer :: i, nif, nbf, nbf1
 integer, intent(in) :: n
 character(len=240) :: chkname
 character(len=240), intent(in) :: fname0(n), fname
 character(len=240), allocatable :: fchname(:)
 real(kind=8), allocatable :: dm(:,:), dm1(:,:)
 logical :: alive
 logical, allocatable :: has_spin_density(:)
 logical, intent(in) :: pos(n) ! negative spin means to switch alpha/beta

 allocate(fchname(n+1))
 fchname(1:n) = fname0
 fchname(n+1) = fname
 ! if fname0(n) and fname are .gjf files, convert filenames into .fch
 call fnames_gjf2fch(n+1, fchname)

 call find_specified_suffix(fname, '.fch', i)
 chkname = fname(1:i-1)//'.chk'
 inquire(file=TRIM(fchname(n+1)), exist=alive)
 if(.not. alive) call formchk(chkname, fchname(n+1))

 call read_nbf_and_nif_from_fch(fchname(n+1), nbf, nif)
 allocate(dm(nbf,nbf), dm1(nbf,nbf))
 dm = 0d0

 do i = 1, n, 1
  call read_nbf_and_nif_from_fch(fchname(i), nbf1, nif)
  if(nbf1 /= nbf) then
   write(6,'(A)') 'ERROR in subroutine sum_frag_dm_in_fch: nbf1 /= nbf.'
   write(6,'(A,I0,A)') 'Inconsistent nbf between fragment ',i,' and the&
                       & total system.'
   write(6,'(2(A,I0))') 'nbf1=', nbf1, ', nbf=', nbf
   stop
  end if
  call read_dm_from_fch(fchname(i), 1, nbf, dm1)
  dm = dm + dm1
 end do ! for i
 call write_dm_into_fch(fchname(n+1), .true., nbf, dm)

 allocate(has_spin_density(n))
 has_spin_density = .false.
 dm = 0d0

 do i = 1, n, 1
  call detect_spin_scf_density_in_fch(fchname(i), has_spin_density(i))
  if(.not. has_spin_density(i)) cycle
  call read_dm_from_fch(fchname(i), 2, nbf, dm1)
  if(pos(i)) then
   dm = dm + dm1
  else
   dm = dm - dm1 ! interchange alpha/beta spin, which is equal to -dm1
  end if
 end do ! for i

 call detect_spin_scf_density_in_fch(fchname(n+1), alive)

 if(ANY(has_spin_density .eqv. .true.)) then
  if(alive) then
   call write_dm_into_fch(fchname(n+1), .false., nbf, dm)
  else ! no Spin SCF Density in the total system .fch file
   write(6,'(/,A)') 'ERROR in subroutine sum_frag_dm_in_fch: some fragment has &
                    &UHF-type wave fun-'
   write(6,'(A)') 'ction, but the total system has RHF-type wave function. Inco&
                  &nsistency detected.'
   stop
  end if
 else ! all fragments has RHF-type wave function
  if(alive) call write_dm_into_fch(fchname(n+1), .false., nbf, dm)
 end if

 deallocate(fchname, dm, dm1, has_spin_density)
end subroutine sum_frag_dm_in_fch

subroutine read_hcore_from_gaulog(logname, nbf, hcore)
 implicit none
 integer :: i, j, k, m, nbatch, fid
 integer, intent(in) :: nbf
 real(kind=8), intent(out) :: hcore(nbf,nbf)
 character(len=240) :: buf
 character(len=240), intent(in) :: logname

 hcore = 0d0
 nbatch = nbf/5
 nbatch = nbatch + 1
 open(newunit=fid,file=TRIM(logname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(9:16) == 'Core Ham') exit
 end do ! for while

 do i = 1, nbatch, 1
  read(fid,'(A)') buf
  k = 5*(i-1)
  do j = k+1, nbf, 1
   read(fid,*) m, hcore(j:j+min(j-k,5)-1,j)
  end do ! for j
 end do ! for i

 close(fid)
 forall(i=1:nbf, j=1:nbf, j>i) hcore(i,j) = hcore(j,i)
end subroutine read_hcore_from_gaulog

! convert multiple .gjf filenames into corresponding .fch filenames
subroutine fnames_gjf2fch(n, fname)
 implicit none
 integer :: i
 integer, intent(in) :: n
 character(len=240), intent(inout) :: fname(n)

 do i = 1, n, 1
  call fname_gjf2fch(fname(i))
 end do ! for i
end subroutine fnames_gjf2fch

! convert a .gjf filename into a .fch filename
subroutine fname_gjf2fch(fname)
 implicit none
 integer :: i
 character(len=240), intent(inout) :: fname

 i = LEN_TRIM(fname)

 select case(fname(i-3:i))
 case('.fch') ! do nothing
 case('.gjf')
  fname = fname(1:i-4)//'.fch'
 case default
  write(6,'(/,A)') 'ERROR in subroutine fname_gjf2fch: suffix cannot be recogni&
                   &zed.'
  write(6,'(A)') 'fname='//TRIM(fname)
  stop
 end select
end subroutine fname_gjf2fch

