! written by jxzou at 20180401
! updated by jxzou at 20200517: consider the parallelism of all MOs in pairs,
!  not just bonding v.s. anti-bonding; simplify code.
! updated by jxzou at 20200518: add logical variable 'gau'

! WHY we need this subroutine: in rare cases(very small CI coefficients of the
!  2nd NOs, the 1st NOs are just core orbitals), GAMESS output(.dat file) MOs
!  are in wrong order, pairs containing 1st NOs as core orbitals will go across
!  other pairs. So I hope this subroutine can correct the order.

program main
 implicit none
 integer :: i, idx1, idx2
 character(len=240) :: buf = ' ', fname = ' '
 logical :: gau
 ! gau = .True., meaning MOs in .fch(k) file are in Gaussian order
 !  (bonding1, bonding2, ..., open, ..., anti-bonding2, anti-bonding1)
 ! gau = .False., meaning MOs in .fch(k) file are in GAMESS order
 !  (open, bonding1, anti-bonding1, bonding2, anti-bonding2, ...)

 idx1 = 0; idx2 = 0

 i = iargc()
 if(i /= 4) then
  write(6,'(/,A)') ' ERROR in gvb_correct_pairs: wrong command line arguments!'
  write(6,'(/,A)') ' Format: gvb_correct_pairs a.fch idx1 idx2 -gms/-gau'
  write(6,'(/,A)') ' Example1: gvb_correct_pairs a.fch 3 12 -gms'
  write(6,'(/,A,/)') ' Example2: gvb_correct_pairs b.fch 3 12 -gau'
  stop
 end if

 call getarg(1,fname)
 call getarg(2,buf)
 read(buf,*) idx1
 call getarg(3,buf)
 read(buf,*) idx2

 call getarg(4,buf)
 buf = ADJUSTL(buf)
 select case(TRIM(buf))
 case('-gau')
  gau = .true.
 case('-gms')
  gau = .false.
 case default
  write(6,'(A)') 'ERROR in gvb_correct_pairs: the 4-th argument is wrong.'
  write(6,'(A)') 'buf = '//TRIM(buf)
  stop
 end select

 if(idx1 >= idx2) then
  write(6,'(A)') 'ERROR in subroutine gvb_correct_pairs: idx1 >= idx2.'
  write(6,'(A,2I5)') 'idx1/idx2 = ', idx1, idx2
  stop
 end if

 call gvb_correct_pairs(fname, idx1, idx2, gau)
 stop
end program main

! gvb_correct_pairs: pair the disordered pair orbitals in a .fch(k) file
subroutine gvb_correct_pairs(fchname, idx1, idx2, gau)
 implicit none
 integer :: i, j, m, n, tmp_idx(1), iter
 integer :: nif, nbf, nmo, npair, nopen
 ! nif: the number of independent functions, i.e., the number of MOs
 ! nbf: the number of basis functions
 ! nmo: the number of target orbitals
 integer, intent(in) :: idx1, idx2
 ! idx1: the beginning index of target MOs (Fortran convention, 1,...)
 ! idx2: the final index of target MOs (Fortran convention, 1,...)
 integer, parameter :: iter_max = 1000
 integer, allocatable :: target_pair_idx(:), ideal_idx(:)
 real(kind=8), allocatable :: coeff(:,:), coeff1(:,:)
 real(kind=8), allocatable :: angle(:,:)
 real(kind=8), allocatable :: tmp_coeff(:), norm(:)
 real(kind=8) :: tmp_norm
 character(len=240), intent(in) :: fchname
 logical :: same_idx
 logical, intent(in) :: gau
 logical, allocatable :: paired(:)

 write(6,'(/,A)') 'Enter subroutine gvb_correct_pairs (check if paired MOs are correct)...'

 if(idx2 - idx1 < 2) return
 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(coeff(nbf,nif), source=0.0d0)
 call read_mo_from_fch(fchname, nbf, nif, 'a', coeff)

 ! Task: sort the angles between any two orbitals by ascent order and pair them
 ! 1st step: turn all target MO coefficients into positive and calculate their norms
 nmo = idx2 - idx1 + 1
 allocate(coeff1(nbf,nmo), source=0.0d0)
 allocate(norm(nmo), source=0.0d0)
 forall(i = 1:nmo, j = 1:nbf) coeff1(j,i) = DABS(coeff(j,idx1-1+i))
 forall(i = 1:nmo) norm(i) = DOT_PRODUCT(coeff1(:,i), coeff1(:,i))
 forall(i = 1:nmo) norm(i) = DSQRT(norm(i))

 ! 2nd step: calculate the angles between any two active orbitals
 ! the angle here and below are in fact in its COSINE value
 allocate(angle(nmo, nmo), source=0.0d0)
 allocate(tmp_coeff(nbf), source=0.0d0)
 do i = 1, nmo-1, 1
  tmp_norm = norm(i)
  tmp_coeff = coeff1(:,i)

  do j = i+1, nmo, 1
   angle(j,i) = DOT_PRODUCT(coeff1(:,j),tmp_coeff)/(norm(j)*tmp_norm)
   angle(i,j) = angle(j,i)
  end do ! for j
 end do ! for i

 deallocate(coeff1, norm, tmp_coeff)
 write(6,'(A)') 'Angle matrix:'
 do i = 1, nmo, 1
  do j = 1, nmo, 1
   write(6,'(2I4,F10.6)') j, i, angle(j,i)
  end do ! for j
 end do ! for i

 ! 3rd step: find the most parallel orbital of an orbital
 allocate(target_pair_idx(nmo), source=0)
 do i = 1, nmo, 1
  tmp_idx = MAXLOC(angle(:,i))
  target_pair_idx(i) = tmp_idx(1)
 end do ! for i

 allocate(paired(nmo))
 paired = .false.
 forall(i=1:nmo, target_pair_idx(i)>0) paired(i) = .true.

 ! 4th step: iterate to achieve consistency
 tmp_idx = 0; iter = 0
 same_idx = .true.
 write(*,'(6I5)') target_pair_idx
 do while(same_idx)
  same_idx = .false.

  do i = 1, nmo-1, 1
   m = target_pair_idx(i)

   do j = i+1, nmo, 1
    n = target_pair_idx(j)

    if(n==m .and. n/=0) then
     if(angle(n,j) <= angle(n,i)) then
      tmp_norm = angle(n,j)
      angle(n,j) = -0.1d0
      tmp_idx = MAXLOC(angle(:,j))
      angle(n,j) = tmp_norm
      paired(j) = .true.; paired(tmp_idx(1)) = .true.
      target_pair_idx(j) = tmp_idx(1)
     else
      tmp_norm = angle(n,i)
      angle(n,i) = -0.1d0
      tmp_idx = MAXLOC(angle(:,i))
      angle(n,i) = tmp_norm
      paired(i) = .true.; paired(tmp_idx(1)) = .true.
      target_pair_idx(i) = tmp_idx(1)
      m = tmp_idx(1) ! udpate m
      write(*,'(3I5,2F10.6)') i, j, n, angle(n,j), angle(n,i)
      write(*,'(6I5)') target_pair_idx
     end if
     same_idx = .true.
    end if

   end do ! for j
  end do ! for i

  iter = iter + 1
  if(iter >= iter_max) exit
  if(iter == 4) stop
 end do ! for while

 write(6,'(/,A,I0)') 'iter = ', iter
 if(iter >= iter_max) then
  write(6,'(A)') 'ERROR in subroutine gvb_correct_pairs: iter exceeds max_cycle 1000.'
  stop
 end if

 nopen = COUNT(paired .eqv. .false.)
 write(6,'(A,I0)') 'detected nopen = ', nopen
 if(MOD(nmo-nopen,2) /= 0) then
  write(6,'(A)') 'ERROR in subroutine gvb_correct_pairs: I tried my best to&
                   & re-pair these slightly disordered orbitals.'
  write(6,'(A)') 'But nmo-nopen is not an even integer, which means possibly&
             & (1) the input idx1,idx2 may be incorrect; (2) re-pairing failed.'
  stop
 end if
 npair = (nmo - nopen)/2
 write(6,'(A,I0)') 'detected npair = ', npair

 allocate(ideal_idx(nmo), source=0)
 if(gau) then ! in Gaussian order
  forall(i=1:npair)
   ideal_idx(i) = nmo - i + 1
   ideal_idx(nmo-i+1) = i
  end forall
  if(nopen > 0) ideal_idx(npair+1:npair+nopen) = 0
 else ! in GAMESS order
  if(nopen > 0) ideal_idx(1:nopen) = 0
  forall(i=1:npair)
   ideal_idx(2*i-1+nopen) = 2*i + nopen
   ideal_idx(2*i + nopen) = 2*i-1 + nopen
  end forall
 end if
 if(ANY(ideal_idx/=target_pair_idx)) then
  if(gau) then
   write(6,'(A)') 'Warning: subroutine gvb_correct_pairs detected these&
                    & orbtials are not strictly in Gaussian GVB MO order.'
  else
   write(6,'(A)') 'Warning: subroutine gvb_correct_pairs detected these&
                    & orbtials are not strictly in GAMES GVB MO order.'
  end if
  write(6,'(A)') 'Trying to re-pair...'
 else
  write(6,'(A)') 'Congratulations! subroutine gvb_correct_pairs detected&
                   & the order of these GVB orbtials is probably correct.'
  return
 end if
 deallocate(ideal_idx)

 write(6,'(/,A)') 'Final pairs:'
 do i = 1, nmo, 1
  j = target_pair_idx(i)

  if(j /= 0) then
   write(6,'(3(I5,1X),F10.6)') i, i, j, angle(j,i)
  else
   write(6,'(2(I5,1X),A)') i, i, 'singly occupied'
  end if
 end do ! for i
 deallocate(angle)

 ! put new MO into the array coeff1
 allocate(coeff1(nbf,nif), source=0.0d0)
 coeff1(:,1:idx1-1) = coeff(:,1:idx1-1)
 coeff1(:,idx2+1:nif) = coeff(:,idx2+1:nif)
 m = idx1
 paired = .false.

 if(gau) then
  ! assign paired MOs into array coeff1 in Gaussian order
  do i = 1, nmo, 1
   if(paired(i)) cycle
   n = target_pair_idx(i)
   if(n == 0) cycle
   coeff1(:,m) = coeff(:,i+idx1-1)
   coeff1(:,idx2+idx1-m) = coeff(:,n+idx1-1)
   m = m + 1
   paired(i) = .true.; paired(n) = .true.
  end do ! for i
  ! assign singly occupied MOs into array coeff1
  do i = 1, nmo, 1
   if(target_pair_idx(i) == 0) then
    coeff1(:,m) = coeff(:,i+idx1-1)
    m = m + 1
   end if
  end do ! for i

 else ! in GAMESS order

  ! assign singly occupied MOs into array coeff1
  do i = 1, nmo, 1
   if(target_pair_idx(i) == 0) then
    coeff1(:,m) = coeff(:,i+idx1-1)
    m = m + 1
   end if
  end do ! for i
  ! assign paired MOs into array coeff1 in GAMESS order
  do i = 1, nmo, 1
   if(paired(i)) cycle
   n = target_pair_idx(i)
   if(n == 0) cycle
   coeff1(:,m) = coeff(:,i+idx1-1)
   coeff1(:,m+1) = coeff(:,n+idx1-1)
   m = m + 2
   paired(i) = .true.; paired(n) = .true.
  end do ! for i
 end if

 deallocate(paired, coeff, target_pair_idx)

 ! write the paired MOs into the new .fchk file
 call write_mo_into_fch(fchname, nbf, nif, 'a', coeff1)
 deallocate(coeff1)

 write(6,'(A)') 'Leave subroutine gvb_correct_pairs.'
 return
end subroutine gvb_correct_pairs

