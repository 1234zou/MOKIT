! written by jxzou at 20180612
! updated by jxzou at 20200529: (1) change input ncore, npair, and nopen into idx1, idx2;
!                               (2) MOLCAS .INPORB file supported

! This is a subroutine/program used to align two sets of orbitals.
! It is useful when you have two sets of similar localized MOs, which are in different order.

program main
 implicit none
 integer :: i, idx1, idx2
 character(len=10) :: str
 character(len=240) :: fname1, fname2

 i = iargc()
 if(i /= 4) then
  write(6,'(/,A)') ' ERROR in subroutine align_orbitals: wrong command line arguments!'
  write(6,'(/,A)') ' Example 1(Gaussian): align_orbitals a.fch b.fchk 5 12'
  write(6,'(/,A)') ' Example 2(MOLCAS)  : align_orbitals a.INPORB b.INPORB 5 12'
  stop
 end if

 fname1 = ' '
 fname2 = ' '
 call getarg(1,fname1)
 call getarg(2,fname2)
 call getarg(3,str)
 read(str,*) idx1
 call getarg(4,str)
 read(str,*) idx2

 call align_orbitals(fname1, fname2, idx1, idx2)
end program main

subroutine align_orbitals(fname1, fname2, idx1, idx2)
 implicit none
 integer :: i, j, k, m, tmp_idx(1)
 integer :: nbf, nif, ncoeff
 integer, intent(in) :: idx1, idx2
 character(len=240),intent(in) :: fname1, fname2
 real(kind=8) :: tempv
 real(kind=8), allocatable :: coeff1(:,:), coeff2(:,:)
 real(kind=8), allocatable :: coeff(:), tmp_coeff(:)
 real(kind=8), allocatable :: diff(:)
 logical :: gau

 i = INDEX(fname1, '.fch', back=.true.)
 if(i /= 0) then
  gau = .true.
 else
  gau = .false.
 end if

 if(gau) then
  call read_nbf_and_nif_from_fch(fname1, nbf, nif)
 else
  call read_nbf_and_nif_from_orb(fname1, nbf, nif)
 end if

 allocate(coeff1(nbf,nif), coeff2(nbf,nif))
 if(gau) then
  call read_mo_from_fch(fname1, nbf, nif, 'a', coeff1)
  call read_mo_from_fch(fname2, nbf, nif, 'a', coeff2)
 else
  call read_mo_from_orb(fname1, nbf, nif, 'a', coeff1)
  call read_mo_from_orb(fname2, nbf, nif, 'a', coeff2)
 end if

 ncoeff = nbf*nif
 allocate(coeff(ncoeff), tmp_coeff(nbf))

 ! permute the order of MOs in array coeff2
 write(6,'(A)') 'In occ subspace:'
 do i = idx1, idx2, 1
  m = k - i + 1
  allocate(diff(m))
  diff = 0d0
  call get_mo_diff(nbf, coeff1(:,i), m, coeff2(:,i:k), diff)
  tmp_idx = MINLOC(diff)
  j = tmp_idx(1)
  tempv = diff(j)
  write(6,'(A,I4,2X,A,F11.5)') 'i=', i, 'diff=', tempv
  deallocate(diff)
  if(j /= 1) then
   tmp_coeff = coeff2(:,i)
   coeff2(:,i) = coeff2(:,j+i-1)
   coeff2(:,j+i-1) = tmp_coeff
  end if
 end do ! for i

 write(6,'(A)') 'In vir subspace:'
 k = nocc + npair
 do i = nocc+1, k, 1
  m = k - i + 1
  allocate(diff(m))
  diff = 0d0
  call get_mo_diff(nbf, coeff1(:,i), m, coeff2(:,i:k), diff)
  tmp_idx = MINLOC(diff)
  j = tmp_idx(1)
  tempv = diff(j)
  write(6,'(A,I4,2X,A,F11.5)') 'i=', i, 'diff=', tempv
  deallocate(diff)
  if(j /= 1) then
   tmp_coeff = coeff2(:,i)
   coeff2(:,i) = coeff2(:,j+i-1)
   coeff2(:,j+i-1) = tmp_coeff
  end if
 end do ! for i
 ! permute done

 ! output the new coeff2 into a new .fchk file

 deallocate(coeff1, coeff2)
 deallocate(coeff, tmp_coeff)
end subroutine align_orbitals

subroutine get_mo_diff(nbf, coeff1, k, coeff2, diff)
 implicit none
 integer :: i, j
 integer, intent(in) :: nbf, k
 real(kind=8) :: tempv1, tempv2
 real(kind=8) :: tmp_coeff(nbf)
 real(kind=8), intent(in) :: coeff1(nbf), coeff2(nbf,k)
 real(kind=8), intent(out) :: diff(k)

 do i = 1, k, 1
  tmp_coeff = coeff1 - coeff2(:,i)
  forall(j = 1:nbf)
   tmp_coeff(j) = tmp_coeff(j)*tmp_coeff(j)
  end forall
  tempv1 = SUM(tmp_coeff)
  tmp_coeff = coeff1 + coeff2(:,i)
  forall(j = 1:nbf)
   tmp_coeff(j) = tmp_coeff(j)*tmp_coeff(j)
  end forall
  tempv2 = SUM(tmp_coeff)
  diff(i) = min(tempv1,tempv2)
 end do ! for i
end subroutine get_mo_diff

