! written by jxzou at 20180406: construct part of or all virtual orbitals using the PAO (projected atomic orbitals)
! updated by jxzou at 20191215: simplify code

! This subroutine is designed to be imported as a module by Python.

subroutine construct_vir(nbf, nif, idx, coeff, ovlp, new_coeff)
 implicit none
 integer :: i, j, nvir
 integer :: nbf, nif, idx
!f2py intent(in) :: nbf, nif, idx
 ! nbf: the number of basis functions
 ! nif: the number of MOs
 ! idx: the beginning index (Fortran convention) of the virtual MOs

 real(kind=8) :: coeff(nbf,nif), ovlp(nbf,nbf), new_coeff(nbf,nif)
!f2py depend(nbf,nif) coeff, new_coeff
!f2py depend(nbf) ovlp
!f2py intent(in) :: ovlp
!f2py intent(in,copy) :: coeff
!f2py intent(out) :: new_coeff

 real(kind=8), allocatable :: p(:,:), v(:,:), s1(:,:), ev(:), x(:,:)
 ! V: projected atomic orbitals (PAO)
 ! P: density matrix of atomic basis functions, sigma_i(Cui*Cvi)
 ! Note that the index i can be larger than nocc, in which case we only
 !  construct part of virtual orbitals.

 if(idx == nif+1) then
  write(6,'(/,A)') 'Warning in subroutine construct_vir: idx=nif+1 found.'
  write(6,'(A)') 'No need to construct virtual MOs.'
  new_coeff = coeff
  return
 end if

 ! Step 1: P = sigma_i(Cui*Cvi)
 allocate(p(nbf,nbf), source=0d0)
 call dgemm('N', 'T', nbf, nbf, idx-1, 1d0, coeff(1:nbf,1:idx-1), nbf, coeff(1:nbf,1:idx-1), nbf, 0d0, p, nbf)

 ! Step 2: V = 1 - PS
 allocate(v(nbf, nbf), source=0d0)
 forall(i = 1:nbf) v(i,i) = 1d0
 call dsymm('R', 'U', nbf, nbf, -1d0, ovlp, nbf, p, nbf, 1d0, v, nbf)
 deallocate(p)

 ! Step 3: S1 = (VT)SV
 allocate(s1(nbf,nbf), source=0d0)
 call calc_CTSC(nbf, nbf, v, ovlp, s1)

 ! Step 4: diagonalize S1 (note that S1 is symmetric) and get X
 allocate(ev(nbf))
 call diag_get_e_and_vec(nbf, s1, ev)
 nvir = nif - idx + 1
 forall(i = nbf-nvir+1:nbf) ev(i) = 1d0/DSQRT(ev(i))
 allocate(x(nbf,nvir), source=0d0)
 forall(i=1:nbf, j=1:nvir) x(i,j) = s1(i,nbf-nvir+j)*ev(nbf-nvir+j)
 deallocate(ev, s1)

 ! Step 5: get new virtual MO coefficients
 call dgemm('N', 'N', nbf, nvir, nbf, 1d0, v, nbf, x, nbf, 0d0, coeff(:,idx:nif), nbf)
 deallocate(x, v)
 new_coeff = coeff

 ! Step 6: check orthonormality
 allocate(x(nif,nif), source=0d0)
 call calc_CTSC(nbf, nif, coeff, ovlp, x)

 forall(i = 1:nif) x(i,i) = x(i,i) - 1d0
 x = DABS(x)
 write(6,'(/,A)') 'The orthonormality of Alpha MO after PAO construction:'
 write(6,'(A,F16.10)') 'maxv=', MAXVAL(x)
 write(6,'(A,F16.10)') 'abs_mean=', SUM(x)/DBLE(nif*nif)
 deallocate(x)
end subroutine construct_vir

