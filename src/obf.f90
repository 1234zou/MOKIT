! written by jxzou at 20221206: orbital-based fragmentation method
! This can be viewed as the 2nd generation of GEBF
! Warning: still some bug in this file

module obf   ! orbital-based fragmentation
 use population, only: nbf, nif, nmo, i1, i2, mo_cluster, mo_dis
 implicit none
 integer :: n_old, n_new
 integer :: n_prim ! number of primitive MO clusters
 integer :: n_tot  ! total number of MO clusters (primitive + derivative)
 integer :: na, nb, nopen ! number of alpha, beta and singly occupied orbitals
 integer, allocatable :: icoeff0(:), icoeff1(:), icoeff(:)
 integer, allocatable :: label0(:,:), label1(:,:)
 real(kind=8) :: e_tot = 0d0       ! total electronic energy
 real(kind=8), parameter :: dis_thres0 = 4d0  ! Angstrom
 real(kind=8), allocatable :: cluster_e(:) ! size n_tot
 character(len=240) :: fchname = ' '
 type(mo_cluster), allocatable :: cluster0(:), cluster1(:), cluster(:)

contains

! '=' cannot be used for copying type, so use this subroutine
subroutine copy_type_mo_cluster(old_c, new_c)
 implicit none
 integer :: nocc, nvir
 type(mo_cluster), intent(in) :: old_c
 type(mo_cluster), intent(out) :: new_c

 nocc = old_c%nocc
 new_c%nocc = nocc
 if(nocc > 0) allocate(new_c%occ_idx(nocc), source=old_c%occ_idx)

 nvir = old_c%nvir
 new_c%nvir = nvir
 if(nvir > 0) allocate(new_c%vir_idx(nvir), source=old_c%vir_idx)
end subroutine copy_type_mo_cluster

! '=' cannot be used for copying type, so use this subroutine
subroutine copy_type_mo_clusters(n, old_c, new_c)
 implicit none
 integer :: i
 integer, intent(in) :: n
 type(mo_cluster), intent(in) :: old_c(n)
 type(mo_cluster), intent(out) :: new_c(n)

 do i = 1, n, 1
  call copy_type_mo_cluster(old_c(i), new_c(i))
 end do ! for i
end subroutine copy_type_mo_clusters

! generate the MO-cluster for each occupied MO
subroutine gen_mo_cluster_per_mo(dis_thres)
 implicit none
 integer :: i, j, k, m, nocc
 real(kind=8), intent(in) :: dis_thres

 if(dis_thres<1d-2 .or. dis_thres>99d0) then
  write(6,'(/,A,F7.3)') 'ERROR in subroutine gen_mo_cluster_per_mo: invalid dis&
                        &_thres=', dis_thres
  stop
 end if

 nmo = i2 - i1 + 1
 allocate(cluster0(nmo))

 do i = 1, nmo, 1
  cluster0(i)%nvir = 0
  k = i - 1 + i1
  nocc = COUNT(mo_dis(:,k) < dis_thres)
  cluster0(i)%nocc = nocc
  allocate(cluster0(i)%occ_idx(nocc))
  m = 1
  cluster0(i)%occ_idx(1) = k

  do j = i1, i2, 1
   if(j == k) cycle
   if(mo_dis(j,k) < dis_thres) then
    m = m + 1
    cluster0(i)%occ_idx(m) = j
   end if
  end do ! for j
 end do ! for i

 deallocate(mo_dis)

 ! sort orbital indices in type cluster0
 do i = 1, nmo, 1
  nocc = cluster0(i)%nocc
  call sort_int_array(nocc, cluster0(i)%occ_idx, .true.)
 end do ! for i

 !write(6,'(A)') 'cluster0:'
 !do i = 1, nmo, 1
 ! write(6,'(I2,A1,23I4)') cluster0(i)%nocc,':',cluster0(i)%occ_idx
 !end do ! for i

 if(ANY(cluster0(:)%nocc == 0)) then
  write(6,'(/,A)') 'ERROR in subroutine gen_mo_cluster_per_mo: some nocc=0.'
  write(6,'(A)') 'This is impossible.'
  write(6,'(A,/,23I4)') 'cluster0(:)%nocc=', cluster0(:)%nocc
  stop
 end if
end subroutine gen_mo_cluster_per_mo

! delete primitive MO clusters which are embraced by other primitive clusters
! cluster0 -> cluster1, icoeff0 -> icoeff1
subroutine delete_embraced_cluster()
 implicit none
 integer :: i, j, k, m, n, p, q, nocc1, nocc2
 integer, allocatable :: idx1(:), idx2(:)
 logical :: deleted
 logical, allocatable :: del(:)

 if(n_old == 0) then
  write(6,'(/,A)') 'ERROR in subroutine delete_embraced_cluster: n_old=0.'
  write(6,'(A)') 'Something must be wrong.'
  stop
 end if
 allocate(del(n_old), source=.false.)

 do i = 1, n_old-1, 1
  if(del(i)) cycle
  nocc1 = cluster0(i)%nocc
  if(nocc1 == 0) then
   del(i) = .true.
   cycle
  end if

  do j = i+1, n_old, 1
   if(del(j)) cycle
   nocc2 = cluster0(j)%nocc
   if(nocc2 == 0) then
    del(j) = .true.
    cycle
   end if

   if(nocc1 < nocc2) then
    k = nocc1; m = nocc2
    allocate(idx1(k), source=cluster0(i)%occ_idx)
    allocate(idx2(m), source=cluster0(j)%occ_idx)
   else ! nocc1 >= nocc2
    k = nocc2; m = nocc1
    allocate(idx1(k), source=cluster0(j)%occ_idx)
    allocate(idx2(m), source=cluster0(i)%occ_idx)
   end if
   deleted = .true.

   do p = 1, k, 1 ! k >= m
    n = idx1(p)
    do q = 1, m, 1
     if(ALL(idx2 /= n)) then
      deleted = .false.
      exit
     end if
    end do ! for q
   end do ! for p

   deallocate(idx1, idx2)

   if(deleted) then
    if(nocc1 < nocc2) then
     del(i) = .true.
     exit
    else ! nocc1 >= nocc2
     del(j) = .true.
     if(nocc1 == nocc2) then
      icoeff0(i) = icoeff0(i) + icoeff0(j)
      icoeff0(j) = 0
     end if
    end if
   end if

  end do ! for j
 end do ! for i

 n_new = COUNT(del .eqv. .false.)
 if(n_new == 0) then
  write(6,'(/,A)') 'ERROR in subroutine delete_embraced_cluster: all MO cluster&
                   &s are deleted.'
  write(6,'(A)') 'Something must be wrong.'
  stop
 end if

 j = 0
 allocate(cluster1(n_new))
 allocate(icoeff1(n_new), source=0)
 p = size(label0,1)
 allocate(label1(p,n_new), source=0)

 do i = 1, n_old, 1
  if(del(i)) cycle
  j = j + 1
  icoeff1(j) = icoeff0(i)
  label1(:,j) = label0(:,i)
  call copy_type_mo_cluster(cluster0(i), cluster1(j))
 end do ! for i

 deallocate(del, cluster0, icoeff0, label0)

 write(6,'(A)') 'merged MO clusters:'
 do i = 1, n_new, 1
  write(6,'(I4,A1,23I4)') icoeff1(i),':',cluster1(i)%occ_idx
 end do ! for i
end subroutine delete_embraced_cluster

! generate primitive MO clusters
subroutine gen_prim_cluster(dis_thres)
 implicit none
 integer :: i
 real(kind=8), intent(in) :: dis_thres

 call gen_mo_cluster_per_mo(dis_thres)
 n_old = nmo
 allocate(icoeff0(n_old), source=1)
 allocate(label0(1,n_old))
 forall(i = 1:n_old) label0(1,i) = i

 call delete_embraced_cluster()
 deallocate(icoeff1, label1) ! useless for primitive subsystems

 if(n_new == 1) then
  write(6,'(/,A)') 'ERROR in subroutine gen_prim_cluster: only one primitive&
                  & subsystem is generated'
  write(6,'(A)') 'and it is equal to the whole system. So there is no need to&
                & use obf.'
  stop
 end if

 n_prim = n_new
 n_tot = n_new
 allocate(cluster(n_new))
 call copy_type_mo_clusters(n_new, cluster1, cluster)
 deallocate(cluster1)

 allocate(icoeff(n_tot), source=1)
end subroutine gen_prim_cluster

! generate derivative MO clusters from primitive MO clusters
subroutine gen_deri_cluster()
 implicit none
 integer :: i, j, k, m, n, p, n1, n2, n_deri
 integer, external :: comb, label2idx
 integer, allocatable :: itmp(:), itmp1(:)
 logical :: alive, cycle_k

 ! copy primitive MO cluster from type cluster to cluster2
 n_old = n_tot
 allocate(cluster1(n_old), icoeff1(n_old), label1(1,n_old))
 call copy_type_mo_clusters(n_old, cluster, cluster1)
 icoeff1 = icoeff
 forall(i = 1:n_old) label1(1,i) = i
 write(6,'(A)') 'derivative MO clusters generation:'

! According to the inclusion-exclusion principle, the upper limit of loop is
! n_prim. But it usually terminates after several loops, so we need to check
! whether each MO occurs exactly once at the end of each loop. If yes, then
! exit the loop, i.e. all derivative MO clusters are generated.
 do i = 2, n_prim, 1
  write(6,'(A,I0)') 'i=', i
  n_deri = n_old*(n_old-1)/2
  ! the number of derivative clusters of the current generation

  allocate(cluster0(n_deri), label0(i,n_deri), icoeff0(n_deri))
  allocate(itmp(2*(i-1)), itmp1(i-1))
  label0 = 0; icoeff0 = 0; m = 0

  do j = 1, n_old-1, 1
   n1 = cluster1(j)%nocc
   do k = j+1, n_old, 1
    n2 = cluster1(k)%nocc

    ! if the intersection does not belong to this generation, cycle
    ! e.g. (AxB)x(CxD) = AxBxCxD is the 4-th generation, not the 3rd
    call find_union(i-1, label1(:,j), label1(:,k), itmp)
    n = COUNT(itmp > 0)
    if(n > i) then
     cycle
    else if(n < i) then
     write(6,'(A)') 'ERROR in subroutine gen_deri_cluster: n<i is impossible.'
     write(6,'(A,3I3,A,10I3)') 'j,k,n=', j,k,n, ', itmp=', itmp
     write(6,'(A,10I3)') 'label1(:,j)=',label1(:,j)
     write(6,'(A,10I3)') 'label1(:,k)=',label1(:,k)
     stop
    end if

    ! if the intersection has appeared before, cycle
    cycle_k = .false.
    do n = 1, m, 1
     if(ALL(label0(:,n) == itmp(1:i))) then
      cycle_k = .true.
      exit
     end if
    end do ! for n
    if(cycle_k) cycle

    ! if any parent set of this intersection is null, cycle
    do n = 1, i, 1
     cycle_k = .true.

     if(n == 1) then
      itmp1 = itmp(2:i)
     else if(n == i) then
      itmp1 = itmp(1:i-1)
     else ! 1 < n < i
      itmp1(1:n-1) = itmp(1:n-1)
      itmp1(n:i-1) = itmp(n+1:i)
     end if

     do p = 1, n_old, 1
      if(ALL(label1(:,p) == itmp1)) then
       cycle_k = .false.
       exit
      end if
     end do ! for p

     if(cycle_k) exit
    end do ! for n
    if(cycle_k) cycle
    ! done check parent sets

    m = m + 1
    label0(:,m) = itmp(1:i)
    call find_intersec(n1, cluster1(j)%occ_idx, n2, cluster1(k)%occ_idx, &
                       cluster0(m))
    if(cluster0(m)%nocc > 0) then
     icoeff0(m) = icoeff1(j)*icoeff1(k)
     if(MOD(i,2) == 0) icoeff0(m) = -icoeff0(m)
    end if
   end do ! for k
  end do ! for j

  deallocate(itmp, itmp1, cluster1, icoeff1, label1)
  if(m == 0) then ! no intersection in this generation
   deallocate(cluster0, icoeff0, label0)
   exit
  end if

  n_old = m   ! will be used in empty_embraced_cluster
  call delete_embraced_cluster()
  ! cluster0 -> cluster1, icoeff0 -> icoeff1, label0 -> label1

  call append_cluster1_to_cluster()

  n_old = n_new   ! update n_old for next cycle
  if(ALL(cluster1(:)%nocc < 2)) exit
 end do ! for i

 call check_occur_once(alive)

 if(.not. alive) then
  write(6,'(/,A)') 'ERROR in subroutine gen_deri_cluster: some MO cluster does&
                   & not occur once.'
  call merge_mo_cluster()

  do i = 1, n_tot, 1
   write(6,'(I3,A1,23I4)') icoeff(i),':',cluster(i)%occ_idx
  end do ! for i
  stop
 end if
end subroutine gen_deri_cluster

! find the union set of two integer arrays
! Note: all elements in integer arrays a1 and a2 must be positive
subroutine find_union(n, a1, a2, a3)
 implicit none
 integer :: i, j
 integer, intent(in) :: n
 integer, intent(in) :: a1(n), a2(n)
 integer, intent(out) :: a3(2*n)

 a3(1:n) = a1; a3(n+1:) = 0; j = n

 do i = 1, n, 1
  if(ANY(a1 == a2(i))) cycle
  j = j + 1
  a3(j) = a2(i)
 end do ! for i

 ! sort a3(1:j)
 call sort_int_array(j, a3(1:j), .true.)
end subroutine find_union

! append type cluster1 to type cluster, the latter will be enlarged
! also enlarge the integer array icoeff
subroutine append_cluster1_to_cluster()
 implicit none
 integer :: n
 integer, allocatable :: tmp_i(:)
 type(mo_cluster), allocatable :: tmp_c(:)

 allocate(tmp_c(n_tot))
 call copy_type_mo_clusters(n_tot, cluster, tmp_c)
 deallocate(cluster)

 n = n_tot + n_new
 allocate(cluster(n))
 call copy_type_mo_clusters(n_tot, tmp_c, cluster(1:n_tot))
 deallocate(tmp_c)

 allocate(tmp_i(n_tot), source=icoeff)
 deallocate(icoeff)
 allocate(icoeff(n))
 icoeff(1:n_tot) = tmp_i
 deallocate(tmp_i)

 icoeff(n_tot+1:) = icoeff1
 call copy_type_mo_clusters(n_new, cluster1, cluster(n_tot+1:))

 n_tot = n
end subroutine append_cluster1_to_cluster

! find intersection of two integer arrays, save the result into type clus
subroutine find_intersec(n1, a1, n2, a2, clus)
 implicit none
 integer :: i, j, k, m
 integer, intent(in) :: n1, n2
 integer, intent(in) :: a1(n1), a2(n2)
 integer, allocatable :: occ_idx(:)
 logical, allocatable :: skip1(:), skip2(:)
 type(mo_cluster), intent(out) :: clus

 if(n1==0 .or. n2==0) return ! the intersection must be empty
 allocate(skip1(n1), source=.false.)
 allocate(skip2(n2), source=.false.)
 k = 0

 do i = 1, n1, 1
  if(skip1(i)) cycle
  m = a1(i)

  do j = 1, n2, 1
   if(skip2(j)) cycle

   if(m == a2(j)) then
    skip1(i) = .true.
    skip2(j) = .true.
    k = k + 1
    if(k > 1) then
     occ_idx = [occ_idx, m]
    else
     allocate(occ_idx(1))
     occ_idx(1) = m
    end if
    exit
   end if
  end do ! for j
 end do ! for i

 deallocate(skip1, skip2)
 clus%nocc = k
 if(k > 0) then
  allocate(clus%occ_idx(k), source=occ_idx)
  deallocate(occ_idx)
 end if
end subroutine find_intersec

subroutine check_occur_once(alive)
 implicit none
 integer :: i, j, k
 logical, intent(out) :: alive

 alive = .true.

 do i = i1, i2, 1
  k = 0

  do j = 1, n_tot, 1
   k = k + icoeff(j)*COUNT(cluster(j)%occ_idx == i)
  end do ! for j

  if(k /= 1) then
   alive = .false.
   return
  end if
 end do ! for i
end subroutine check_occur_once

! merger primitive and derivative MO clusters, n_tot and type cluster will be
! updated if there exists any derivative MO cluster which is identical to some
! primitive MO cluster
subroutine merge_mo_cluster()
 implicit none
 integer :: i, j, nocc
 integer, allocatable :: tmp_o(:)
 logical, allocatable :: del(:)
 type(mo_cluster), allocatable :: tmp_c(:)

 if(n_tot < 2) return
 allocate(del(n_tot), source=.false.)

 do i = 1, n_tot-1, 1
  if(del(i)) cycle
  nocc = cluster(i)%nocc
  allocate(tmp_o(nocc), source=cluster(i)%occ_idx)

  do j = i+1, n_tot, 1
   if(del(j)) cycle
   if(cluster(j)%nocc /= nocc) cycle

   if(ALL(tmp_o == cluster(j)%occ_idx)) then
    icoeff(i) = icoeff(i) + icoeff(j)
    del(j) = .true.

    if(icoeff(i) == 0) then
     del(i) = .true.
     exit
    end if
   end if
  end do ! for j

  deallocate(tmp_o)
 end do ! for i

 j = COUNT(del .eqv. .false.)
 if(j == 0) then
  write(6,'(/,A)') 'ERROR in subroutine merge_mo_cluster: all MO clusters are&
                  & deleted.'
  write(6,'(A)') 'Something must be wrong.'
  stop
 end if

 if(j == n_tot) then ! nothing to be updated
  deallocate(del)
  return
 end if

 ! update type cluster
 allocate(tmp_c(j))
 j = 0
 do i = 1, n_tot, 1
  if(del(i)) cycle
  j = j + 1
  call copy_type_mo_cluster(cluster(i), tmp_c(j))
 end do ! for i

 deallocate(cluster)
 allocate(cluster(j))
 call copy_type_mo_clusters(j, tmp_c, cluster)
 deallocate(tmp_c)

 ! update icoeff
 allocate(tmp_o(j), source=0)
 j = 0
 do i = 1, n_tot, 1
  if(del(i)) cycle
  j = j + 1
  tmp_o(j) = icoeff(i)
 end do ! for i

 deallocate(icoeff, del)
 allocate(icoeff(j), source=tmp_o)
 deallocate(tmp_o)
 n_tot = j
end subroutine merge_mo_cluster

! Add paired active virtual MO indices into the type cluster
! This is because we always assume the input fchname includes paired active or-
!  bitals, e.g. paired localized UNO(i.e. associated rotatied UNO), or GVB orb-
!  itals
subroutine add_paired_vir2cluster()
 implicit none
 integer :: i, j, k, nocc
 integer, allocatable :: idx(:)

 call read_na_and_nb_from_fch(fchname, na, nb)
 nopen = na - nb
 k = 2*i2 + 1 + nopen

 do i = 1, n_tot, 1
  nocc = cluster(i)%nocc
  allocate(idx(nocc), source=cluster(i)%occ_idx)
  deallocate(cluster(i)%occ_idx)
  cluster(i)%nocc = 2*nocc
  allocate(cluster(i)%occ_idx(2*nocc))
  cluster(i)%occ_idx(1:nocc) = idx
  deallocate(idx)
  forall(j = 1:nocc)
   cluster(i)%occ_idx(nocc+j) = k - cluster(i)%occ_idx(nocc-j+1)
  end forall
 end do ! for i
end subroutine add_paired_vir2cluster

! generate all .fch files with active orbitals permuted near HONO or LUNO
subroutine gen_permute_fch()
 implicit none
 integer :: i
 character(len=240) :: proname, new_fch
 real(kind=8), allocatable :: mo(:,:), mo1(:,:)

 proname = ' '; new_fch = ' '
 i = index(fchname,'.fch', back=.true.)
 proname = fchname(1:i-1)

 allocate(mo(nbf,nif), mo1(nbf,nif))
 call read_mo_from_fch(fchname, nbf, nif, 'a', mo)

 do i = 1, n_tot, 1
  write(new_fch,'(A,I0,A)') TRIM(proname), i, '.fch'
  call copy_file(fchname, new_fch, .false.)
  call permute_mo_in_sub_cluster(nbf, nif, mo, cluster(i), mo1)
  call write_mo_into_fch(new_fch, nbf, nif, 'a', mo1)
 end do ! for i

 deallocate(mo, mo1)
end subroutine gen_permute_fch

end module obf

program main
 use population, only: mulliken_pop_of_mo
 use obf, only: dis_thres0, n_tot, icoeff, cluster, fchname, gen_prim_cluster,&
  gen_deri_cluster, merge_mo_cluster, add_paired_vir2cluster, gen_permute_fch
 implicit none
 integer :: i, ibegin, iend
 ! ibegin: the beginning index of active occupied MOs
 ! iend: the final index of active occupied MOs, singly occupied not included
 integer(kind=4) :: hostnm
 character(len=5) :: str = ' '
 character(len=24) :: hostname, data_string
 real(kind=8) :: dis_thres

 i = iargc()
 if(i<3 .or. i>4) then
  write(6,'(/,A)') 'ERROR in program obf: wrong command line arguments.'
  write(6,'(A)') 'Example 1: obf tetracene_uno_asrot.fch 52 60'
  write(6,'(A)') 'Example 2: obf tetracene_uno_asrot.fch 52 60 4.0'
  write(6,'(A,/)') 'Note: do not include any singly occupied orbital.'
  stop
 end if

 fchname = ' '
 call getarg(1, fchname)
 call getarg(2, str)
 read(str,*) ibegin
 call getarg(3, str)
 read(str,*) iend
 if(i == 4) then
  call getarg(4, str)
  read(str,*) dis_thres
 else
  dis_thres = dis_thres0
 end if

 call fdate(data_string)
 write(6,'(A)') 'Obf program begins at '//TRIM(data_string)
 i = hostnm(hostname)
 write(6,'(A)') 'HOST '//TRIM(hostname)
 write(6,'(2(A,I0),A,F5.2)') 'ibegin=',ibegin,', iend=',iend,', dis_thres=',&
                              dis_thres

 ! perform Mulliken population for each MO, find the centers of each MO,
 ! and calculate the distances between any two MOs.
 call mulliken_pop_of_mo(fchname, ibegin, iend)

 call gen_prim_cluster(dis_thres) ! generate primitive MO clusters

 call gen_deri_cluster() ! generate derivative MO clusters

 ! Note: here all icoeff(i) can only be either -1 or 1
 write(6,'(A,I3)') 'Before merge, n_tot=', n_tot
 do i = 1, n_tot, 1
  write(6,'(I2,A1,23I4)') icoeff(i),':',cluster(i)%occ_idx
 end do ! for i

 ! merge primitive MO clusters and derivative MO clusters
 ! now some icoeff(i) may be -2, 2, etc
 call merge_mo_cluster()

 call add_paired_vir2cluster()

 write(6,'(A,I3)') 'After merge and paired, n_tot=', n_tot
 do i = 1, n_tot, 1
  write(6,'(I2,A1,23I4)') icoeff(i),':',cluster(i)%occ_idx
 end do ! for i
 stop

 ! generate all .fch files with active orbitals permuted near HONO or LUNO
 call gen_permute_fch()

 ! generate all automr input files(.gjf) and submit one by one
 call gen_automr_gjf_and_submit()
 deallocate(cluster)

 call read_cluster_e_from_out()

 call fdate(data_string)
 write(6,'(/,A)') 'Obf program ends at '//TRIM(data_string)
end program main

! Compute Cn,k. E.g. C4,2 = 6
function comb(n,k)
 implicit none
 integer :: i, k1
 integer :: numerator, denominator, comb
 integer, intent(in) :: n, k

 numerator = 1
 denominator = 1
 if(2*k <= n) then
  k1 = k     ! Cn,k
 else
  k1 = n - k ! Cn,(n-k)
 end if

 do i = n, n-k1+1, -1
  numerator = numerator*i
 end do ! for i

 do i = 1, k1, 1
  denominator = denominator*i
 end do ! for i

 comb = numerator/denominator
end function comb

! initialize the integer array label
! Example: when np=4, m=2, label contains 1,2, 1,3, 1,4, 2,3, 2,4, 3,4, n=6
subroutine init_label(np, m, n, label)
 implicit none
 integer, intent(in) :: np, m, n
 integer, intent(out) :: label(m,n)

 select case(m)
 case(2)
  call init_label2(np, n, label)
 case(3)
  call init_label3(np, n, label)
 case(4)
  call init_label4(np, n, label)
 case(5)
  call init_label5(np, n, label)
 case(6)
  call init_label6(np, n, label)
 case default
  write(6,'(/,A,I0)') 'ERROR in subroutine init_label: m=', m
  write(6,'(A)') 'Currently only m=2,3,4,5,6 are supported.'
  stop
 end select
end subroutine init_label

subroutine init_label2(np, n, label)
 implicit none
 integer :: i, j, k
 integer, intent(in) :: np, n
 integer, intent(out) :: label(2,n)

 k = 0
 do i = 1, np-1, 1
  do j = i+1, np, 1
   k = k + 1
   label(:,k) = [i,j]
  end do ! for j
 end do ! for i
end subroutine init_label2

subroutine init_label3(np, n, label)
 implicit none
 integer :: i, j, k, m
 integer, intent(in) :: np, n
 integer, intent(out) :: label(3,n)

 m = 0
 do i = 1, np-2, 1
  do j = i+1, np-1, 1
   do k = j+1, np, 1
    m = m + 1
    label(:,m) = [i,j,k]
   end do ! for k
  end do ! for j
 end do ! for i
end subroutine init_label3

subroutine init_label4(np, n, label)
 implicit none
 integer :: i, j, k, m, p
 integer, intent(in) :: np, n
 integer, intent(out) :: label(4,n)

 p = 0
 do i = 1, np-3, 1
  do j = i+1, np-2, 1
   do k = j+1, np-1, 1
    do m = k+1, np, 1
     p = p + 1
     label(:,p) = [i,j,k,m]
    end do ! for m
   end do ! for k
  end do ! for j
 end do ! for i
end subroutine init_label4

subroutine init_label5(np, n, label)
 implicit none
 integer :: i, j, k, m, p, q
 integer, intent(in) :: np, n
 integer, intent(out) :: label(5,n)

 q = 0
 do i = 1, np-4, 1
  do j = i+1, np-3, 1
   do k = j+1, np-2, 1
    do m = k+1, np-1, 1
     do p = m+1, np, 1
      q = q + 1
      label(:,q) = [i,j,k,m,p]
     end do ! for p
    end do ! for m
   end do ! for k
  end do ! for j
 end do ! for i
end subroutine init_label5

subroutine init_label6(np, n, label)
 implicit none
 integer :: i, i1, i2, i3, i4, i5, i6
 integer, intent(in) :: np, n
 integer, intent(out) :: label(6,n)

 i = 0
 do i1 = 1, np-5, 1
  do i2 = i1+1, np-4, 1
   do i3 = i2+1, np-3, 1
    do i4 = i3+1, np-2, 1
     do i5 = i4+1, np-1, 1
      do i6 = i5+1, np, 1
       i = i + 1
       label(:,i) = [i1,i2,i3,i4,i5,i6]
      end do ! for i6
     end do ! for i5
    end do ! for i4
   end do ! for i3
  end do ! for i2
 end do ! for i1
end subroutine init_label6

! given a set of integers, return its corresponding index
! e.g. np=6, a=[2,3,4] => idx=11
! see subroutine init_label and you will know why the index is calculated in
! this way
function label2idx(np, n, a) result(idx)
 implicit none
 integer :: i, j, k, m, p, q, r, idx
 integer, intent(in) :: np, n
 integer, intent(in) :: a(n)

 i = a(1)
 if(n > 1) j = a(2)

 select case(n)
 case(1)
  idx = i
 case(2)
  idx = (2*np-i)*(i-1)/2 + j - i
 case(3)
  k = np - i
  idx = (np*(np-1)*(np-2) - (k+1)*k*(k-1))/6 + (2*np-i-j)*(j-i-1)/2 + a(3) - j
 case(4)
  k = a(3); m = np - i; p = np - j
  idx = (np*(np-1)*(np-2)*(np-3) - (m+1)*m*(m-1)*(m-2))/24 + &
        (m*(m-1)*(m-2) - (p+1)*p*(p-1))/6 + (2*np-j-k)*(k-j-1)/2 + a(4) - k
 case(5)
  k = a(3); m = a(4); p = np - i; q = np - j; r = np - k
  idx = (np*(np-1)*(np-2)*(np-3)*(np-4) - (p+1)*p*(p-1)*(p-2)*(p-3))/120 + &
        (p*(p-1)*(p-2)*(p-3) - (q+1)*q*(q-1)*(q-2))/24 + &
        (q*(q-1)*(q-2) - (r+1)*r*(r-1))/6 + (2*np-k-m)*(m-k-1)/2 + a(5) - m
 case default
  write(6,'(/,A,I0)') 'ERROR in function label2idx: n=', n
  write(6,'(A)') 'Currently only n=1,2,3,4,5 are supported.'
  stop
 end select
end function label2idx

! sort an integer array by ascending order
subroutine sort_int_array(n, a, ascending)
 implicit none
 integer :: i, j, k, m
 integer, intent(in) :: n
 integer, intent(inout) :: a(n)
 logical, intent(in) :: ascending

 if(n == 1) return

 if(ascending) then
  do i = 1, n-1, 1
   k = a(i)
   do j = i+1, n, 1
    if(k > a(j)) then
     a(i) = a(j); m = a(j); a(j) = k; k = m
    end if
   end do ! for j
  end do ! for i

 else ! descending order
  do i = 1, n-1, 1
   k = a(i)
   do j = i+1, n, 1
    if(k < a(j)) then
     a(i) = a(j); m = a(j); a(j) = k; k = m
    end if
   end do ! for j
  end do ! for i
 end if
end subroutine sort_int_array

! permute active MOs of a sub-cluster
subroutine permute_mo_in_sub_cluster(nbf, nif, mo, clus, mo1)
 use obf, only: na, nb, mo_cluster
 implicit none
 integer :: i, j, m, ndb, ivir
 integer, intent(in) :: nbf, nif
 real(kind=8), intent(in) :: mo(nbf,nif)
 real(kind=8), intent(out) :: mo1(nbf,nif)
 type(mo_cluster), intent(in) :: clus

 mo1 = 0d0
 if(nb < na) mo1(:,nb+1:na) = mo(:,nb+1:na) ! copy singly occupied MOs

 ! copy doubly occupied MOs and unmoved active occupied MOs
 ndb = MINVAL(clus%occ_idx) - 1
 if(ndb > 0) mo1(:,1:ndb) = mo(:,1:ndb)

 ivir = MAXVAL(clus%occ_idx) + 1
 if(ivir <= nif) mo1(:,ivir:) = mo(:,ivir:) ! copy inactive virtual MOs

 j = ndb
 m = clus%nocc/2

 do i = ndb+1, nb, 1
  if(ANY(clus%occ_idx(1:m)==i)) cycle
  j = j + 1
  mo1(:,j) = mo(:,i)
 end do ! for i

 ! permute active occupied and unoccupied MOs, near HONO and LUNO
 forall(i = 1:m)
  mo1(:,j+i) = mo(:,clus%occ_idx(i))
  mo1(:,na+i) = mo(:,clus%occ_idx(m+i))
 end forall

 j = na + m

 do i = na+1, ivir-1, 1
  if(ANY(clus%occ_idx==i)) cycle
  j = j + 1
  mo1(:,j) = mo(:,i)
 end do ! for i
end subroutine permute_mo_in_sub_cluster

! generate all automr input files(.gjf) and submit one by one
! TODO: if file 'hosts' exists, then run on multiple nodes.
subroutine gen_automr_gjf_and_submit()
 use obf, only: fchname, n_tot, cluster
 implicit none
 integer :: i, ne, fid
 integer, parameter :: nproc = 48, mem = 200
 character(len=240) :: proname, gjfname, new_fch

 i = index(fchname, '.fch', back=.true.)
 proname = fchname(1:i-1)

 do i = 1, n_tot, 1
  write(gjfname,'(A,I0,A)') TRIM(proname),i,'.gjf'
  write(new_fch,'(A,I0,A)') TRIM(proname),i,'.fch'
  open(newunit=fid,file=TRIM(gjfname),status='replace')
  write(fid,'(A,I0)') '%nprocshared=',nproc
  write(fid,'(A,I0,A)') '%mem=',mem,'GB'
  ne = cluster(i)%nocc
  write(fid,'(2(A,I0),A)') '#p CASCI(',ne,',',ne,')/cc-pVDZ'
  write(fid,'(/,A)') "mokit{ist=5,readno='"//TRIM(new_fch)//"'}"
  close(fid)
 end do ! for i

 do i = 1, n_tot, 1
  write(gjfname,'(A,I0,A)') TRIM(proname),i,'.gjf'
  call submit_automr_job(gjfname)
 end do ! for i
end subroutine gen_automr_gjf_and_submit

! read the electronic energy of each MO cluster from output files
subroutine read_cluster_e_from_out()
 use obf, only: fchname, n_tot, e_tot, cluster_e, icoeff
 implicit none
 integer :: i, k, fid, system
 character(len=240) :: buf, proname, proname1, outname
 character(len=500) :: longbuf = ' '

 i = index(fchname, '.fch', back=.true.)
 proname = fchname(1:i-1)
 allocate(cluster_e(n_tot), source=0d0)

 do i = 1, n_tot, 1
  write(proname1,'(A,I0)') TRIM(proname),i
  outname = TRIM(proname1)//'_CASCI.out'

  open(newunit=fid,file=TRIM(outname),status='old',position='append')
  BACKSPACE(fid)
  read(fid,'(A)') buf
  close(fid)

  if(buf(1:7) /= 'CASCI E') then
   write(6,'(A)') "ERROR in subroutine gen_automr_gjf_and_submit: 'CASCI E' not&
                 & found at the end of "//TRIM(outname)
   stop
  end if

  k = index(buf, '=')
  read(buf(k+1:),*) cluster_e(i)
  longbuf = 'tar -zcf '//TRIM(proname1)//'.tar.gz '//TRIM(proname1)//'.* '//&
            TRIM(proname1)//'_* --remove-files'
  k = system(longbuf)
  if(k /= 0) write(6,'(A)') 'ERROR in subroutine read_cluster_e_from_out: faile&
                            &d to create package '//TRIM(proname1)//'.tar.gz'
 end do ! for i

 e_tot = DOT_PRODUCT(cluster_e, DBLE(icoeff))
 write(6,'(/,A,F18.8,A)') 'E_tot = ',e_tot,' a.u.'
 deallocate(cluster_e)
end subroutine read_cluster_e_from_out

