! written by jxzou at 20221206: orbital-based fragmentation method
! This can be viewed as the 2nd generation of GEBF
! Warning: still some bug in this file

module obf   ! orbital-based fragmentation
 use population, only: nbf, nif, nmo, i1, i2, mo_cluster, mo_dis
 implicit none
 integer :: n_old, n_new
 integer :: n_prim ! number of primitive MO clusters
 integer :: n_tot  ! total number of MO clusters (primitive + derivative)
 integer :: na, nb, nopen ! number of alpha/beta/singly occupied orbitals
 integer, allocatable :: icoeff0(:), icoeff1(:), icoeff(:)
 integer, allocatable :: label0(:,:), label1(:,:), del_lab(:,:)
 real(kind=8) :: e_tot = 0d0   ! total electronic energy
 real(kind=8), parameter :: dis_thres0 = 4d0 ! Angstrom
 real(kind=8), allocatable :: cluster_e(:) ! size n_tot
 character(len=240) :: fchname = ' '
 logical :: calc_no = .false. ! calculate density and NOs
 type(mo_cluster), allocatable :: cluster(:), cluster0(:), cluster1(:)

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

!$omp parallel do schedule(dynamic) default(private) shared(n,old_c,new_c)
 do i = 1, n, 1
  call copy_type_mo_cluster(old_c(i), new_c(i))
 end do ! for i
!$omp end parallel do
end subroutine copy_type_mo_clusters

! Generate the MO-cluster for each occupied MO
! 1) doubly occupied MOs and active occupied MOs are considered
! 2) singly occupied MOs are not considered here since they should be included in each MO-cluster
! 3) active unoccupied MOs are not considered here since they are not included in any MO-cluster
subroutine gen_mo_cluster_per_mo(dis_thres)
 implicit none
 integer :: i, j, k, m, nocc
 integer, allocatable :: idx(:)
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
  allocate(idx(nocc))
  call sort_int_array(nocc, cluster0(i)%occ_idx, .true., idx)
  deallocate(idx)
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

! delete primitive MO clusters which are embraced by other primitive MO clusters
! cluster0 -> cluster1
subroutine delete_embraced_cluster0()
 implicit none
 integer :: i, j, n1, n2
 integer, allocatable :: idx(:)
 logical :: alive
 logical, external :: check_subset
 logical, allocatable :: del(:)

 if(n_old == 0) then
  write(6,'(/,A)') 'ERROR in subroutine delete_embraced_cluster0: n_old=0.'
  write(6,'(A)') 'Something must be wrong.'
  stop
 end if
 allocate(del(n_old), source=.false.)

 do i = 1, n_old-1, 1
  if(del(i)) cycle
  n1 = cluster0(i)%nocc
  allocate(idx(n1), source=cluster0(i)%occ_idx)

  do j = i+1, n_old, 1
   if(del(j)) cycle
   n2 = cluster0(j)%nocc

   if(n1 < n2) then
    alive = check_subset(n1, idx, n2, cluster0(j)%occ_idx)
    if(alive) then
     del(i) = .true.
     exit
    end if
   else ! n1 >= n2
    alive = check_subset(n2, cluster0(j)%occ_idx, n1, idx)
    if(alive) del(j) = .true.
   end if
  end do ! for j

  deallocate(idx)
 end do ! for i

 n_new = COUNT(del .eqv. .false.)
 if(n_new == 0) then
  write(6,'(/,A)') 'ERROR in subroutine delete_embraced_cluster: all MO cluster&
                   &s are deleted.'
  write(6,'(A)') 'Something must be wrong.'
  stop
 end if

 allocate(cluster1(n_new))
 j = 0

 do i = 1, n_old, 1
  if(del(i)) cycle
  j = j + 1
  call copy_type_mo_cluster(cluster0(i), cluster1(j))
 end do ! for i

 deallocate(del, cluster0)

 write(6,'(A)') 'merged primitive MO clusters:'
 do i = 1, n_new, 1
  write(6,'(I4,A1,23I4)') i,':',cluster1(i)%occ_idx
 end do ! for i
end subroutine delete_embraced_cluster0

! delete MO clusters which are embraced by other MO clusters
! Note: there may exists empty MO clusters, which will also be deleted
! cluster0 -> cluster1, icoeff0 -> icoeff1
subroutine delete_embraced_cluster()
 implicit none
 integer :: i, j, k, m, n1, n2, n3, n4
 integer, allocatable :: idx(:), idx1(:), idx2(:), itmp(:), del_lab1(:,:)
 logical :: alive, cycle_j
 logical, external :: check_subset
 logical, allocatable :: del(:)

 if(n_old == 0) then
  write(6,'(/,A)') 'ERROR in subroutine delete_embraced_cluster: n_old=0.'
  write(6,'(A)') 'Something must be wrong.'
  stop
 end if
 if(allocated(del_lab)) deallocate(del_lab)
 allocate(del(n_old), source=.false.)
 forall(i=1:n_old, cluster0(i)%nocc==0) del(i) = .true.

 do i = 1, n_old-1, 1
  if(del(i)) cycle
  n1 = cluster0(i)%nocc
  allocate(idx(n1), source=cluster0(i)%occ_idx)

  do j = i+1, n_old, 1
   if(del(j)) cycle
   n2 = cluster0(j)%nocc

   if(n1 < n2) then
    n3 = n1; n4 = n2
    allocate(idx1(n3), source=idx)
    allocate(idx2(n4), source=cluster0(j)%occ_idx)
   else ! n1 >= n2
    n3 = n2; n4 = n1
    allocate(idx1(n3), source=cluster0(j)%occ_idx)
    allocate(idx2(n4), source=idx)
   end if
   alive = check_subset(n3, idx1, n4, idx2)
   deallocate(idx1, idx2)
   if(.not. alive) cycle
   ! now alive is .true., i.e. idx1 is subset of idx2

   k = size(label0, 1)
   allocate(itmp(2*k))
   call find_union(k, label0(:,i), label0(:,j), itmp)
   if(COUNT(itmp>0) > k+1) then ! not the next generation, skip
    deallocate(itmp)
    cycle
   end if
   ! now COUNT(itmp > 0) = k+1

   if(allocated(del_lab)) then
    m = size(del_lab, 2)
    cycle_j = .false.
    do n3 = 1, m, 1
     if(ALL(del_lab(:,n3) == itmp)) then
      cycle_j = .true.
      exit
     end if
    end do ! for n3
    if(cycle_j) then
     deallocate(itmp)
     cycle
    end if

    allocate(del_lab1(k+1,m), source=del_lab)
    deallocate(del_lab)
    allocate(del_lab(k+1,m+1))
    del_lab(:,1:m) = del_lab1
    deallocate(del_lab1)
    del_lab(:,m+1) = itmp
   else ! not allocated
    allocate(del_lab(k+1,1))
    del_lab(:,1) = itmp
   end if

   deallocate(itmp)

   if(n1 < n2) then
    del(i) = .true.
    exit
   else
    del(j) = .true.
   end if
  end do ! for j

  deallocate(idx)
 end do ! for i

 n_new = COUNT(del .eqv. .false.)
 if(n_new == 0) then
  write(6,'(/,A)') 'ERROR in subroutine delete_embraced_cluster: all MO cluster&
                   &s are deleted.'
  write(6,'(A)') 'Something must be wrong.'
  stop
 end if

 k = size(label0, 1)
 allocate(label1(k,n_new), source=0)
 allocate(icoeff1(n_new), source=0)
 allocate(cluster1(n_new))
 j = 0

 do i = 1, n_old, 1
  if(del(i)) cycle
  j = j + 1
  label1(:,j) = label0(:,i)
  icoeff1(j) = icoeff0(i)
  call copy_type_mo_cluster(cluster0(i), cluster1(j))
 end do ! for i

 deallocate(del, label0, icoeff0, cluster0)

 write(6,'(A)') 'merged MO clusters:'
 do i = 1, n_new, 1
  write(6,'(I4,A1,23I4)') icoeff1(i),':',cluster1(i)%occ_idx
 end do ! for i
end subroutine delete_embraced_cluster

! generate primitive MO clusters
subroutine gen_prim_cluster(dis_thres)
 implicit none
 real(kind=8), intent(in) :: dis_thres

 call gen_mo_cluster_per_mo(dis_thres)
 n_old = nmo

 call delete_embraced_cluster0()

 if(n_new == 1) then
  write(6,'(/,A)') 'ERROR in subroutine gen_prim_cluster: only one primitive su&
                   &bsystem is generated'
  write(6,'(A)') 'and it is equal to the whole system. So there is no need to u&
                 &se obf.'
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
 integer, allocatable :: idx(:), itmp(:)
 logical :: cycle_k, alive
 logical, external :: check_subset

 ! copy primitive MO cluster from type cluster to cluster1
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

  allocate(cluster0(n_deri), label0(i,n_deri), icoeff0(n_deri), itmp(2*(i-1)))
  label0 = 0; icoeff0 = 0; m = 0

  do j = 1, n_old-1, 1
   n1 = cluster1(j)%nocc
   allocate(idx(n1), source=cluster1(j)%occ_idx)

   do k = j+1, n_old, 1
    n2 = cluster1(k)%nocc

    ! if the intersection does not belong to this generation, cycle
    ! e.g. (AxB)x(CxD) = AxBxCxD is the 4-th generation, not the 3rd
    call find_union(i-1, label1(:,j), label1(:,k), itmp)
    n = COUNT(itmp > 0)
    if(n > i) then
     cycle
    else if(n < i) then
     write(6,'(/,A)') 'ERROR in subroutine gen_deri_cluster: n<i is impossible.'
     write(6,'(A,3I3,A,10I3)') 'j,k,n=', j,k,n, ', itmp=', itmp
     write(6,'(A,10I3)') 'label1(:,j)=',label1(:,j)
     write(6,'(A,10I3)') 'label1(:,k)=',label1(:,k)
     stop
    end if

    ! if the intersection pattern has appeared before, cycle
    ! e.g. if (AxB)x(CxD) is reserved, (AxC)x(BxD) should be discarded
    cycle_k = .false.
    do n = 1, m, 1
     if(ALL(label0(:,n) == itmp(1:i))) then
      cycle_k = .true.
      exit
     end if
    end do ! for n
    if(cycle_k) cycle

    ! if the intersection pattern has been determined as 'deleted' in the
    ! previous cycle, skip
    if(allocated(del_lab)) then
     cycle_k = .false.
     p = size(del_lab, 2)
     do n = 1, p, 1
      if(ALL(del_lab(:,n) == itmp)) then
       cycle_k = .true.
       exit
      end if
     end do ! for n
     if(cycle_k) cycle
    end if

    m = m + 1
    label0(:,m) = itmp(1:i)
    call find_intersec(n1, idx, n2, cluster1(k)%occ_idx, cluster0(m))
   end do ! for k

   deallocate(idx)
  end do ! for j

  deallocate(itmp, label1, cluster1, icoeff1)
  if(m == 0) then ! no intersection in this generation
   deallocate(label0)
   exit
  end if
  if(MOD(i,2) == 0) then
   forall(j=1:m, cluster0(j)%nocc>0) icoeff0(j) = -1
  else
   forall(j=1:m, cluster0(j)%nocc>0) icoeff0(j) = 1
  end if

  n_old = m   ! will be used in delete_embraced_cluster
  call delete_embraced_cluster()
  ! cluster0 -> cluster1, icoeff0 -> icoeff1, label0 -> label1

  call append_cluster1_to_cluster()

  n_old = n_new   ! update n_old for next cycle
  if(n_old == 1) exit
  !if(ALL(cluster1(:)%nocc < 2)) exit
 end do ! for i

 if(allocated(del_lab)) deallocate(del_lab)
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
 call copy_type_mo_clusters(n_new, cluster1, cluster(n_tot+1:))

 allocate(tmp_i(n_tot), source=icoeff)
 deallocate(icoeff)
 allocate(icoeff(n))
 icoeff(1:n_tot) = tmp_i
 deallocate(tmp_i)
 icoeff(n_tot+1:) = icoeff1

 n_tot = n
end subroutine append_cluster1_to_cluster

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

! merge primitive and derivative MO clusters, n_tot and type cluster will be
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

! Add paired active virtual MO indices into the type cluster.
! This is because we always assume the input fchname includes paired active orb-
! itals, e.g. associated rotated UNO, or GVB orbitals.
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
 i = INDEX(fchname,'.fch', back=.true.)
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
 use population, only: get_mo_dis_from_fch
 use obf, only: calc_no, dis_thres0, n_tot, icoeff, cluster, fchname, &
  gen_prim_cluster, gen_deri_cluster, merge_mo_cluster, add_paired_vir2cluster,&
  gen_permute_fch
 implicit none
 integer :: i, ibegin, iend
 ! ibegin: the beginning index of active occupied MOs
 ! iend: the final index of active occupied MOs, singly occupied not included
 integer(kind=4) :: hostnm
 real(kind=8) :: dis_thres
 character(len=5) :: str = ' '
 character(len=24) :: hostname, data_string

 i = iargc()
 if(i<3 .or. i>4) then
  write(6,'(/,A)') ' ERROR in program obf: wrong command line arguments!'
  write(6,'(A)') ' Example 1: obf tetracene_uno_asrot.fch 52 60'
  write(6,'(A)') ' Example 2: obf tetracene_uno_asrot.fch 52 60 4.0'
  write(6,'(A)') ' Note: do not include singly occupied orbitals since they can&
                 & be recognized'
  write(6,'(A,/)') ' from the .fch file.'
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
 write(6,'(/,A)') 'fchname='//TRIM(fchname)
 write(6,'(2(A,I0),A,F5.2)') 'ibegin=',ibegin,', iend=',iend,', dis_thres=',&
                              dis_thres

 ! perform Mulliken population for each MO, find the centers of each MO,
 ! and calculate the distances between any two MOs.
 call get_mo_dis_from_fch(fchname, ibegin, iend)

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

 write(6,'(A)') 'Debug checkpoint. STOP'
 stop
 ! generate all .fch files with active orbitals permuted near HONO or LUNO
 call gen_permute_fch()

 ! generate all automr input files(.gjf) and submit one by one
 call gen_automr_gjf_and_submit(calc_no)
 deallocate(cluster)

 call read_cluster_e_from_out()

 call fdate(data_string)
 write(6,'(/,A)') 'Obf program ends at '//TRIM(data_string)
end program main

! find the union set of two integer arrays
! Note: all elements in integer arrays a1 and a2 must be positive
subroutine find_union(n, a1, a2, a3)
 implicit none
 integer :: i, j
 integer, intent(in) :: n
 integer, intent(in) :: a1(n), a2(n)
 integer, intent(out) :: a3(2*n)
 integer, allocatable :: idx(:)

 a3(1:n) = a1; a3(n+1:) = 0; j = n

 do i = 1, n, 1
  if(ANY(a1 == a2(i))) cycle
  j = j + 1
  a3(j) = a2(i)
 end do ! for i

 ! sort a3(1:j)
 allocate(idx(j))
 call sort_int_array(j, a3(1:j), .true., idx)
 deallocate(idx)
end subroutine find_union

! find intersection of two integer arrays, save the result into type clus
! Note: assuming all elements in arrays a1 and a2 are positive integers
subroutine find_intersec(n1, a1, n2, a2, clus)
 use population, only: mo_cluster
 implicit none
 integer :: i, j, k, m
 integer, intent(in) :: n1, n2
 integer, intent(in) :: a1(n1), a2(n2)
 integer, allocatable :: occ_idx(:)
 logical, allocatable :: skip1(:), skip2(:)
 type(mo_cluster), intent(out) :: clus

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

! check whether the integer array a1 is a subset of a2
! Note: n1 <= n2 is required
function check_subset(n1, a1, n2, a2) result(subset)
 implicit none
 integer :: i, j, k
 integer, intent(in) :: n1, n2
 integer, intent(in) :: a1(n1), a2(n2)
 logical :: subset

 subset = .true.

 do i = 1, n1, 1
  k = a1(i)

  do j = 1, n2, 1
   if(ALL(a2 /= k)) then
    subset = .false.
    return
   end if
  end do ! for j

 end do ! for i
end function check_subset

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
subroutine gen_automr_gjf_and_submit(calc_no)
 use obf, only: fchname, n_tot, cluster
 implicit none
 integer :: i, ne, fid
 integer, parameter :: nproc = 24, mem = 170
 character(len=240) :: proname, gjfname, new_fch
 logical, intent(in) :: calc_no

 i = INDEX(fchname, '.fch', back=.true.)
 proname = fchname(1:i-1)

 do i = 1, n_tot, 1
  write(gjfname,'(A,I0,A)') TRIM(proname),i,'.gjf'
  write(new_fch,'(A,I0,A)') TRIM(proname),i,'.fch'
  open(newunit=fid,file=TRIM(gjfname),status='replace')

  write(fid,'(A,I0)') '%nprocshared=',nproc
  write(fid,'(A,I0,A)') '%mem=',mem,'GB'
  ne = cluster(i)%nocc
  write(fid,'(2(A,I0),A)') '#p CASCI(',ne,',',ne,')/cc-pVDZ'
  write(fid,'(/,A)',advance='no') "mokit{ist=5,readno='"//TRIM(new_fch)//"'"
  if(.not. calc_no) write(fid,'(A)',advance='no') ',noDMRGNO'
  write(fid,'(A)') '}'
  close(fid)
 end do ! for i

 do i = 1, n_tot, 1
  write(gjfname,'(A,I0,A)') TRIM(proname),i,'.gjf'
  call submit_automr_job(gjfname)
 end do ! for i
end subroutine gen_automr_gjf_and_submit

! 1) read the electronic energy of each MO cluster from output files
! 2) calculate the total density using linear combinations of the total density
! of various MO clusters
subroutine read_cluster_e_from_out()
 use obf, only: fchname, n_tot, e_tot, cluster_e, icoeff
 implicit none
 integer :: i, k, nbf, nif, fid, system
 real(kind=8), allocatable :: den0(:,:), den1(:,:)
 character(len=240) :: buf, proname, proname1, outname, fchname1, fchname2
 character(len=500) :: longbuf = ' '

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(den0(nbf,nbf), source=0d0)
 allocate(den1(nbf,nbf), source=0d0)

 i = INDEX(fchname, '.fch', back=.true.)
 proname = fchname(1:i-1)
 fchname2 = fchname(1:i-1)//'_NO.fch'
 call copy_file(fchname, fchname2, .false.)
 allocate(cluster_e(n_tot), source=0d0)

 do i = 1, n_tot, 1
  write(proname1,'(A,I0)') TRIM(proname),i
  outname = TRIM(proname1)//'_CASCI.out'
  fchname1 = TRIM(proname1)//'_CASCI_NO.fch'

  open(newunit=fid,file=TRIM(outname),status='old',position='append')
  BACKSPACE(fid)
  read(fid,'(A)') buf
  close(fid)

  if(buf(1:7) /= 'CASCI E') then
   write(6,'(/,A)') "ERROR in subroutine gen_automr_gjf_and_submit: 'CASCI E' n&
                    &ot found at the"
   write(6,'(A)') "end of "//TRIM(outname)
   close(fid)
   stop
  end if
  k = INDEX(buf, '=')
  read(buf(k+1:),*) cluster_e(i)

  call read_dm_from_fch(fchname1, 1, nbf, den1)
  den0 = den0 + DBLE(icoeff(i))*den1

  longbuf = 'tar -zcf '//TRIM(proname1)//'.tar.gz '//TRIM(proname1)//'.* '//&
            TRIM(proname1)//'_* --remove-files'
  k = SYSTEM(longbuf)
  if(k /= 0) then
   write(6,'(/,A)') 'ERROR in subroutine read_cluster_e_from_out: failed to cre&
                    &ate package'
   write(6,'(A)') TRIM(proname1)//'.tar.gz'
  end if
 end do ! for i

 deallocate(den1)
 e_tot = DOT_PRODUCT(cluster_e, DBLE(icoeff))
 deallocate(cluster_e)
 write(6,'(/,A,F18.8,A)') 'E_tot = ',e_tot,' a.u.'

 call write_dm_into_fch(fchname2, .true., nbf, den0)
 deallocate(den0)
 call gen_no_using_density_in_fch(fchname2, 1) 
end subroutine read_cluster_e_from_out

