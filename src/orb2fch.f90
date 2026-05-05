! written by jxzou at 20200907: adjust the orders of d,f,g, etc functions in
!  OpenMolcas, to that in Gaussian .fch(k) file
! Originally copied from fch2inporb.f90, some modifications are made
! updated by jxzou at 20210413: remove '-uhf', add automatic determination

program main
 implicit none
 integer :: narg, npair
 character(len=3) :: str = ' '
 character(len=26), parameter :: error_warn = 'ERROR in program orb2fch: '
 character(len=240) :: fchname, orbname
 logical :: gp, prt_no

 narg = iargc()
 if(narg<2 .or. narg>4) then
  write(6,'(/,1X,A)')error_warn//'wrong command line arguments!'
  write(6,'(A)')   ' Example 1 (RHF)   : orb2fch a.ScfOrb a.fch'
  write(6,'(A)')   ' Example 2 (UHF)   : orb2fch a.UhfOrb a.fch'
  write(6,'(A)')   ' Example 3 (CAS)   : orb2fch a.RasOrb a.fch'
  write(6,'(A)')   ' Example 4 (UNO)   : orb2fch a.UnaOrb a.fch -no'
  write(6,'(A)')   ' Example 5 (CAS NO): orb2fch a.RasOrb.1 a.fch -no'
  write(6,'(A,/)') ' Example 6 (GP NO) : orb2fch ben.RasOrb.1 ben.fch -gp 3'
  stop
 end if

 npair = 0; gp = .false.; prt_no = .false.; fchname = ' '; orbname = ' '
 call getarg(1, orbname)
 call require_file_exist(orbname)

 call getarg(2, fchname)
 call require_file_exist(fchname)

 if(narg > 2) then
  call getarg(3, str)
  select case(TRIM(str))
  case('-no') ! copy natural orbital occupation numbers into .fch
   prt_no = .true.
   if(narg > 3) then
    write(6,'(/,A)') error_warn//'no more argument is accepted when'
    write(6,'(A)') '"-no" is specified.'
    stop
   end if
  case('-gp') ! transform GP NOs from .xxxOrb to .fch
   gp = .true.; prt_no = .true.
   if(narg == 4) then
    call getarg(4, str)
    read(str,*) npair ! < 1000 pairs due to length of str
    if(npair < 0) then
     write(6,'(/,A)') error_warn//'npair>=0 is required!'
     write(6,'(A,I0)') 'But got npair=', npair
     stop
    else if(npair == 0) then
     gp = .false. ! degrade to a R(O)HF calculation
    end if
   else ! narg /= 4
    write(6,'(/,A)') error_warn//'when "-gp" is specified, the number of'
    write(6,'(A)') 'pairs must also be specified. For example,'
    write(6,'(A)') ' orb2fch ben.RasOrb.1 ben.fch -gp 3'
    write(6,'(A)') ' orb2fch xxx.RasOrb.1 xxx.fch -gp 0'
    stop
   end if
  case default
   write(6,'(/,A)') error_warn//'the 3rd argument is wrong! Currently'
   write(6,'(A)') 'only "-no" or "-gp" is accepted. But got "'//str//'"'
   stop
  end select
 end if

 call orb2fch(orbname, fchname, prt_no, gp, npair)
end program main

! read the MOs in orbital file of OpenMolcas and adjust its d,f,g,h functions
!  order to that of Gaussian
subroutine orb2fch(orbname, fchname, prt_no, gp, npair)
 implicit none
 integer :: i, j, k, m, length
 integer :: na, nb, nbf, nif, nbf0, nbf1
 integer :: n6dmark, n10fmark, n15gmark, n21hmark, n28imark
 integer :: n5dmark, n7fmark, n9gmark, n11hmark, n13imark
 integer, allocatable :: shell_type(:), shell2atom_map(:)
 integer, allocatable :: idx(:), idx2(:), d_mark(:), f_mark(:), g_mark(:), &
  h_mark(:), i_mark(:)
 integer, intent(in) :: npair
 real(kind=8) :: rtmp
 real(kind=8), allocatable :: coeff(:,:), coeff2(:,:), occ_num(:), norm(:)
 character(len=29), parameter :: error_warn = 'ERROR in subroutine orb2fch: '
 character(len=240), intent(in) :: orbname, fchname
 ! orbname is one of .ScfOrb/.RasOrb/.RasOrb.1/.UnaOrb/.UhfOrb file of OpenMolcas
 logical :: uhf, sph, alter
 logical, intent(in) :: prt_no, gp

 call check_uhf_in_fch(fchname, uhf)
 if(uhf .and. gp) then
  write(6,'(/,A)') error_warn//'UHF and GP are both activated.'
  write(6,'(A)') 'This is not allowed since Unrestricted GP is not supported.'
  write(6,'(A)') 'orbname='//TRIM(orbname)
  write(6,'(A)') 'fchname='//TRIM(fchname)
  stop
 end if
 ! TODO: gp=.true. is not tested with mult > 1

 call read_na_and_nb_from_fch(fchname, na, nb)
 if(gp .and. npair>nb) then
  write(6,'(/,A)') error_warn//'the number of pairs is too large.'
  write(6,'(2(A,I0))') 'nb=', nb, ', npair=', npair
  write(6,'(A)') 'fchname='//TRIM(fchname)
  stop
 end if

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 nbf0 = nbf ! make a copy of nbf

 ! read MO Coefficients from .ScfOrb, .RasOrb, .RasOrb.1, .UnaOrb, or .UhfOrb
 ! file of OpenMolcas
 if(uhf) then
  allocate(coeff(nbf,2*nif))
  call read_mo_from_orb(orbname, nbf, nif, 'a', coeff(:,1:nif))
  call read_mo_from_orb(orbname, nbf, nif, 'b', coeff(:,nif+1:2*nif))
  nif = 2*nif   ! double the size
 else
  if(prt_no) then
   allocate(occ_num(nif), coeff(nbf,nif))
   call read_on_from_orb(orbname, nif, 'a', occ_num)
   call read_mo_from_orb(orbname, nbf, nif, 'a', coeff)
  else
   allocate(coeff(nbf,nif))
   call read_mo_from_orb(orbname, nbf, nif, 'a', coeff)
  end if
 end if

 call read_ncontr_from_fch(fchname, k)
 allocate(shell_type(2*k), source=0)
 allocate(shell2atom_map(2*k), source=0)
 call read_shltyp_and_shl2atm_from_fch(fchname, k, shell_type, shell2atom_map)

 ! check if any spherical functions
 if(ANY(shell_type<-1) .and. ANY(shell_type>1)) then
  write(6,'(/,A)') error_warn//'mixed spherical harmonic/Cartesian functions'
  write(6,'(A)') 'detected. You probably used a basis set like 6-31G(d) in Gaus&
                 &sian. Its default'
  write(6,'(A)') 'setting is (6D,7F). You need to add `5D 7F` (recommended) or &
                 &`6D 10F` keywords'
  write(6,'(A)') 'into gjf.'
  stop
 else if( ANY(shell_type>1) ) then
  sph = .false.
 else
  sph = .true.
 end if

! first we adjust the basis functions in each MO according to the Shell to atom map
 ! 1) split the 'L' into 'S' and 'P', this is to ensure that D comes after L functions
 call split_L_func(k, shell_type, shell2atom_map, length)
 allocate(idx(nbf))
 forall(i = 1:nbf) idx(i) = i
 allocate(norm(nbf), source=1d0)

 ! 2) sort the shell_type and shell2atom_map by ascending order,
 !  the indices of MOs will be adjusted accordingly
 call sort_shell_and_mo_idx(length, shell_type, shell2atom_map, nbf, idx)
! adjust done

! then we adjust the basis functions in each MO according to the type of basis functions
 k = length  ! update k
 allocate(d_mark(k), f_mark(k), g_mark(k), h_mark(k), i_mark(k))

 ! adjust the order of d, f, etc. functions
 if(sph) then ! spherical harmonic
  call read_mark_from_shltyp_sph(k, shell_type, n5dmark, n7fmark, n9gmark, &
                 n11hmark, n13imark, d_mark, f_mark, g_mark, h_mark, i_mark)
  call fch2inporb_permute_sph(n5dmark, n7fmark, n9gmark, n11hmark, k, d_mark, &
                              f_mark, g_mark, h_mark, nbf, idx)
 else ! Cartesian-type basis
  call read_mark_from_shltyp_cart(k, shell_type, n6dmark, n10fmark, n15gmark, &
                    n21hmark, n28imark, d_mark, f_mark, g_mark, h_mark, i_mark)
  call fch2inporb_permute_cart(n6dmark, n10fmark, n15gmark, n21hmark, k, &
                               d_mark, f_mark, g_mark, h_mark, nbf, idx, norm)
 end if
! adjustment finished
 deallocate(d_mark, f_mark, g_mark, h_mark, i_mark)

! move the 2nd, 3rd, ... Zeta basis functions forward
 i = 0
 nbf = 0
 do while(i < k)
  i = i + 1
  j = shell2atom_map(i)
  m = shell_type(i)
  nbf1 = nbf
  select case(m)
  case( 0)   ! S
   nbf = nbf + 1
  case( 1)   ! 3P
   nbf = nbf + 3
  case(-2)   ! 5D
   nbf = nbf + 5
  case( 2)   ! 6D
   nbf = nbf + 6
  case(-3)   ! 7F
   nbf = nbf + 7
  case( 3)   ! 10F
   nbf = nbf + 10
  case(-4)   ! 9G
   nbf = nbf + 9
  case( 4)   ! 15G
   nbf = nbf + 15
  case(-5)   ! 11H
   nbf = nbf + 11
  case( 5)   ! 21H
   nbf = nbf + 21
  end select
  if(m == 0) cycle

  length = 1
  do while(i < k)
   i = i + 1
   if(shell_type(i) /= m) exit
   if(shell2atom_map(i) /= j) exit
   length = length + 1
  end do ! for while

  if(i < k) i = i - 1
  if(length > 1) then
   call zeta_mv_forwd_idx(nbf1, m, length, nbf0, idx, norm)
   nbf = nbf1 + length*(nbf-nbf1)
  end if
 end do ! for while

 deallocate(shell_type, shell2atom_map)
! move done

 nbf = nbf0
 allocate(idx2(nbf), coeff2(nbf,nif))
 forall(i = 1:nbf) idx2(idx(i)) = i
 deallocate(idx)
 if(sph) then
  forall(i=1:nbf, j=1:nif) coeff2(i,j) = coeff(idx2(i),j)
 else
  forall(i=1:nbf, j=1:nif) coeff2(i,j) = coeff(idx2(i),j)/norm(idx2(i))
 end if
 deallocate(idx2, norm, coeff)

 ! Permute GP NOs and NOONs from GAMESS order to Gaussian order.
 ! GAMESS order: bonding1, antibonding1, bonding2, antibonding2, ...
 ! Gaussian order: bonding1, bonding2, antibonding2, antibonding1, ...
 ! And make N(bonding1) > N(bonding2), N(antibonding1) < N(antibonding2)
 if(gp .and. npair>0) then
  k = 2*npair
  allocate(idx(k))
  do i = 1, npair, 1
   idx(i) = 2*i - 1
   idx(npair+i) = k - 2*i + 2
  end do ! for i
  m = nb - npair
  allocate(norm(k))
  do i = 1, k, 1
   norm(i) = occ_num(m+idx(i))
  end do ! for i
  do while(.true.)
   alter = .false.
   do i = 1, npair-1, 1
    if(norm(i) < norm(i+1)) then
     rtmp=norm(i); norm(i)=norm(i+1); norm(i+1)=rtmp
     rtmp=norm(k-i+1); norm(k-i+1)=norm(k-i); norm(k-i)=rtmp
     j=idx(i); idx(i)=idx(i+1); idx(i+1)=j
     j=idx(k-i+1); idx(k-i+1)=idx(k-i); idx(k-i)=j
     alter = .true.
    end if
   end do ! for i
   if(.not. alter) exit
  end do ! for while
  occ_num(m+1:m+k) = norm
  deallocate(norm)
  allocate(coeff(nbf,k), source=coeff2(:,m+1:m+k))
  do i = 1, k, 1
   coeff2(:,m+i) = coeff(:,idx(i))
  end do ! for i
  deallocate(coeff)
 end if

 if(uhf) then ! UHF-type
  nif = nif/2
  call write_mo_into_fch(fchname, nbf, nif, 'a', coeff2(:,1:nif))
  call write_mo_into_fch(fchname, nbf, nif, 'b', coeff2(:,nif+1:2*nif))
 else         ! not UHF-type
  if(prt_no) then
   call write_eigenvalues_to_fch(fchname, nif, 'a', occ_num, .true.)
   call write_mo_into_fch(fchname, nbf, nif, 'a', coeff2)
  else
   call write_mo_into_fch(fchname, nbf, nif, 'a', coeff2)
  end if
 end if

 deallocate(coeff2)
 if(allocated(occ_num)) deallocate(occ_num)
end subroutine orb2fch

