! written by jxzou at 20230519: transfer MOs from Gaussian -> AMESP

! Note: You should use AMESP >= 2023 May 28. Older AMESP versions have bugs in
!  ECP/PP

! Current limitations:
! 1) not supported for GHF
! 2) not supported for different basis set for the same element

program main
 use util_wrapper, only: formchk
 implicit none
 integer :: i
 character(len=240) :: fchname

 i = iargc()
 if(i /= 1) then
  write(6,'(/,A)') ' ERROR in subroutine fch2amo: wrong command line arguments!'
  write(6,'(A,/)') ' Example: fch2amo water.fch'
  stop
 end if

 fchname = ' '
 call getarg(1, fchname)
 call require_file_exist(fchname)

 ! if .chk file provided, convert into .fch file automatically
 i = LEN_TRIM(fchname)
 if(fchname(i-3:i) == '.chk') then
  call formchk(fchname)
  fchname = fchname(1:i-3)//'fch'
 end if

 call fch2amo(fchname)
end program main

subroutine fch2amo(fchname)
 use fch_content
 implicit none
 integer :: i, j, k, m, n, n1, n2, p, q, nelmtyp, length, fid
 integer, parameter :: shltyp2nbas(-5:5) = [21,15,10,6,4,1,3,6,10,15,21]
 integer, allocatable :: frozen_e(:) ! size natom, frozen core electrons
 integer, allocatable :: natmbas(:) ! the number of basis functions of each atom
 integer, allocatable :: atmnshl(:) ! the number of shells of each atom
 integer, allocatable :: end_idx(:), itmp(:), jtmp(:,:), ktmp(:,:)
 real(kind=8), allocatable :: bmnl(:), ecpexp(:,:), ecpcof(:,:)
 character(len=1) :: str = ' '
 character(len=2) :: str2 = '  '
 character(len=1), parameter :: am_type(-1:6) = ['L','S','P','D','F','G','H','I']
 character(len=1), parameter :: am_type1(0:6) = ['s','p','d','f','g','h','i']
 character(len=3) :: sph_str
 character(len=240) :: inpname, amoname
 character(len=240), intent(in) :: fchname
 logical :: uhf, sph, has_sp, ecp
 logical, allocatable :: skip_elem(:) ! .True. for skipping this element

 i = index(fchname, '.fch', back=.true.)
 if(i == 0) then
  write(6,'(/,A)') "ERROR in subroutine fch2amo: '.fch' suffix not found in fil&
                   &ename "//TRIM(fchname)
  stop
 end if

 inpname = fchname(1:i-1)//'.aip'
 amoname = fchname(1:i-1)//'.amo'
 call check_nobasistransform_in_fch(fchname)
 call check_nosymm_in_fch(fchname)

 uhf = .false.; sph = .false.; has_sp = .false.; ecp = .false.
 call check_uhf_in_fch(fchname, uhf) ! determine whether UHF
 call read_fch(fchname, uhf)
 if(LenNCZ > 0) ecp = .true.

 ! check if any spherical functions
 if(ANY(shell_type<-1) .and. ANY(shell_type>1)) then
  write(6,'(A)') 'ERROR in subroutine fch2qchem: mixed spherical harmonic/&
                 &Cartesian functions detected.'
  write(6,'(A)') 'You probably used a basis set like 6-31G(d) in Gaussian. Its&
                & default setting is (6D,7F).'
  write(6,'(A)') "You need to add '5D 7F' or '6D 10F' keywords in Gaussian."
  stop
 else if( ANY(shell_type < -1) ) then
  sph = .true.
 end if

! Firstly, generate the input file (.in)
 open(newunit=fid,file=TRIM(inpname),status='replace')
 write(fid,'(A)') '% npara 4'
 write(fid,'(A)') '% maxcore 1000'

 write(fid,'(A)',advance='no') '! '
 if(uhf) then
  write(fid,'(A)',advance='no') 'uhf'
 else ! RHF, ROHF
  if(mult == 1) then
   write(fid,'(A)',advance='no') 'hf'
  else
   write(fid,'(A)',advance='no') 'rohf'
  end if
 end if
 write(fid,'(A)') ' define'

 if(ANY(shell_type>1)) write(fid,'(A,/,A,/,A)') '>ope',' inttype car','end'
 write(fid,'(A)') '>scf'
 write(fid,'(A)') ' guess read'
 write(fid,'(A)') ' scfmode direct'
 write(fid,'(A)') 'end'
 write(fid,'(A,I0,1X,I0)') '>xyz ', charge, mult

 do i = 1, natom, 1
  write(fid,'(A2,1X,3(1X,F18.8))') elem(i), coor(:,i)
 end do ! for i

 write(fid,'(A,/,A)') 'end','>basis'
 allocate(skip_elem(natom))
 skip_elem = .false.
 nelmtyp = 1

 do i = 2, natom, 1
  if( ANY(elem(1:i-1) == elem(i)) ) then
   skip_elem(i) = .true.
  else
   nelmtyp = nelmtyp + 1
  end if
 end do ! for i

 do i = 1, natom, 1
  if(skip_elem(i)) cycle

  if(allocated(LPSkip)) then
   if(LPSkip(i) == 0) then
    write(fid,'(1X,A)') TRIM(elem(i))//' readecp'
   else
    write(fid,'(1X,A)') TRIM(elem(i))//' read'
   end if
  else
   write(fid,'(1X,A)') TRIM(elem(i))//' read'
  end if
 end do ! for i
 write(fid,'(A,//,A)') 'end','$basis'

 ! print basis sets into the .in file
 write(fid,'(A,/,A2,/,A)') ' *',elem(1), ' *'
 k = 0
 do i = 1, ncontr, 1
  m = shell2atom_map(i)

  if(skip_elem(m)) then
   k = k + prim_per_shell(i) ! remember to update k
   cycle
  end if

  if(m > 1) then
   if(shell2atom_map(i-1) == m-1) then
    write(fid,'(A,/,A2,/,A)') ' *', elem(m), ' *'
   end if
  end if

  m = shell_type(i); n = prim_per_shell(i)
  if(m < -1) m = -m
  if(m == -1) then
   str2 = 'S'
  else
   str2 = am_type(m)//' '
  end if
  write(fid,'(I2,1X,A)') n, str2

  has_sp = .false.
  if(allocated(contr_coeff_sp)) then
   if(ANY(contr_coeff_sp(k+1:k+n) > 1d-6)) has_sp = .true.
  end if

  if(has_sp) then
   do j = k+1, k+n, 1
    write(fid,'(3(2X,ES15.8))') prim_exp(j), contr_coeff(j)
   end do ! for j
   write(fid,'(I2,1X,A)') n, 'P'
   do j = k+1, k+n, 1
    write(fid,'(3(2X,ES15.8))') prim_exp(j), contr_coeff_sp(j)
   end do ! for j
  else ! no SP in this paragraph
   do j = k+1, k+n, 1
    write(fid,'(2(2X,ES15.8))') prim_exp(j), contr_coeff(j)
   end do ! for j
  end if

  k = k + n
 end do ! for i

 write(fid,'(A,/,A)') ' *','$end'

 if(ecp) then
  allocate(frozen_e(natom))
  forall(i = 1:natom) frozen_e(i) = INT(RNFroz(i))
  deallocate(RNFroz)
  write(fid,'(/,A)') '$ecp'

  do i = 1, natom, 1
   if(LPSkip(i) /= 0) cycle
   write(fid,'(A,/,A2,/,A)') ' *',elem(i),' *'
   write(fid,'(2(A,I0))') 'ncore = ',frozen_e(i),'  lmax = ',LMax(i)
   str = am_type1(LMax(i))

   do j = 1, 10, 1
    n1 = KFirst(i,j); n2 = KLast(i,j)
    if(n1 == 0) exit
    if(j == 1) then
     write(fid,'(A)',advance='no') str//'-ul'
    else
     write(fid,'(A)',advance='no') am_type1(j-2)//'-ul'
    end if
    write(fid,'(1X,I2)') n2-n1+1
    do n = n1, n2, 1
     write(fid,'(ES15.8,3X,I0,3X,ES15.8)') CLP(n), NLP(n), ZLP(n)
    end do ! for n
   end do ! for j

   write(fid,'(A)') ' *'
  end do ! for i

  write(fid,'(A)') '$end'
 end if

 close(fid)

! Secondly, generate the .amo orbital file
 open(newunit=fid,file=TRIM(amoname),status='replace')
 write(fid,'(A)') "[Amesp's MO File Generated by fch2amo]"
 write(fid,'(/,A,/,A,I0)') '[Atoms]', 'Natm= ', natom
 coor = coor/Bohr_const

 if(ecp) then
  do i = 1, natom, 1
   write(fid,'(A2,2X,I3,3(1X,F18.8))') elem(i), ielem(i)-frozen_e(i), coor(:,i)
  end do ! for i
  deallocate(frozen_e)
 else
  do i = 1, natom, 1
   write(fid,'(A2,2X,I3,3(1X,F18.8))') elem(i), ielem(i), coor(:,i)
  end do ! for i
 end if

 deallocate(coor)
 write(fid,'(A,/,A,I0)') '[GTO]', 'NAtom_Type: ', nelmtyp

 ! Check whether 'SP' or 'L' exists, Pople type basis set
 ! If yes, enlarge arrays prim_per_shell, shell_type, shell2atom_map
 if(ANY(shell_type == -1)) then
  allocate(itmp(ncontr), source=prim_per_shell)
  deallocate(prim_per_shell)
  allocate(prim_per_shell(2*ncontr), source=0)
  prim_per_shell(1:ncontr) = itmp
  call split_prim_per_shell(ncontr, shell_type, prim_per_shell, k)

  k = SUM(prim_per_shell)
  allocate(bmnl(nprim), source=prim_exp)
  deallocate(prim_exp)
  allocate(prim_exp(k), source=0d0)
  prim_exp(1:nprim) = bmnl

  bmnl = contr_coeff
  deallocate(contr_coeff)
  allocate(contr_coeff(k), source=0d0)
  contr_coeff(1:nprim) = bmnl
  deallocate(bmnl)
  call split_exp_and_contr_coeff(ncontr, shell_type, itmp, nprim, contr_coeff_sp,&
                                 k, prim_exp, contr_coeff)
  nprim = k ! update nprim

  itmp = shell_type
  deallocate(shell_type)
  allocate(shell_type(2*ncontr), source=0)
  shell_type(1:ncontr) = itmp

  itmp = shell2atom_map
  deallocate(shell2atom_map)
  allocate(shell2atom_map(2*ncontr), source=0)
  shell2atom_map(1:ncontr) = itmp
  deallocate(itmp)

  call split_L_func(ncontr, shell_type, shell2atom_map, length)
  ncontr = length ! update ncontr
 end if

 allocate(atmnshl(natom), source=0)
 do i = 1, natom, 1
  atmnshl(i) = COUNT(shell2atom_map == i)
 end do ! for i
 write(fid,'(A,I0)') 'maxNShl: ', MAXVAL(atmnshl)
 write(fid,'(A,I0)') 'maxcont: ', MAXVAL(prim_per_shell)

 allocate(natmbas(natom), source=0)
 do i = 1, ncontr, 1
  j = shell2atom_map(i)
  natmbas(j) = natmbas(j) + shltyp2nbas(shell_type(i))
 end do ! for i
 write(fid,'(A,I0)') 'MaxAtmBas: ', MAXVAL(natmbas)

 ! find the end indices of each atom in the array shell2atom_map
 allocate(end_idx(0:natom), source=0)
 do i = 1, ncontr, 1
  end_idx(shell2atom_map(i)) = i
 end do ! for i
 deallocate(shell2atom_map)

 k = 0
 do i = 1, natom, 1
  m = end_idx(i) - end_idx(i-1)
  allocate(itmp(m), source=prim_per_shell(end_idx(i-1)+1:end_idx(i)))

  if(skip_elem(i)) then ! update k then cycle
   k = k + SUM(itmp)
   deallocate(itmp)
   cycle
  end if

  write(fid,'(A2,2X,I0)') elem(i), ielem(i)
  write(fid,'(A,I0)') ' AtmNShl: ', atmnshl(i)
  write(fid,'(A)') ' AngShl:'
  allocate(jtmp(m,1))
  jtmp(:,1) = IABS(shell_type(end_idx(i-1)+1:end_idx(i)))
  write(fid,'(8I8)') jtmp(:,1)
  write(fid,'(A)') ' bas_num:'
  write(fid,'(8I8)') itmp

  write(fid,'(A)') ' a_Basis:'
  n = k ! copy
  do j = 1, m, 1
   write(fid,'(5(1X,ES15.8))') prim_exp(n+1:n+itmp(j))
   n = n + itmp(j)
  end do ! for j

  write(fid,'(A)') ' BasCa:'
  do j = 1, m, 1
   p = itmp(j)
   call normalize_contr_coeff(jtmp(j,1),p,prim_exp(k+1:k+p),contr_coeff(k+1:k+p))
   write(fid,'(5(1X,ES15.8))') contr_coeff(k+1:k+p)
   k = k + p
  end do ! for j
  deallocate(itmp)

  j = natmbas(i)
  write(fid,'(A,I0)') ' Bmnl: ', j
  allocate(bmnl(j))
  call gen_bmnl_from_shltyp(m, jtmp(:,1), j, bmnl)
  write(fid,'(5(1X,ES15.8))') bmnl
  deallocate(jtmp, bmnl)
 end do ! for i

 deallocate(elem, ielem, atmnshl, shell_type, end_idx, prim_exp, contr_coeff, &
            natmbas, skip_elem)

 if(ecp) then
  write(fid,'(A)') 'ECP:   T'
  k = COUNT(LPSkip == 0)
  write(fid,'(A,/,I8)') 'NECPAtm:', k
  write(fid,'(A)') 'lECPAtm:'
  write(fid,'(8L8)') ((LPSkip(i)==0),i=1,natom)

  allocate(itmp(k),source=0)
  j = 0
  do i = 1, natom, 1
   if(LPSkip(i) /= 0) cycle
   j = j + 1
   itmp(j) = LMax(i) + 1
  end do ! for i
  deallocate(LMax)
  write(fid,'(A)') 'ECPBlock:'
  write(fid,'(8I8)') itmp
  deallocate(itmp)

  p = 6*k
  allocate(jtmp(6,k), ktmp(10,p), ecpexp(10,p), ecpcof(10,p))
  jtmp = 0; ktmp = 0; ecpexp = 0d0; ecpcof = 0d0; n = 0; p = 0

  do i = 1, natom, 1
   if(LPSkip(i) /= 0) cycle
   n = n + 1
   m = 0
   do j = 1, 10, 1
    n1 = KFirst(i,j); n2 = KLast(i,j)
    if(n1 == 0) exit
    q = n2-n1+1; m = m+1; p = p+1
    jtmp(m,n) = q
    ktmp(1:q,p) = NLP(n1:n2)
    ecpexp(1:q,p) = ZLP(n1:n2)
    ecpcof(1:q,p) = CLP(n1:n2)
   end do ! for j
   if(MOD(p,6) > 0) p = (p/6 + 1)*6 ! make p as 6*N
  end do ! for i

  deallocate(NLP, ZLP, CLP, CLP2, KFirst, KLast, LPSkip)
  write(fid,'(A)') 'ECPNum:'
  write(fid,'(8I8)') ((jtmp(j,i),j=1,6),i=1,k)
  write(fid,'(A)') 'ECPPow:'
  p = 6*k
  write(fid,'(8I8)') ((ktmp(j,i),j=1,10),i=1,p)
  write(fid,'(A)') 'ECPExp:'
  write(fid,'(5(1X,ES15.8))') ((ecpexp(j,i),j=1,10),i=1,p)
  write(fid,'(A)') 'ECPCof:'
  write(fid,'(5(1X,ES15.8))') ((ecpcof(j,i),j=1,10),i=1,p)
  deallocate(jtmp, ktmp, ecpexp, ecpcof)
 else
  write(fid,'(A)') 'ECP:   F'
 end if

 sph_str = 'car'
 if(sph) sph_str = 'sph'

 if((.not.uhf) .and. mult==1) then
  write(fid,'(A)') '[MO] '//sph_str//' C noRI'
  write(fid,'(A,I0)') 'Nocc= ', na
  write(fid,'(A)') 'En:'
  write(fid,'(5(1X,ES15.8))') eigen_e_a
  write(fid,'(A)') 'MoCu:'
  write(fid,'(5(1X,ES15.8))') alpha_coeff
 else
  write(fid,'(A)') '[MO] '//sph_str//' O noRI'
  write(fid,'(A,I0)') 'NoccA= ', na
  write(fid,'(A)') 'EnA:'
  write(fid,'(5(1X,ES15.8))') eigen_e_a
  write(fid,'(A)') 'MoCuA:'
  write(fid,'(5(1X,ES15.8))') alpha_coeff
  write(fid,'(A,I0)') 'NoccB= ', nb
  write(fid,'(A)') 'EnB:'
  if(uhf) then
   write(fid,'(5(1X,ES15.8))') eigen_e_b
  else
   write(fid,'(5(1X,ES15.8))') eigen_e_a
  end if
  write(fid,'(A)') 'MoCuB:'
  if(uhf) then
   write(fid,'(5(1X,ES15.8))') beta_coeff
  else
   write(fid,'(5(1X,ES15.8))') alpha_coeff
  end if
 end if

 close(fid)
 deallocate(eigen_e_a, alpha_coeff)
 if(uhf) deallocate(eigen_e_b, beta_coeff)
end subroutine fch2amo

! compute the array Bnml using the shell type of an atom
subroutine gen_bmnl_from_shltyp(n, shltyp, nb, bmnl)
 implicit none
 integer :: i, k
 integer, intent(in) :: n, nb
 integer, intent(in) :: shltyp(n)
 real(kind=8), intent(out) :: bmnl(nb)
 real(kind=8), parameter :: r1 = 1d0/DSQRT(3d0)
 real(kind=8), parameter :: r2 = 1d0/DSQRT(15d0)
 real(kind=8), parameter :: r3 = 1d0/DSQRT(105d0)
 real(kind=8), parameter :: r4 = 1d0/3d0
 real(kind=8), parameter :: r5 = 1d0/DSQRT(945d0)
 real(kind=8), parameter :: r6 = 1d0/DSQRT(45d0)
 real(kind=8), parameter :: array1(15) = [r3,r2,r4,r2,r3,r2,r1,r1,r2,r4,r1,r4,&
                                          r2,r2,r3]
 real(kind=8), parameter :: array2(21) = [r5,r3,r6,r6,r3,r5,r3,r2,r4,r2,r3,r6,&
                                          r4,r4,r6,r6,r2,r6,r3,r3,r5]

 bmnl = 1d0 ! most elements are 1.0, just calculate special elements
 k = 0

 do i = 1, n, 1
  select case(shltyp(i))
  case(0) ! S
   k = k + 1
  case(1) ! P
   k = k + 3
  case(2) ! D
   bmnl(k+1:k+3) = r1
   k = k + 6
  case(3) ! F
   bmnl(k+1:k+3) = r2
   bmnl(k+4:k+9) = r1
   k = k + 10
  case(4) ! G
   bmnl(k+1:k+15) = array1
   k = k + 15
  case(5) ! H
   bmnl(k+1:k+21) = array2
   k = k + 21
  case default
   write(6,'(/,A)') 'ERROR in subroutine gen_bmnl_from_shltyp: shell type out o&
                    &f range.'
   write(6,'(2(A,I0))') 'i=', i, ', shltyp(i)=', shltyp(i)
   stop
  end select
 end do ! for i
end subroutine gen_bmnl_from_shltyp

subroutine normalize_contr_coeff(ang, p, prim_exp, contr_coeff)
 implicit none
 integer, intent(in) :: ang, p
 real(kind=8), intent(in) :: prim_exp(p)
 real(kind=8), intent(out) :: contr_coeff(p)
 real(kind=8), parameter :: PI = 4d0*DATAN(1d0)
 real(kind=8), parameter :: tp = 8d0*DATAN(1d0) ! 2*PI
 real(kind=8), parameter :: stp = 2d0*DSQRT(2d0*DATAN(1d0))
 real(kind=8), parameter :: tps = 64d0*DATAN(1d0)*DATAN(1d0)

 select case(ang)
 case(0)
  contr_coeff = contr_coeff*(2d0*prim_exp/PI)**0.75d0
 case(1)
  contr_coeff = contr_coeff*stp*(2d0*prim_exp/PI)**1.25d0
 case(2)
  contr_coeff = contr_coeff*tp*(2d0*prim_exp/PI)**1.75d0
 case(3)
  contr_coeff = contr_coeff*stp*tp*(2d0*prim_exp/PI)**2.25d0
 case(4)
  contr_coeff = contr_coeff*tps*(2d0*prim_exp/PI)**2.75d0
 case(5)
  contr_coeff = contr_coeff*stp*tps*(2d0*prim_exp/PI)**3.25d0
 case default
  write(6,'(A,I0)') 'ERROR in subroutine normalize_contr_coeff: invalid ang=',ang
  stop
 end select
end subroutine normalize_contr_coeff

subroutine split_prim_per_shell(ncontr, shell_type, prim_per_shell, length)
 implicit none
 integer :: i, j, k
 integer, intent(in) :: ncontr
 integer, intent(in) :: shell_type(ncontr)
 integer, intent(inout) :: length
 integer, intent(inout) :: prim_per_shell(2*ncontr)

 j = 1; k = ncontr

 do i = 1, ncontr, 1
  if(shell_type(i) == -1) then
   prim_per_shell(j+2:k+1) = prim_per_shell(j+1:k)
   prim_per_shell(j+1) = prim_per_shell(j)
   j = j + 2
   k = k + 1
  else
   j = j + 1
  end if
 end do ! for i

 length = k
end subroutine split_prim_per_shell

subroutine split_exp_and_contr_coeff(ncontr, shell_type, prim_per_shell, nprim,&
                                     contr_coeff_sp, k1, prim_exp, contr_coeff)
 implicit none
 integer :: i, j, k, m, n
 integer, intent(in) :: ncontr, nprim, k1
 integer, intent(in) :: shell_type(ncontr), prim_per_shell(ncontr)
 real(kind=8), intent(in) :: contr_coeff_sp(nprim)
 real(kind=8), intent(inout) :: prim_exp(k1), contr_coeff(k1)

 j = 1; k = nprim; n = 1

 do i = 1, ncontr, 1
  m = prim_per_shell(i)

  if(shell_type(i) == -1) then
   prim_exp(j+2*m:k+m) = prim_exp(j+m:k)
   prim_exp(j+m:j+2*m-1) = prim_exp(j:j+m-1)
   contr_coeff(j+2*m:k+m) = contr_coeff(j+m:k)
   contr_coeff(j+m:j+2*m-1) = contr_coeff_sp(n:n+m-1)
   j = j + 2*m
   k = k + m
  else
   j = j + m
  end if

  n = n + m
 end do ! for i
end subroutine split_exp_and_contr_coeff

