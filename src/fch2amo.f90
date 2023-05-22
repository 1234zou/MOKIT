! written by jxzou at 20230519: transfer MOs from Gaussian -> AMESP

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
 integer :: i, j, k, m, n, n1, n2, p, nelmtyp, fid
 integer, parameter :: shltyp2nbas(-5:5) = [21,15,10,6,4,1,3,6,10,15,21]
 integer, allocatable :: frozen_e(:) ! size natom, frozen core electrons
 integer, allocatable :: natmbas(:) ! the number of basis functions of each atom
 integer, allocatable :: atmnshl(:) ! the number of shells of each atom
 integer, allocatable :: end_idx(:), itmp(:), jtmp(:,:), ktmp(:,:)
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
                   &e "//TRIM(fchname)
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

 if(.not. sph) write(fid,'(A,/,A,/,A)') '>ope',' inttype car','end'
 write(fid,'(A,/,A,/,A)') '>scf',' guess read','end'
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
   if(LPSkip(i)==0) write(fid,'(1X,A)') TRIM(elem(i))//' readecp'
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
   str2 = 'SP'
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
    write(fid,'(3(2X,ES15.8))') prim_exp(j), contr_coeff(j), contr_coeff_sp(j)
   end do ! for j
  else ! no SP in this paragraph
   do j = k+1, k+n, 1
    write(fid,'(2(2X,ES15.8))') prim_exp(j), contr_coeff(j)
   end do ! for j
  end if

  k = k + n
 end do ! for i

 write(fid,'(A,/,A)') ' *','$end'
 if(allocated(contr_coeff_sp)) deallocate(contr_coeff_sp)

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

! Secondly, permute MO coefficients and generate the orbital file
 open(newunit=fid,file=TRIM(amoname),status='replace')
 write(fid,'(A)') "[Amesp's MO File Generated by fch2amo]"
 write(fid,'(/,A,/,A,I0)') '[Atoms]', 'Natm= ', natom
 coor = coor/Bohr_const

 if(ecp) then
  do i = 1, natom, 1
   write(fid,'(A2,2X,I3,3(1X,F18.8))') elem(i), ielem(i)-frozen_e(i), coor(:,i)
  end do ! for i
 else
  do i = 1, natom, 1
   write(fid,'(A2,2X,I3,3(1X,F18.8))') elem(i), ielem(i), coor(:,i)
  end do ! for i
 end if

 deallocate(coor)
 write(fid,'(A,/,A,I0)') '[GTO]', 'NAtom_Type: ', nelmtyp

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
  write(fid,'(20I8)') IABS(shell_type(end_idx(i-1)+1:end_idx(i)))
  write(fid,'(A)') ' bas_num:'
  write(fid,'(20I8)') itmp

  write(fid,'(A)') ' a_Basis:'
  n = k ! copy
  do j = 1, m, 1
   write(fid,'(5F19.9)') prim_exp(n+1:n+itmp(j))
   n = n + itmp(j)
  end do ! for j

  write(fid,'(A)') ' BasCa:'
  do j = 1, m, 1
   write(fid,'(5F19.9)') contr_coeff(k+1:k+itmp(j))
   k = k + itmp(j)
  end do ! for j
  deallocate(itmp)

  write(fid,'(A,I0)') ' Bmnl: ', natmbas(i)
  write(fid,'(5F19.9)') (1d0,j=1,natmbas(i))
 end do ! for i

 deallocate(elem, ielem, atmnshl, shell_type, end_idx, prim_exp, contr_coeff, &
            natmbas)

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
  write(fid,'(A)') 'ECPBlock:'
  write(fid,'(8I8)') itmp

  p = SUM(itmp)
  allocate(jtmp(6,k), source=0)
  allocate(ktmp(10,p), source=0)
  n = 0
  do i = 1, natom, 1
   if(LPSkip(i) /= 0) cycle
   n = n + 1
   m = 0
   do j = 1, 10, 1
    n1 = KFirst(i,j); n2 = KLast(i,j)
    if(n1 == 0) exit
    jtmp(m+1,n) = n2 - n1 + 1
    m = m + 1
   end do ! for j
  end do ! for i
  write(fid,'(A)') 'ECPNum:'
  write(fid,'(8I8)') ((jtmp(j,i),j=1,6),i=1,k)
  write(fid,'(A)') 'ECPPow:'
  write(fid,'(8I8)') ((ktmp(j,i),j=1,10),i=1,p)

  write(fid,'(A)') 'ECPExp:'
  write(fid,'(A)') 'ECPCof:'
  deallocate(LPSkip)
 else
  write(fid,'(A)') 'ECP:   F'
 end if

 sph_str = 'car'
 if(sph) sph_str = 'sph'

 if((.not.uhf) .and. mult==1) then
  write(fid,'(A)') '[MO] '//sph_str//' C noRI'
  write(fid,'(A,I0)') 'Nocc= ', na
  write(fid,'(A)') 'En:'
  write(fid,'(5(1X,E19.12))') eigen_e_a
  write(fid,'(A)') 'MoCu:'
  write(fid,'(5(1X,E19.12))') alpha_coeff
 else
  write(fid,'(A)') '[MO] '//sph_str//' O noRI'
  write(fid,'(A,I0)') 'NoccA= ', na
  write(fid,'(A)') 'EnA:'
  write(fid,'(5(1X,E19.12))') eigen_e_a
  write(fid,'(A)') 'MoCuA:'
  write(fid,'(5(1X,E19.12))') alpha_coeff
  write(fid,'(A,I0)') 'NoccB= ', nb
  write(fid,'(A)') 'EnB:'
  if(uhf) then
   write(fid,'(5(1X,E19.12))') eigen_e_b
  else
   write(fid,'(5(1X,E19.12))') eigen_e_a
  end if
  write(fid,'(A)') 'MoCuB:'
  if(uhf) then
   write(fid,'(5(1X,E19.12))') beta_coeff
  else
   write(fid,'(5(1X,E19.12))') alpha_coeff
  end if
 end if

 close(fid)
 deallocate(eigen_e_a, alpha_coeff)
 if(uhf) deallocate(eigen_e_b, beta_coeff)
end subroutine fch2amo

