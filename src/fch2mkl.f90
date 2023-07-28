! written by jxzou at 20200113: convert .fch(k) file (Gaussian) to .mkl file (Molekel, ORCA)
! updated by jxzou at 20200202: ECP/PP supported
! updated by jxzou at 20200304: Pople-type basis sets supported
! updated by jxzou at 20200322: move read_fch to read_fch.f90
! updated by jxzou at 20200622: fix the bug in F, G, H (some multiply by -1); add 1 more digit for MOs
! updated by jxzou at 20200802: add $CHARGES section to .mkl file
! updated by jxzou at 20201118: detect DKH/RESC keywords in .fch(k) file
! updated by jxzou at 20201221: add 'DelGTO' for elements Rb~Rn
! updated by jxzou at 20210407: remove '-uhf', use automatic determination
! updated by jxzou at 20210730: add 'sthresh 1e-6' in %scf, in case of linear dependence

! The 'Shell types' array in Gaussian .fch file:
!
!   Spherical     |     Cartesian
! -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5
!  H  G  F  D  L  S  P  D  F  G  H
!
! 'L' is 'SP' in Pople-type basis sets

program main
 use util_wrapper, only: formchk
 implicit none
 integer :: i
 character(len=240) :: fchname

 i = iargc()
 if(i /= 1) then
  write(6,'(/,A)') ' ERROR in subroutine fch2mkl: wrong command line argument!'
  write(6,'(/,A)') ' Example (R(O)HF, UHF, CAS): fch2mkl a.fch'
  write(6,'(A,/)') ' (a_o.inp and a_o.mkl files will be generated)'
  stop
 end if

 call getarg(1, fchname)
 call require_file_exist(fchname)

 ! if .chk file provided, convert into .fch file automatically
 i = LEN_TRIM(fchname)
 if(fchname(i-3:i) == '.chk') then
  call formchk(fchname)
  fchname = fchname(1:i-3)//'fch'
 end if

 call fch2mkl(fchname)
end program main

! convert .fch(k) file (Gaussian) to .mkl file (Molekel, ORCA)
subroutine fch2mkl(fchname)
 use fch_content
 implicit none
 integer :: i, j, k, m, n, n1, n2, am, nf3mark, ng3mark, nh3mark
 integer :: fid1, fid2 ! file id of .mkl/.inp file
 integer :: rel        ! the order of DKH, or RESC
 integer, parameter :: list(10) = [2,3,4,5,6,7,8,9,10,1]
 integer, allocatable :: f3_mark(:), g3_mark(:), h3_mark(:)
 real(kind=8), allocatable :: coeff(:,:)

 ! six types of angular momentum
 character(len=1), parameter :: am_type(0:5) = ['S','P','D','F','G','H']
 character(len=1), parameter :: am_type1(0:5) = ['s','p','d','f','g','h']
 character(len=240) :: mklname, inpname
 character(len=240), intent(in) :: fchname
 logical :: uhf, ecp, X2C

 call find_specified_suffix(fchname, '.fch', i)
 mklname = fchname(1:i-1)//'_o.mkl'
 inpname = fchname(1:i-1)//'_o.inp'

 call check_nobasistransform_in_fch(fchname)
 call check_nosymm_in_fch(fchname)
 call check_DKH_in_fch(fchname, rel)
 call check_uhf_in_fch(fchname, uhf) ! determine whether UHF
 call read_fch(fchname, uhf) ! read content in .fch(k) file
 ecp = .false.
 if(LenNCZ > 0) ecp = .true.

 ! check if any Cartesian functions
 if( ANY(shell_type > 1) ) then
  write(6,'(/,A)') 'ERROR in subroutine fch2mkl: Cartesian functions detected i&
                   &n file '//TRIM(fchname)//'.'
  write(6,'(A)') "ORCA supports only spherical functions. You need to add '5D 7&
                 &F' keywords in Gaussian."
  stop
 else if( ANY(shell_type < -5) ) then
  write(6,'(/,A)') 'ERROR in subroutine fch2mkl: angular momentum too high! not&
                   & supported.'
  stop
 end if
 ! check done

 ! update MO coefficients
 if(uhf) then ! UHF
  k = 2*nif
  allocate(coeff(nbf,k))
  coeff(:,1:nif) = alpha_coeff
  coeff(:,nif+1:) = beta_coeff
 else         ! R(O)HF
  k = nif
  allocate(coeff(nbf,k), source=alpha_coeff)
 end if

 ! find F+3, G+3 and H+3 functions, multiply them by -1
 allocate(f3_mark(nbf), g3_mark(nbf), h3_mark(nbf))
 call get_bas_mark_from_shltyp(ncontr, shell_type, nbf, nf3mark, ng3mark, &
                               nh3mark, f3_mark, g3_mark, h3_mark)
 call update_mo_using_bas_mark(nbf, k, nf3mark, ng3mark, nh3mark, f3_mark, &
                               g3_mark, h3_mark, coeff)
 deallocate(f3_mark, g3_mark, h3_mark)

 if(uhf) then ! UHF
  alpha_coeff = coeff(:,1:nif)
  beta_coeff = coeff(:,nif+1:)
 else         ! R(O)HF
  alpha_coeff = coeff
 end if
 deallocate(coeff)
 ! update MO coefficients done

 ! print elements and coordinates into .mkl file
 open(newunit=fid1,file=TRIM(mklname),status='replace')
 write(fid1,'(A)') '$MKL'
 write(fid1,'(A)') '#'
 write(fid1,'(A)') '# MKL format file produced by MOKIT'
 write(fid1,'(A)') '#'
 write(fid1,'(A)') '$CHAR_MULT'
 write(fid1,'(I0,1X,I0)') charge, mult
 write(fid1,'(A,/)') '$END'

 write(fid1,'(A)') '$COORD'
 do i = 1, natom, 1
  write(fid1,'(I3,1X,3F15.8)') ielem(i), (coor(j,i), j=1,3)
 end do ! for i
 write(fid1,'(A,/)') '$END'

 write(fid1,'(A)') '$CHARGES'
 do i = 1, natom, 1
  write(fid1,'(1X,A3)') '0.0'
 end do ! for i
 write(fid1,'(A,/)') '$END'

 ! print basis sets into .mkl file (Note: mkl file contains no ECP/PP data)
 write(fid1,'(A)') '$BASIS'
 k = 0
 do i = 1, ncontr, 1
  m = shell2atom_map(i)
  if(m > 1) then
   if(shell2atom_map(i-1) == m-1) write(fid1,'(A2)') '$$'
  end if

  m = shell_type(i); n = prim_per_shell(i)

  if(m /= -1) then ! m<-1 or m=0,1

   m = IABS(m)
   write(fid1,'(I2,1X,A)') 2*m+1, am_type(m)//' 1.0'
   do j = k+1, k+n, 1
    write(fid1,'(2(2X,ES15.8))') prim_exp(j), contr_coeff(j)
   end do ! for j

  else ! m = -1, 'L' or 'SP' in Pople-type basis sets

   write(fid1,'(I2,1X,A)') 1, am_type(0)//' 1.0'
   do j = k+1, k+n, 1
    write(fid1,'(2(2X,ES15.8))') prim_exp(j), contr_coeff(j)
   end do ! for j

   write(fid1,'(I2,1X,A)') 3, am_type(1)//' 1.0'
   do j = k+1, k+n, 1
    write(fid1,'(2(2X,ES15.8))') prim_exp(j), contr_coeff_sp(j)
   end do ! for j

  end if

  k = k + n
 end do ! for i

 write(fid1,'(/,A,/)') '$END'

 ! print coordinates, basis sets and ECP/PP data into .inp file
 open(newunit=fid2,file=TRIM(inpname),status='replace')
 write(fid2,'(A)') '%pal nprocs 4 end'
 write(fid2,'(A)') '%maxcore 1000'
 if(uhf) then
  write(fid2,'(A)') '! UHF VeryTightSCF noTRAH'
 else
  if(nopen == 0) then
   write(fid2,'(A)') '! RHF VeryTightSCF noTRAH'
  else ! nopen > 0
   write(fid2,'(A)') '! ROHF VeryTightSCF noTRAH'
  end if
 end if

 select case(rel)
 case(-2) ! no relativistic
  call check_X2C_in_fch(fchname, X2C)
  if(X2C) then
   write(6,'(A)') 'Warning in subroutine fch2mkl: X2C detected.'
   write(6,'(A)') 'But ORCA does not support X2C.'
   write(6,'(A)') 'DKH2 keywords will be printed into ORCA .inp file.'
   rel = 2 ! mimic X2C as DKH2
  end if
 case(-1)
  write(6,'(A)') 'ERROR in subroutine fch2mkl: RESC keyword detected in&
                   & file '//TRIM(fchname)//'.'
  write(6,'(A)') 'But ORCA does not support RESC.'
  stop
 case(0) ! DKH0
  write(6,'(A)') 'Warning in subroutine fch2mkl: DKH0 detected in file '//TRIM(fchname)
  write(6,'(A)') 'But ORCA does not support DKH0. DKH2 keywords will&
                   & be printed into ORCA .inp file.'
  rel = 2
 case(2) ! DKH2
 case(4) ! DKHSO, DKH4 with SO
  write(6,'(A)') 'Warning in subroutine fch2mkl: DKHSO detected in&
                   & file '//TRIM(fchname)//'.'
  write(6,'(A)') 'But ORCA does not support DKHSO. DKH2 keywords will&
                   & be printed into ORCA .inp file.'
  rel = 2
 case default
  write(6,'(A)') 'ERROR in subroutine fch2mkl: rel out of range!'
  write(6,'(A,I0)') 'rel=', rel
  stop
 end select

 if(rel == 2) then
  write(fid2,'(A)') '%rel'
  write(fid2,'(A)') ' method DKH'
  write(fid2,'(A)') ' order 2'
  write(fid2,'(A)') 'end'
 end if
 write(fid2,'(A)') '%scf'
 write(fid2,'(A)') ' Thresh 1e-12'
 write(fid2,'(A)') ' Tcut 1e-14'
 if(nif < nbf)  write(fid2,'(A)') ' sthresh 1e-6'
 write(fid2,'(A)') 'end'
 write(fid2,'(A)') '%coords'
 write(fid2,'(A)') ' Units = angs'
 write(fid2,'(A,I0)') ' Charge = ', charge
 write(fid2,'(A,I0)') ' Mult = ', mult
 write(fid2,'(A)') ' Coords'

 k = 0
 do i = 1, ncontr, 1
  m = shell2atom_map(i)

  if(m == 1) then
   if(i == 1) then
    write(fid2,'(1X,A,3(1X,F16.8))') TRIM(elem(1))//'(1)', coor(:,1)
    if(ielem(1)>36 .and. ielem(1)<87) write(fid2,'(2X,A)') 'DelECP'
    write(fid2,'(2X,A)') 'NewGTO'
   end if

  else ! m > 1
   if(shell2atom_map(i-1) == m-1) then
    write(fid2,'(2X,A)') 'end'   ! print GTO end of last atom

    if(ecp) then   ! print ECP/PP data of last atom
     if(LPSkip(m-1) == 0) then
      write(fid2,'(2X,A)') 'NewECP'
      write(fid2,'(3X,A,1X,I3)') 'N_core', INT(RNFroz(m-1))
      write(fid2,'(3X,A)') 'lmax '//am_type1(LMax(m-1))
      am = 0
      do j = 1, 10, 1
       n1 = KFirst(m-1,list(j)); n2 = KLast(m-1,list(j))
       if(n1 == 0) cycle
       am = am + 1
       write(fid2,'(3X,A1,1X,I1)') am_type1(am-1), n2-n1+1
       do n = n1, n2, 1
        write(fid2,'(3X,I2,2(1X,ES15.8),1X,I1)') n-n1+1, ZLP(n), CLP(n), NLP(n)
       end do ! for n
      end do ! for j

      write(fid2,'(2X,A)') 'end'  ! in accord with 'NewECP'
     end if
    end if         ! print ECP/PP data done

    ! print coordinates of the current atom
    write(fid2,'(1X,A,I0,A1,3(1X,F16.8))') TRIM(elem(m))//'(',m,')',coor(:,m)
    if(ielem(m)>36 .and. ielem(m)<87) write(fid2,'(2X,A)') 'DelECP'
    write(fid2,'(2X,A)') 'NewGTO'
   end if
  end if

  m = shell_type(i); n = prim_per_shell(i)

  if(m /= -1) then
   m = IABS(m)
   write(fid2,'(4X,A1,1X,I3)') am_type(m), n
   do j = k+1, k+n, 1
    write(fid2,'(2X,I3,2(2X,ES15.8))') j-k,prim_exp(j), contr_coeff(j)
   end do ! for j

  else ! m = -1, 'L' or 'SP'
   write(fid2,'(4X,A1,1X,I3)') 'S', n
   do j = k+1, k+n, 1
    write(fid2,'(2X,I3,3(2X,ES15.8))') j-k, prim_exp(j), contr_coeff(j)
   end do ! for j
   write(fid2,'(4X,A1,1X,I3)') 'P', n
   do j = k+1, k+n, 1
    write(fid2,'(2X,I3,3(2X,ES15.8))') j-k, prim_exp(j), contr_coeff_sp(j)
   end do ! for j
  end if

  k = k + n
 end do ! for i

 write(fid2,'(2X,A)') 'end'  ! in accord with 'NewGTO'
 k = shell2atom_map(ncontr)
 if(k < natom) then
  do i = k+1, natom, 1
   write(fid2,'(1X,A,I0,A,3(1X,F16.8))') 'H:(',i,')',coor(:,i)
   write(fid2,'(2X,A)') 'NewGTO S 1 1 1e6 1 end'
  end do ! for i
 end if

 if(ecp) then   ! print ECP/PP data of the last atom
  if(LPSkip(natom) == 0) then
   write(fid2,'(2X,A)') 'NewECP'
   write(fid2,'(3X,A,1X,I3)') 'N_core', INT(RNFroz(natom))
   write(fid2,'(3X,A)') 'lmax '//am_type1(LMax(natom))
   am = 0
   do j = 1, 10, 1
    n1 = KFirst(natom,list(j)); n2 = KLast(natom,list(j))
    if(n1 == 0) cycle
    am = am + 1
    write(fid2,'(3X,A1,1X,I1)') am_type1(am-1), n2-n1+1
    do n = n1, n2, 1
     write(fid2,'(3X,I2,2(1X,ES15.8),1X,I1)') n-n1+1, ZLP(n), CLP(n), NLP(n)
    end do ! for n
   end do ! for j

   write(fid2,'(2X,A)') 'end'  ! in accord with 'NewECP'
  end if

  deallocate(KFirst, KLast, Lmax, LPSkip, NLP, RNFroz, CLP, ZLP)
 end if         ! print ECP/PP data done

 write(fid2,'(1X,A)') 'end'  ! in accord with ' Coords'
 write(fid2,'(A,/)') 'end'   ! in accord with '%coords'
 close(fid2)

 deallocate(ielem, elem, coor)
 deallocate(shell_type, prim_per_shell, shell2atom_map, prim_exp, contr_coeff)
 if(allocated(contr_coeff_sp)) deallocate(contr_coeff_sp)

 ! print Alpha MO and corresponding energies into .mkl file
 write(fid1,'(A)') '$COEFF_ALPHA'
 k = 0
 do while(.true.)
  if(k+1 > nif) exit
  if(k+5 > nif) then
   j = nif
  else
   j = k + 5
  end if
  write(fid1,'(5(A4,1X))') (' a1g', i=k+1,j)
  write(fid1,'(5(F14.8,1X))') (eigen_e_a(i),i=k+1,j)
  do i = 1, nbf, 1
   write(fid1,'(5(ES15.8,1X))') (alpha_coeff(i,m),m=k+1,j)
  end do ! for i
  k = j
 end do ! for while

 write(fid1,'(A,/)') '$END'
 deallocate(alpha_coeff, eigen_e_a)

 ! print Alpha orbital occupation numbers into .mkl file
 write(fid1,'(A)') '$OCC_ALPHA'
 allocate(eigen_e_a(nif), source=0d0)
 if(uhf) then
  forall(i = 1:na) eigen_e_a(i) = 1d0
 else ! .not. uhf
  if(nopen == 0) then
   forall(i = 1:na) eigen_e_a(i) = 2d0
  else ! nopen > 0
   forall(i = 1:nb)    eigen_e_a(i) = 2d0
   forall(i = nb+1:na) eigen_e_a(i) = 1d0
  end if
 end if
 write(fid1,'(5(F12.7,1X))') (eigen_e_a(i), i=1,nif)
 deallocate(eigen_e_a)
 write(fid1,'(A,/)') '$END'

 if(uhf) then
  ! print Beta MO (if any) and corresponding energies into .mkl file
  write(fid1,'(A)') '$COEFF_BETA'
  k = 0
  do while(.true.)
   if(k+1 > nif) exit
   if(k+5 > nif) then
    j = nif
   else
    j = k + 5
   end if
   write(fid1,'(5(A4,1X))') (' a1g', i=k+1,j)
   write(fid1,'(5(F14.8,1X))') (eigen_e_b(i),i=k+1,j)
   do i = 1, nbf, 1
    write(fid1,'(5(ES15.8,1X))') (beta_coeff(i,m),m=k+1,j)
   end do ! for i
   k = j
  end do ! for while
  write(fid1,'(A,/)') '$END'
  deallocate(beta_coeff, eigen_e_b)

  ! print Beta orbital occupation numbers (if any) into .mkl file
  write(fid1,'(A)') '$OCC_BETA'
  allocate(eigen_e_b(nif), source=0d0)
  forall(i = 1:nb) eigen_e_b(i) = 1d0
  write(fid1,'(5(F12.7,1X))') (eigen_e_b(i), i=1,nif)
  deallocate(eigen_e_b)
  write(fid1,'(A,/)') '$END'
 end if

 close(fid1)
end subroutine fch2mkl

