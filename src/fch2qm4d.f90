! written by jxzou at 20210328: adjust the orders of d,f,g, etc functions in
!  .fch(k) file, to the orders in QM4D

! Note:
! 1) QM4D supports only Cartesian-type basis set (6D 10F)
! 2) For ECP integrals, QM4D use spherical harmonic functions. For example,
!    energies of RHF/LANL2DZ for Zn atom is different in QM4D vs Gaussian,
!    but MOs are the same.

! Restrictions:
! 1) Atoms belong to an element must have only one basis set
!    C1 cc-pVTZ with C2 cc-pVDZ is not supported
! 2) Background charges are not stored in .fch(k) file, so the generated .inp
!    file will not contain background charges, neither
! 3) Currently only HF/DFT wavefunction is supported

! The 'Shell types' array in Gaussian .fch file:
!
!   Spherical     |     Cartesian
! -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5
!  H  G  F  D  L  S  P  D  F  G  H
!
! 'L' is 'SP' in Pople-type basis sets

program main
 implicit none
 integer :: i
 integer, parameter :: iout = 6
 character(len=4) :: str
 character(len=240) :: fname
 logical :: binary

 i = iargc()
 if(i<1 .or. i>2) then
  write(iout,'(/,1X,A)') 'ERROR in subroutine fch2qm4d: wrong command line arguments!'
  write(iout,'(1X,A)') 'Example 1 (for HF/DFT): fch2qm4d a.fch'
  write(iout,'(1X,A)') 'Example 2 (for HF/DFT): fch2qm4d a.fch -xml'
  write(iout,'(1X,A,/)') 'Example 3 (for HF/DFT): fch2qm4d a.fch -bin'
  stop
 end if

 fname = ' '
 call getarg(1, fname)
 call require_file_exist(fname)

 binary = .false.
 ! False: use .xml ASCII text file
 ! True : use .bin binary file

 if(i == 2) then
  str = ' '
  call getarg(2, str)
  select case(str)
  case('-xml')
  case('-bin')
   binary = .true.
  case default
   write(iout,'(/,A)') 'ERROR in subroutine fch2qm4d: the 2nd command line&
                      & argument is wrong.'
   write(iout,'(A)') "Only '-xml' or '-bin' is accepted. But got '"//str//"'"
   stop
  end select
 end if

 call fch2qm4d(fname, binary)
 stop
end program main

! read the MOs in .fch(k) file and adjust its d,f,g, etc. functions order
!  of Gaussian to that of QM4D
subroutine fch2qm4d(fchname, binary)
 use fch_content, only: check_uhf_in_fch, read_fch, nif, nbf, ncontr, shell_type,&
  natom, coor, elem, na, nb, alpha_coeff, beta_coeff
 implicit none
 integer :: i, j, k, fid
 integer :: n6dmark, n10fmark
 integer, parameter :: iout = 6
 integer, allocatable :: d_mark(:), f_mark(:)
 character(len=240) :: inpname ! QM4D input file
 character(len=240) :: xyzname ! .xyz file to hold Cartesian coordinates
 character(len=240) :: xmlname ! .xml file to hold density matrix
 character(len=240) :: binname ! .bin file to hold MOs (binary)
 character(len=240), intent(in) :: fchname
 real(kind=8), allocatable :: coeff(:,:), alpha_occ(:,:), beta_occ(:,:)
 real(kind=8), allocatable :: alpha_dm(:,:), beta_dm(:,:)
 logical :: uhf
 logical, intent(in) :: binary

 if(binary) then
  write(iout,'(A)') 'ERROR in subroutine fch2qm4d: currently binary=.True.&
                   & is not supported.'
  write(iout,'(A)') "Please use '-xml'."
  stop
 end if

 call check_uhf_in_fch(fchname, uhf) ! determine whether UHF
 call read_fch(fchname, uhf)         ! read content in .fch

 ! check if any Cartesian functions
 if( ANY(shell_type < -1) ) then
  write(iout,'(A)') 'ERROR in subroutine fch2qm4d: spherical harmonic functions&
                   & detected in file '//TRIM(fchname)//'.'
  write(iout,'(A)') "QM4D currently supports only Cartesian functions. You need to add&
                  & '6D 10F' keywords in Gaussian."
  stop
 end if

 ! check any angular momentum higher than F
 if( ANY(shell_type > 3) ) then
  write(iout,'(A)') 'ERROR in subroutine fch2qm4d: it seems that QM4D does not&
                   & support G, H functions.'
  write(iout,'(A)') 'You can use a smaller basis set (e.g. cc-pVQZ -> cc-pVTZ).'
  stop
 end if

 i = index(fchname,'.fch',back=.true.)
 inpname = fchname(1:i-1)//'.inp'
 xyzname = fchname(1:i-1)//'.xyz'
 xmlname = fchname(1:i-1)//'_0.xml'
 binname = fchname(1:i-1)//'_0.bin'

 ! create the .xyz file
 call write_xyzfile(natom, coor, elem, xyzname)
 deallocate(coor)

 ! create the QM4D input file
 call write_qm4d_inp(inpname, uhf, binary)

 ! create basis set files *.gen for each element
 call write_basis_files_for_each_elem()

 if(uhf) then
  allocate(coeff(nbf,2*nif), source=0d0)
  coeff(1:nbf,1:nif) = alpha_coeff
  coeff(1:nbf,nif+1:2*nif) = beta_coeff
  nif = 2*nif
 else
  allocate(coeff(nbf,nif), source=alpha_coeff)
 end if
 deallocate(alpha_coeff)

! then we adjust the basis functions in each MO according to the type of basis functions
 n6dmark = 0
 n10fmark = 0
 allocate(d_mark(ncontr), f_mark(ncontr))
 d_mark = 0
 f_mark = 0
 nbf = 0
 do i = 1, ncontr, 1
  select case(shell_type(i))
  case( 0)   ! S
   nbf = nbf + 1
  case( 1)   ! 3P
   nbf = nbf + 3
  case(-1)   ! SP or L
   nbf = nbf + 4
  case( 2)   ! 6D
   n6dmark = n6dmark + 1
   d_mark(n6dmark) = nbf + 1
   nbf = nbf + 6
  case( 3)   ! 10F
   n10fmark = n10fmark + 1
   f_mark(n10fmark) = nbf + 1
   nbf = nbf + 10
  case default
   write(iout,'(A)') 'ERROR in subroutine fch2qm4d: shell_type out of range!'
   write(iout,'(2(A,I0))') 'i=', i, ', shell_type(i)=', shell_type(i)
   stop
  end select
 end do ! for i

 deallocate(shell_type)

 ! adjust the order of d, f, etc. functions
 do i = 1, n6dmark, 1
  call fch2qm4d_permute_6d(nif,coeff(d_mark(i):d_mark(i)+5,:))
 end do
 do i = 1, n10fmark, 1
  call fch2qm4d_permute_10f(nif,coeff(f_mark(i):f_mark(i)+9,:))
 end do
! adjustment finished

 deallocate(d_mark, f_mark)

 ! initilize occupation number arrays
 if(uhf) then
  nif = nif/2
  allocate(alpha_occ(nif,nif), source=0d0)
  forall(i=1:na) alpha_occ(i,i) = 1d0
  allocate(beta_occ(nif,nif), source=0d0)
  forall(i=1:nb) beta_occ(i,i) = 1d0
 else
  allocate(alpha_occ(nif,nif), source=0d0)
  forall(i=1:nb) alpha_occ(i,i) = 2d0
  if(na > nb) forall(i=nb+1:na) alpha_occ(i,i) = 1d0
 end if

 ! calculate alpha density (and beta density, if any)
 alpha_coeff = coeff(:,1:nif)
 allocate(alpha_dm(nbf,nbf))
 alpha_dm = MATMUL(MATMUL(alpha_coeff,alpha_occ),TRANSPOSE(alpha_coeff))
 deallocate(alpha_coeff)
 if(uhf) then
  beta_coeff = coeff(:,nif+1:2*nif)
  allocate(beta_dm(nbf,nbf))
  beta_dm = MATMUL(MATMUL(beta_coeff,beta_occ),TRANSPOSE(beta_coeff))
  deallocate(beta_coeff)
 end if
 deallocate(coeff)

 ! print alpha density and beta density into _0.xml file
 open(newunit=fid,file=TRIM(xmlname),status='replace')
 write(fid,'(A)') "<?xml version=""1.0"" encoding=""ISO-8859-1""?>"
 write(fid,'(A)') '<QM4D>'
 write(fid,'(A)') ' <!--The symmetric matrix is saved as the lower triangle format.-->'
 write(fid,'(2(A,I6),A)') " <Alpha_Density_Matrix Dimension=""", nbf, """ Data=""",&
                          nbf*(nbf+1)/2,""">"
 k = 0
 do i = 1, nbf, 1
  do j = 1, i, 1
   k = k + 1
   write(fid,'(A,I6,A,1X,ES23.16,A)') "  <Data id=""",k,""">",alpha_dm(j,i),'</Data>'
  end do ! for j
 end do ! for i

 deallocate(alpha_dm)
 write(fid,'(A)') ' </Alpha_Density_Matrix>'

 write(fid,'(A)') ' <!--The full matrix is saved as nRow*nCol.-->'
 write(fid,'(2(A,I6),A)') " <Alpha_Occup_No Row=""     1"" Col=""",nif,""" Data=""",nif,""">"
 do i = 1, nif, 1
  write(fid,'(A,I6,A,1X,ES23.16,A)') "  <Data id=""",i,""">",alpha_occ(i,i),'</Data>'
 end do ! for i
 deallocate(alpha_occ)
 write(fid,'(A)') '  </Alpha_Occup_No>'

 if(uhf) then
  write(fid,'(A)') ' <!--The symmetric matrix is saved as the lower triangle format.-->'
  write(fid,'(2(A,I6),A)') " <Beta_Density_Matrix Dimension=""", nbf, """ Data=""",&
                           nbf*(nbf+1)/2,""">"
  k = 0
  do i = 1, nbf, 1
   do j = 1, i, 1
    k = k + 1
    write(fid,'(A,I6,A,1X,ES23.16,A)') "  <Data id=""",k,""">",beta_dm(j,i),'</Data>'
   end do ! for j
  end do ! for i
  deallocate(beta_dm)
  write(fid,'(A)') ' </Beta_Density_Matrix>'

  write(fid,'(A)') ' <!--The full matrix is saved as nRow*nCol.-->'
  write(fid,'(2(A,I6),A)') " <Beta_Occup_No Row=""     1"" Col=""",nif,""" Data=""",nif,""">"
  do i = 1, nif, 1
   write(fid,'(A,I6,A,1X,ES23.16,A)') "  <Data id=""",i,""">",beta_occ(i,i),'</Data>'
  end do ! for i
  deallocate(beta_occ)
  write(fid,'(A)') '  </Beta_Occup_No>'
 end if

 write(fid,'(A)') '</QM4D>'
 close(fid)
 return
end subroutine fch2qm4d

! write/create a QM4D .inp file
subroutine write_qm4d_inp(inpname, uhf, binary)
 use fch_content, only: charge, mult, natom, elem, LPSkip
 implicit none
 integer :: i, fid
 character(len=240) :: xml0name, xmlname, binname, xyzname
 character(len=240), intent(in) :: inpname
 logical, intent(in) :: uhf, binary

 i = index(inpname, '.inp', back=.true.)
 xml0name = inpname(1:i-1)//'_0.xml'
 xmlname = inpname(1:i-1)//'.xml'
 binname = inpname(1:i-1)//'.bin'
 xyzname = inpname(1:i-1)//'.xyz'

 open(newunit=fid,file=TRIM(inpname),status='replace')
 write(fid,'(A)') '$qm'
 write(fid,'(A)') 'xyz    '//TRIM(xyzname)
 if(uhf) then
  write(fid,'(A)') 'spin   2'
 else
  write(fid,'(A)') 'spin   1'
 end if
 write(fid,'(A,I0)') 'charge ', charge
 write(fid,'(A,I0)') 'mult   ', mult
 write(fid,'(A)') 'iter   128'
 write(fid,'(A)') 'method hf'
 write(fid,'(A)') 'etol   1.e-08'
 if(binary) then
  write(fid,'(A)') 'guess  read '//TRIM(binname)
 else
  write(fid,'(A)') 'guess  read '//TRIM(xml0name)
 end if
 write(fid,'(A)') 'diis   6 1000.0'

 do i = 1, natom, 1
  if(i > 1) then
   if(ANY(elem(1:i-1) == elem(i))) cycle
  end if

  if(allocated(LPSkip)) then
   if(LPSkip(i) == 0) then
    write(fid,'(A)') 'ecpbasis '//TRIM(elem(i))//' '//TRIM(elem(i))//'.gen'
   else
    write(fid,'(A)') 'basis  '//TRIM(elem(i))//' '//TRIM(elem(i))//'.gen'
   end if
  else ! not allocated
    write(fid,'(A)') 'basis  '//TRIM(elem(i))//' '//TRIM(elem(i))//'.gen'
  end if
 end do ! for i

 write(fid,'(A)') 'print  normal'
 write(fid,'(A)') 'end'

 write(fid,'(/,A)') '$doqm'

 write(fid,'(/,A)') '$xml'
 write(fid,'(A)') 'ow '//TRIM(xmlname)
 write(fid,'(A)') 'denmatw'
 write(fid,'(A)') 'cw'
 write(fid,'(A)') 'end'
 close(fid)
 return
end subroutine write_qm4d_inp

! write/create basis files (and ECP, if any) for each element
subroutine write_basis_files_for_each_elem()
 use fch_content, only: natom, elem, ncontr, prim_exp, contr_coeff, contr_coeff_sp,&
  shell_type, shell2atom_map, prim_per_shell, LPSkip, LMax, RNFroz, KFirst, &
  KLast, NLP, CLP, ZLP
 implicit none
 integer :: i, j, k0, k, m, n, p, q, fid
 character(len=1), parameter :: am_type(0:6)=['S','P','D','F','G','H','I']
 character(len=2) :: str
 character(len=30) :: basname

 ! print basis sets into the .inp file
 k = 0 ! current index in array shell_type
 p = 0 ! current index in array prim_exp

 do i = 1, natom, 1
  if(i > 1) then
   if(ANY(elem(1:i-1) == elem(i))) then
    if(i == natom) exit
    k0 = k
    do j = k0+1, ncontr, 1
     if(shell2atom_map(j) == i+1) exit
     k = k + 1
     p = p + prim_per_shell(j)
    end do ! for j
    cycle
   end if
  end if

  basname = TRIM(elem(i))//'.gen'
  open(newunit=fid,file=TRIM(basname),status='replace')
  write(fid,'(A)') '****'
  write(fid,'(A)') elem(i)//' 0'

  k = k + 1
  k0 = k
  do j = k0, ncontr, 1
   if(shell2atom_map(j) == i+1) exit
   m = shell_type(j); n = prim_per_shell(j)
   if(m == -1) then
    str = 'SP'
   else
    str = am_type(m)//' '
   end if
   write(fid,'(A2,1X,I2,A)') str, n, ' 1.00'

   do q = p+1, p+n, 1
    if(allocated(contr_coeff_sp)) then
     if(DABS(contr_coeff_sp(q)) > 1e-6) then
      write(fid,'(3(2X,ES15.8))') prim_exp(q), contr_coeff(q), contr_coeff_sp(q)
     else
      write(fid,'(2(2X,ES15.8))') prim_exp(q), contr_coeff(q)
     end if
    else ! no 'SP'
     write(fid,'(2(2X,ES15.8))') prim_exp(q), contr_coeff(q)
    end if
   end do ! for q

   p = p + n ! update p
  end do ! for j

  write(fid,'(A)') '****'
  k = j - 1 ! update k

  if(allocated(LPSkip)) then
   if(LPSkip(i) == 0) call prt_ecp_in_gau_format(fid, i)
  end if
  close(fid)
 end do ! for i

 deallocate(elem, shell2atom_map, prim_per_shell, prim_exp, contr_coeff)
 if(allocated(contr_coeff_sp)) deallocate(contr_coeff_sp)
 if(allocated(LPSkip)) deallocate(LPSkip)
 if(allocated(LMax)) deallocate(LMax)
 if(allocated(RNFroz)) deallocate(RNFroz)
 if(allocated(KFirst)) deallocate(KFirst)
 if(allocated(KLast)) deallocate(KLast)
 if(allocated(NLP)) deallocate(NLP)
 if(allocated(CLP)) deallocate(CLP)
 if(allocated(ZLP)) deallocate(ZLP)
 return
end subroutine write_basis_files_for_each_elem

! print ECP data in Gaussian format
subroutine prt_ecp_in_gau_format(fid, i)
 use fch_content, only: elem, LMax, RNFroz, KFirst, KLast, NLP, CLP, ZLP
 implicit none
 integer :: j, n, n1, n2
 integer, intent(in) :: fid, i
 character(len=1) :: str
 character(len=1), parameter :: am_type1(0:6)=['s','p','d','f','g','h','i']

 write(fid,'(/,A)') elem(i)//' 0'
 write(fid,'(A,I3,1X,I4)') TRIM(elem(i))//'-ECP ',LMax(i),INT(RNFroz(i))
 str = am_type1(LMax(i))

 do j = 1, 10, 1
  n1 = KFirst(i,j); n2 = KLast(i,j)
  if(n1 == 0) exit
  if(j == 1) then
   write(fid,'(A)') str//' potential'
  else
   write(fid,'(A)') am_type1(j-2)//'-'//str//' potential'
  end if
  write(fid,'(I3)') n2-n1+1
  do n = n1, n2, 1
   write(fid,'(I0,2(3X,ES15.8))') NLP(n), ZLP(n), CLP(n)
  end do ! for n
 end do ! for j

 return
end subroutine prt_ecp_in_gau_format

subroutine fch2qm4d_permute_6d(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(6) = [1, 4, 5, 2, 6, 3]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(6,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian d functions in Gaussian
! To: the order of Cartesian d functions in QM4D
! 1  2  3  4  5  6
! XX,YY,ZZ,XY,XZ,YZ
! XX,XY,XZ,YY,YZ,ZZ

 allocate(coeff2(6,nif), source=0d0)
 forall(i = 1:6) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
 return
end subroutine fch2qm4d_permute_6d

subroutine fch2qm4d_permute_10f(nif,coeff)
 implicit none
 integer :: i
 integer, parameter :: order(10) = [1, 5, 6, 4, 10, 7, 2, 9, 8, 3]
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: coeff(10,nif)
 real(kind=8), allocatable :: coeff2(:,:)
! From: the order of Cartesian f functions in Gaussian
! To: the order of Cartesian f functions in QM4D
! 1   2   3   4   5   6   7   8   9   10
! XXX,YYY,ZZZ,XYY,XXY,XXZ,XZZ,YZZ,YYZ,XYZ
! XXX,XXY,XXZ,XYY,XYZ,XZZ,YYY,YYZ,YZZ,ZZZ

 allocate(coeff2(10,nif), source=0d0)
 forall(i = 1:10) coeff2(i,:) = coeff(order(i),:)
 coeff = coeff2
 deallocate(coeff2)
 return
end subroutine fch2qm4d_permute_10f

