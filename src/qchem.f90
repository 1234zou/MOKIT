! written by jxzou at 20230217

! Q-Chem .fch(k) -> AMESP (.aip, .amo)
subroutine qchem2amesp(fchname, aipname)
 implicit none
 integer :: i, system, RENAME
 character(len=240) :: std_fch, std_aip, std_orb, orbname
 character(len=240), intent(in) :: fchname, aipname

 i = index(aipname, '.aip', back=.true.)
 if(i == 0) then
  write(6,'(/,A)') "ERROR in subroutine qchem2amesp: aipname must include '.aip&
                   &' as suffix!"
  stop
 end if
 orbname = aipname(1:i-1)//'.amo'

 i = index(fchname, '.', back=.true.)
 std_fch = fchname(1:i-1)//'_std.fch'
 std_aip = fchname(1:i-1)//'_std.aip'
 std_orb = fchname(1:i-1)//'_std.amo'

 call standardize_fch(fchname)
 i = system('fch2amo '//TRIM(std_fch))
 call delete_file(std_fch)
 i = RENAME(TRIM(std_aip), TRIM(aipname))
 i = RENAME(TRIM(std_orb), TRIM(orbname))
end subroutine qchem2amesp

! Q-Chem .fch(k) -> BDF (.inp, .scforb, .BAS)
subroutine qchem2bdf(fchname, inpname)
 implicit none
 integer :: i, system, RENAME
 character(len=240) :: std_fch, std_inp, std_orb, orbname
 character(len=240), intent(in) :: fchname, inpname

 i = index(inpname, '.inp', back=.true.)
 if(i == 0) then
  write(6,'(/,A)') "ERROR in subroutine qchem2bdf: inpname must include '.inp'&
                   & as suffix!"
  stop
 end if
 orbname = inpname(1:i-1)//'.scforb'

 i = index(fchname, '.', back=.true.)
 std_fch = fchname(1:i-1)//'_std.fch'
 std_inp = fchname(1:i-1)//'_std_bdf.inp'
 std_orb = fchname(1:i-1)//'_std_bdf.scforb'

 call standardize_fch(fchname)
 i = system('fch2bdf '//TRIM(std_fch))
 call delete_file(std_fch)
 i = RENAME(TRIM(std_inp), TRIM(inpname))
 i = RENAME(TRIM(std_orb), TRIM(orbname))
end subroutine qchem2bdf

! Q-Chem .fch(k) -> CFOUR (ZMAT, OLDMOS, GENBAS, ECPDATA)
subroutine qchem2cfour(fchname)
 implicit none
 integer :: i, system
 character(len=240) :: std_fch
 character(len=240), intent(in) :: fchname

 i = index(fchname, '.', back=.true.)
 std_fch = fchname(1:i-1)//'_std.fch'
 call standardize_fch(fchname)
 i = system('fch2cfour '//TRIM(std_fch))
 call delete_file(std_fch)
end subroutine qchem2cfour

! Q-Chem .fch(k) -> Dalton (.dal, .mol)
subroutine qchem2dalton(fchname, dalname)
 implicit none
 integer :: i, system, RENAME
 character(len=240) :: std_fch, std_dal, std_mol, molname
 character(len=240), intent(in) :: fchname, dalname

 i = index(dalname, '.dal', back=.true.)
 if(i == 0) then
  write(6,'(/,A)') "ERROR in subroutine qchem2dalton: dalname must include '.d&
                   &al' as suffix!"
  stop
 end if
 molname = dalname(1:i-1)//'.mol'

 i = index(fchname, '.', back=.true.)
 std_fch = fchname(1:i-1)//'_std.fch'
 std_dal = fchname(1:i-1)//'_std.dal'
 std_mol = fchname(1:i-1)//'_std.mol'

 call standardize_fch(fchname)
 i = system('fch2dal '//TRIM(std_fch))
 call delete_file(std_fch)
 i = RENAME(TRIM(std_dal), TRIM(dalname))
 i = RENAME(TRIM(std_mol), TRIM(molname))
end subroutine qchem2dalton

! Q-Chem .fch(k) -> GAMESS (.inp)
subroutine qchem2gms(fchname, inpname)
 implicit none
 integer :: i, system, RENAME
 character(len=240) :: std_fch, std_inp
 character(len=240), intent(in) :: fchname, inpname

 i = index(fchname, '.', back=.true.)
 std_fch = fchname(1:i-1)//'_std.fch'
 std_inp = fchname(1:i-1)//'_std.inp'
 call standardize_fch(fchname)
 i = system('fch2inp '//TRIM(std_fch))
 call delete_file(std_fch)
 i = RENAME(TRIM(std_inp), TRIM(inpname))
end subroutine qchem2gms

! Q-Chem .fch(k) -> (Open)Molcas (.inporb, .INPORB)
subroutine qchem2molcas(fchname, inpname)
 implicit none
 integer :: i, system, RENAME
 character(len=240) :: std_fch, std_inp
 character(len=240), intent(in) :: fchname, inpname

 i = index(inpname, '.input', back=.true.)
 if(i == 0) then
  write(6,'(/,A)') "ERROR in subroutine qchem2molcas: inpname must include '.in&
                   &p' as suffix!"
  stop
 end if

 i = index(fchname, '.', back=.true.)
 std_fch = fchname(1:i-1)//'_std.fch'
 std_inp = fchname(1:i-1)//'_std.input'

 call standardize_fch(fchname)
 i = system('fch2inporb '//TRIM(std_fch))
 call delete_file(std_fch)
 i = RENAME(TRIM(std_inp), TRIM(inpname))
end subroutine qchem2molcas

! Q-Chem .fch(k) -> Molpro (.com, .a, .b)
subroutine qchem2molpro(fchname, inpname)
 implicit none
 integer :: i, system, RENAME
 character(len=240) :: std_fch, std_inp
 character(len=240), intent(in) :: fchname, inpname

 i = index(inpname, '.com', back=.true.)
 if(i == 0) then
  write(6,'(/,A)') "ERROR in subroutine qchem2molpro: inpname must include '.c&
                   &om' as suffix!"
  stop
 end if

 i = index(fchname, '.', back=.true.)
 std_fch = fchname(1:i-1)//'_std.fch'
 std_inp = fchname(1:i-1)//'_std.com'

 call standardize_fch(fchname)
 i = system('fch2com '//TRIM(std_fch))
 call delete_file(std_fch)
 i = RENAME(TRIM(std_inp), TRIM(inpname))
end subroutine qchem2molpro

! Q-Chem .fch(k) -> PSI4 (.inp, .A, .B)
subroutine qchem2psi(fchname, inpname)
 implicit none
 integer :: i, system, RENAME
 character(len=240) :: std_fch, std_inp
 character(len=240), intent(in) :: fchname, inpname

 i = index(inpname, '.inp', back=.true.)
 if(i == 0) then
  write(6,'(/,A)') "ERROR in subroutine qchem2psi: inpname must include '.inp'&
                   & as suffix!"
  stop
 end if

 i = index(fchname, '.', back=.true.)
 std_fch = fchname(1:i-1)//'_std.fch'
 std_inp = fchname(1:i-1)//'_std_psi.inp'

 call standardize_fch(fchname)
 i = system('fch2psi '//TRIM(std_fch))
 call delete_file(std_fch)
 i = RENAME(TRIM(std_inp), TRIM(inpname))
end subroutine qchem2psi

! Q-Chem .fch(k) -> ORCA (.inp, .mkl, .gbw)
subroutine qchem2orca(fchname, inpname)
 use util_wrapper, only: mkl2gbw
 implicit none
 integer :: i, system, RENAME
 character(len=240) :: std_fch, std_inp, std_mkl, mklname, gbwname
 character(len=240), intent(in) :: fchname, inpname

 i = index(inpname, '.inp', back=.true.)
 if(i == 0) then
  write(6,'(/,A)') "ERROR in subroutine qchem2orca: inpname must include '.inp'&
                   & as suffix!"
  stop
 end if
 mklname = inpname(1:i-1)//'.mkl'
 gbwname = inpname(1:i-1)//'.gbw'

 i = index(fchname, '.', back=.true.)
 std_fch = fchname(1:i-1)//'_std.fch'
 std_inp = fchname(1:i-1)//'_std_o.inp'
 std_mkl = fchname(1:i-1)//'_std_o.mkl'

 call standardize_fch(fchname)
 i = system('fch2mkl '//TRIM(std_fch))
 call delete_file(std_fch)
 i = RENAME(TRIM(std_inp), TRIM(inpname))
 i = RENAME(TRIM(std_mkl), TRIM(mklname))
 call mkl2gbw(mklname, gbwname)
end subroutine qchem2orca

! standardize a .fch file, usually for Q-Chem generated .fchk file
subroutine standardize_fch(fchname)
 use fch_content
 implicit none
 integer :: i, ne1, ne2
 character(len=240) :: new_fch
 character(len=240), intent(in) :: fchname
 logical, external :: has_pople_sp

 new_fch = ' '
 i = index(fchname, '.', back=.true.)
 new_fch = fchname(1:i-1)//'_std.fch'

 call read_charge_and_mult_from_fch(fchname, charge, mult)
 call read_na_and_nb_from_fch(fchname, na, nb)
 call read_natom_from_fch(fchname, natom)
 allocate(ielem(natom), coor(3,natom))
 call read_nuc_from_fch(fchname, natom, ielem)
 call read_coor_from_fch(fchname, natom, coor)

 call read_ncontr_from_fch(fchname, ncontr)
 allocate(shell_type(ncontr), prim_per_shell(ncontr), shell2atom_map(ncontr))
 call read_shltyp_and_shl2atm_from_fch(fchname, ncontr, shell_type, &
                                       shell2atom_map)
 call read_prim_per_shell(fchname, ncontr, prim_per_shell)

 call read_nprim_from_fch(fchname, nprim)
 allocate(prim_exp(nprim), contr_coeff(nprim))
 call read_prim_and_contr_coeff_from_fch(fchname, nprim, prim_exp, contr_coeff)
 if(has_pople_sp(fchname)) then
  allocate(contr_coeff_sp(nprim))
  call read_contr_coeff_sp_from_fch(fchname, nprim, contr_coeff_sp)
 end if

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(eigen_e_a(nif), alpha_coeff(nbf,nif), tot_dm(nbf,nbf))
 call read_eigenvalues_from_fch(fchname, nif, 'a', eigen_e_a)
 call read_mo_from_fch(fchname, nbf, nif, 'a', alpha_coeff)
 call read_density_from_fch(fchname, 1, nbf, tot_dm)

 call check_uhf_in_fch(fchname, is_uhf)
 if(is_uhf) then
  allocate(eigen_e_b(nif), beta_coeff(nbf,nif), spin_dm(nbf,nbf))
  call read_eigenvalues_from_fch(fchname, nif, 'b', eigen_e_b)
  call read_mo_from_fch(fchname, nbf, nif, 'b', beta_coeff)
  call read_density_from_fch(fchname, 2, nbf, spin_dm)
 end if

 call write_fch(new_fch)

 ! check whether ECP/PP is used
 call calc_ne_from_fch(fchname, ne1)
 call read_ne_from_fch(fchname, ne2)
 if(ne1 > ne2) then
  write(6,'(A)') REPEAT('-',79)
  write(6,'(A)') 'Warning from subroutine standardize_fch: ECP/PP is probably u&
                 &sed during calculation.'
  write(6,'(A)') 'The .fch(k) file generated by Q-Chem does not include any ECP&
                 &/PP data. So you'
  write(6,'(A)') 'need to write/add the ECP/PP information in the generated inp&
                 &ut file. For example:'
  write(6,'(/,A)') '#p RHF genecp nosymm int=nobasistransform guess=read geom=a&
                   &llcheck'
  write(6,'(/,A,/,A,/,A)') 'H C 0','cc-pVDZ','****'
  write(6,'(A,/,A,/,A)') 'I 0','SDD','****'
  write(6,'(/,A,/,A,//)') 'I 0', 'SDD'
  write(6,'(A)') REPEAT('-',79)
 end if
end subroutine standardize_fch

subroutine read_prim_per_shell(fchname, ncontr, prim_per_shell)
 implicit none
 integer :: i, fid
 integer, intent(in) :: ncontr
 integer, intent(out) :: prim_per_shell(ncontr)
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname

 prim_per_shell = 0
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:24) == 'Number of primitives per') exit
 end do ! for while

 read(fid,'(6(1X,I11))') (prim_per_shell(i), i=1,ncontr)
 close(fid)
end subroutine read_prim_per_shell

! Check whether there exists 'P(S=P) Contraction coefficients' section in a
! specified .fch file
! If this section exists but contains all zero elements, has_sp is still .T.
function has_pople_sp(fchname) result(has_sp)
 implicit none
 integer :: fid
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
 logical :: has_sp

 has_sp = .false.
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:6) == 'P(S=P)') then
   has_sp = .true.
   exit
  end if
  if(buf(1:19) == 'Coordinates of each') exit
 end do ! for while

 close(fid)
end function has_pople_sp

subroutine read_contr_coeff_sp_from_fch(fchname, nprim, contr_coeff_sp)
 implicit none
 integer :: i, fid
 integer, intent(in) :: nprim
 real(kind=8), intent(out) :: contr_coeff_sp(nprim)
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:6) == 'P(S=P)') exit
 end do ! for while

 if(i /= 0) then
  close(fid)
  write(6,'(A)') "ERROR in subroutine read_contr_coeff_sp_from_fch: no 'P(S=P)&
                 &' found in file "//TRIM(fchname)
  stop
 end if

 read(fid,'(5(1X,ES15.8))') (contr_coeff_sp(i),i=1,nprim)
 close(fid)
end subroutine read_contr_coeff_sp_from_fch

subroutine read_nprim_from_fch(fchname, nprim)
 implicit none
 integer :: fid
 integer, intent(out) :: nprim
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname

 nprim = 0
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:14) == 'Number of prim') exit
 end do ! for while

 BACKSPACE(fid)
 read(fid,'(A49,2X,I10)') buf, nprim
 close(fid)
end subroutine read_nprim_from_fch

subroutine read_prim_and_contr_coeff_from_fch(fchname, nprim, prim_exp, contr_coeff)
 implicit none
 integer :: i, fid
 integer, intent(in) :: nprim
 real(kind=8), intent(out) :: prim_exp(nprim), contr_coeff(nprim)
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:13) == 'Primitive exp') exit
 end do ! for while

 read(fid,'(5(1X,ES15.8))') (prim_exp(i),i=1,nprim)
 read(fid,'(A)') buf
 read(fid,'(5(1X,ES15.8))') (contr_coeff(i),i=1,nprim)
 close(fid)
end subroutine read_prim_and_contr_coeff_from_fch

! Calculate the total number of electrons using elements in a .fch(k) file.
! This value is not affected by ECP/PP.
subroutine calc_ne_from_fch(fchname, ne)
 implicit none
 integer :: charge, mult, natom
 integer, intent(out) :: ne
 integer, allocatable :: nuc(:)
 character(len=240), intent(in) :: fchname

 ne = 0
 call read_charge_and_mult_from_fch(fchname, charge, mult)
 call read_natom_from_fch(fchname, natom)
 allocate(nuc(natom))
 call read_nuc_from_fch(fchname, natom, nuc)
 ne = SUM(nuc) - charge
end subroutine calc_ne_from_fch

! read the type of HF from a .fch file
subroutine read_hf_type_from_fch(fchname, hf_type)
 use fch_content, only: check_ghf_in_fch, check_uhf_in_fch
 implicit none
 integer :: fid, na, nb
 integer, intent(out) :: hf_type
!f2py intent(out) :: hf_type
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname
 logical :: ghf, uhf

 hf_type = 1 ! default: real RHF, IOP(3/116=1)

 call check_ghf_in_fch(fchname, ghf)
 if(ghf) then
  hf_type = 7 ! complex GHF, IOP(3/116=7)
  return
 end if

 call check_uhf_in_fch(fchname, uhf)
 if(uhf) then
  hf_type = 2 ! real UHF, IOP(3/116=2)
  return
 end if

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:15) == 'Number of alpha') exit
 end do ! for while

 BACKSPACE(fid)
 read(fid,'(A49,2X,I10)') buf, na
 read(fid,'(A49,2X,I10)') buf, nb
 close(fid)

 if(na /= nb) hf_type = 101 ! real ROHF, IOP(3/116=101)
end subroutine read_hf_type_from_fch

