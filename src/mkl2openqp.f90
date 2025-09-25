! written by jxzou at 20250922: a wrapper of mkl2fch and fch2openqp for ORCA ->
! OpenQP
! Note: use with caution since .mkl file does not have ECP/PP infomation

program mkl2openqp
 use util_wrapper, only: mkl2fch_wrap, fch_sph2cart_wrap, fch2openqp_wrap
 implicit none
 integer :: i, narg, irel, sf_type, charge, mult
 character(len=8) :: str8
 character(len=240) :: mklname, sph_fch, cart_fch, inpname

 narg = iargc()
 if(narg<1 .or. narg>2) then
  write(6,'(/,A)')' ERROR in program mkl2openqp: wrong command line argument!'
  write(6,'(A)')  ' Example 1 (R(O)HF, UHF): mkl2inporb h2o.mkl'
  write(6,'(A)')  ' Example 2 (SF-CIS)     : mkl2inporb h2o.mkl -sfcis'
  write(6,'(A)')  ' Example 3 (SF-TDDFT)   : mkl2inporb h2o.mkl -sf'
  write(6,'(A)')  ' Example 4 (MRSF-CIS)   : mkl2inporb h2o.mkl -mrsfcis'
  write(6,'(A,/)')' Example 5 (MRSF-TDDFT) : mkl2inporb h2o.mkl -mrsf'
  stop
 end if

 mklname = ' '
 call getarg(1, mklname)
 call require_file_exist(mklname)

 sf_type = 0
 if(narg == 2) then
  call getarg(2, str8)
  select case(TRIM(str8))
  case('-sfcis')
   sf_type = 1
  case('-sf')
   sf_type = 2
  case('-mrsfcis')
   sf_type = 3
  case('-mrsf')
   sf_type = 4
  case default
   write(6,'(/,A)') 'ERROR in program mkl2openqp: the 2nd command line argument&
                    & can only be'
   write(6,'(A)') "'-sfcis', '-sf', '-mrsfcis, '-mrsf'. But got "//TRIM(str8)
   stop
  end select

 end if

 call read_charge_and_mult_from_mkl(mklname, charge, mult)
 if((sf_type==1 .or. sf_type==2) .and. mult<3) then
  write(6,'(/,A)') 'ERROR in program mkl2openqp: spin-flip methods requires the&
                   & spin multiplicity'
  write(6,'(A,I0,A)') 'to be >=3. But got mult=',mult,' in file '//TRIM(mklname)
  stop
 end if
 if((sf_type==3 .or. sf_type==4) .and. mult/=3) then
  write(6,'(/,A)') 'ERROR in program mkl2openqp: MRSF methods requires the spin&
                   & multiplicity to'
  write(6,'(A,I0,A)') 'be 3. But got mult=',mult,' in file '//TRIM(mklname)
  stop
 end if

 call find_specified_suffix(mklname, '.mkl', i)
 inpname = mklname(1:i-1)//'.inp'

 irel = -1
 call find_specified_suffix(mklname, '.mkl', i)
 sph_fch = mklname(1:i-1)//'_sph.fch'
 cart_fch = mklname(1:i-1)//'_cart.fch'

 call mkl2fch_wrap(mklname=mklname,fchname=sph_fch,irel=irel)
 write(6,'(/,A)') REPEAT('-',79)
 write(6,'(A)') 'Remark from program mkl2openqp: ORCA supports only spherical h&
                &armonic type'
 write(6,'(A)') 'basis functions. OpenQP supports only pure Cartesian type basi&
                &s functions.'
 write(6,'(A)') 'Therefore the sph -> cart conversion is needed. This will be i&
                &nvoked automatically.'
 write(6,'(A)') REPEAT('-',79)
 call fch_sph2cart_wrap(sph_fch, cart_fch)
 write(6,'(A)') 'Done conversion. Now generate OpenQP .inp and .json files...'
 call fch2openqp_wrap(cart_fch, sf_type, inpname)
 call delete_files(2, [sph_fch, cart_fch])
end program mkl2openqp

