! written by jxzou at 20230729: a wrapper of mkl2fch and fch2cfour for ORCA ->
! CFOUR

program mkl2cfour
 use util_wrapper, only: mkl2fch_wrap
 implicit none
 integer :: i, k, irel, SYSTEM
 character(len=6) :: str
 character(len=240) :: mklname, fchname

 i = iargc()
 if(i<1 .or. i>2) then
  write(6,'(/,A)') ' ERROR in program mkl2cfour: wrong command line argument!'
  write(6,'(A)')   ' Example 1: mkl2cfour a.mkl'
  write(6,'(A,/)') ' Example 2: mkl2cfour a.mkl -sfx2c'
  stop
 end if

 mklname = ' '
 call getarg(1, mklname)
 call require_file_exist(mklname)

 irel = -1
 if(i > 1) then
  call getarg(2, str)
  if(TRIM(str) == '-sfx2c') then
   irel = -3
  else
   write(6,'(/,A)') 'ERROR in program mkl2cfour: wrong command line argument!'
   write(6,'(A)')   'Example: mkl2cfour a.mkl -sfx2c'
   stop
  end if
 end if

 fchname = ' '
 call find_specified_suffix(mklname, '.mkl', i)
 call get_a_random_int(k)
 write(fchname,'(A,I0,A)') mklname(1:i-1)//'_',k,'.fch'

 call mkl2fch_wrap(mklname=mklname,fchname=fchname,irel=irel)

 i = SYSTEM('fch2cfour '//TRIM(fchname))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in program mkl2cfour: failed to call utility fch2cfour.'
  write(6,'(A)') 'mkl2cfour is a wrapper of mkl2fch and fch2cfour, so fch2cfour&
                 & is required'
  write(6,'(A)') 'to be called successfully.'
  stop
 end if

 call delete_file(fchname)
end program mkl2cfour

