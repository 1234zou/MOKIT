! written by jxzou at 20250411: a wrapper of mkl2fch and fch2mrcc for ORCA ->
! Gaussian

program mkl2mrcc
 use util_wrapper, only: mkl2fch_wrap, fch2mrcc_wrap
 implicit none
 integer :: i, k
 character(len=240) :: mklname, fchname

 i = iargc()
 if(i /= 1) then
  write(6,'(/,A)') ' ERROR in program mkl2mrcc: wrong command line argument!'
  write(6,'(A,/)') ' Example: mkl2mrcc h2o.mkl'
  stop
 end if

 mklname = ' '
 call getarg(1, mklname)
 call require_file_exist(mklname)

 fchname = ' '
 call find_specified_suffix(mklname, '.mkl', i)
 call get_a_random_int(k)
 write(fchname,'(A,I0,A)') mklname(1:i-1)//'_',k,'.fch'

 call mkl2fch_wrap(mklname, fchname)
 call fch2mrcc_wrap(fchname, 0, .false.)
 call delete_file(fchname)
end program mkl2mrcc

