! written by jxzou at 20230730: a wrapper of mkl2fch and fch2qchem for ORCA ->
! Q-Chem

program mkl2qchem
 use util_wrapper, only: mkl2fch_wrap, fch2qchem_wrap
 implicit none
 integer :: i, j, k
 character(len=240) :: mklname, fchname, inpname

 i = iargc()
 if(i<1 .or. i>2) then
  write(6,'(/,A)') ' ERROR in program mkl2qchem: wrong command line argument!'
  write(6,'(A)')   ' Example 1: mkl2qchem a.mkl'
  write(6,'(A,/)') ' Example 2: mkl2qchem a.mkl b.in'
  stop
 end if

 mklname = ' '
 call getarg(1, mklname)
 call require_file_exist(mklname)

 call find_specified_suffix(mklname, '.mkl', j)
 if(i == 2) then
  call getarg(2, inpname)
  call find_specified_suffix(inpname, '.in', k)
 else
  inpname = mklname(1:j-1)//'.in'
 end if
 call get_a_random_int(k)
 write(fchname,'(A,I0,A)') mklname(1:j-1)//'_',k,'.fch'

 call mkl2fch_wrap(mklname, fchname)
 call fch2qchem_wrap(fchname, 0, 0, inpname)
 call delete_file(fchname)
end program mkl2qchem

