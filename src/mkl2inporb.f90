! written by jxzou at 20230728: a wrapper of mkl2fch and fch2inporb for ORCA ->
! (Open)Molcas

program mkl2inporb
 use util_wrapper, only: mkl2fch_wrap, fch2inporb_wrap
 implicit none
 integer :: i, k
 character(len=240) :: mklname, fchname, inpname

 i = iargc()
 if(i<1 .or. i>2) then
  write(6,'(/,A)') 'ERROR in program mkl2inporb: wrong command line argument!'
  write(6,'(A)')   'Example 1: mkl2inporb a.mkl'
  write(6,'(A,/)') 'Example 2: mkl2inporb a.mkl b.input'
  stop
 end if

 mklname = ' '
 call getarg(1, mklname)
 call require_file_exist(mklname)

 inpname = ' '
 if(i == 2) then
  call getarg(2, inpname)
 else
  call find_specified_suffix(mklname, '.mkl', i)
  inpname = mklname(1:i-1)//'.input'
 end if

 fchname = ' '
 call find_specified_suffix(mklname, '.mkl', i)
 call get_a_random_int(k)
 write(fchname,'(A,I0,A)') mklname(1:i-1)//'_',k,'.fch'

 call mkl2fch_wrap(mklname, fchname)
 call fch2inporb_wrap(fchname, .false., inpname)
 call delete_file(TRIM(fchname))
end program mkl2inporb

