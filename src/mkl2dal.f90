! written by jxzou at 20230728: a wrapper of mkl2fch and fch2dal for ORCA ->
! Dalton

program mkl2dal
 use util_wrapper, only: mkl2fch_wrap
 implicit none
 integer :: i, k, SYSTEM, RENAME
 character(len=240) :: mklname, fchname, inpname0, inpname, molfile0, molfile

 i = iargc()
 if(.not. (i==1 .or. i==2)) then
  write(6,'(/,A)') 'ERROR in program mkl2dal: wrong command line argument!'
  write(6,'(A)')   'Example 1: mkl2dal a.mkl'
  write(6,'(A,/)') 'Example 2: mkl2dal a.mkl b.dal'
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
  inpname = mklname(1:i-1)//'.dal'
 end if

 fchname = ' '
 call find_specified_suffix(mklname, '.mkl', i)
 call get_a_random_int(k)
 write(fchname,'(A,I0,A)') mklname(1:i-1)//'_',k,'.fch'

 call mkl2fch_wrap(mklname, fchname, .false.)

 i = SYSTEM('fch2dal '//TRIM(fchname))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in program mkl2dal: failed to call utility fch2dal.'
  write(6,'(A)') 'mkl2dal is a wrapper of mkl2fch and fch2dal, so fch2dal is re&
                 &quired to be'
  write(6,'(A)') 'called successfully.'
  stop
 end if

 call delete_file(fchname)
 call find_specified_suffix(fchname, '.fch', i)
 inpname0 = fchname(1:i-1)//'.dal'
 molfile0 = fchname(1:i-1)//'.mol'

 call find_specified_suffix(inpname, '.dal', i)
 molfile = inpname(1:i-1)//'.mol'

 i = RENAME(TRIM(inpname0), TRIM(inpname))
 i = RENAME(TRIM(molfile0), TRIM(molfile))
end program mkl2dal

