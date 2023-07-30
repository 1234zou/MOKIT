! written by jxzou at 20230728: a wrapper of mkl2fch and fch2amo for ORCA ->
! AMESP

program mkl2amo
 use util_wrapper, only: mkl2fch_wrap
 implicit none
 integer :: i, k, SYSTEM, RENAME
 character(len=240) :: mklname, fchname, inpname0, inpname, orbfile0, orbfile

 i = iargc()
 if(.not. (i==1 .or. i==2)) then
  write(6,'(/,A)') 'ERROR in program mkl2amo: wrong command line argument!'
  write(6,'(A)')   'Example 1: mkl2amo a.mkl'
  write(6,'(A,/)') 'Example 2: mkl2amo a.mkl b.aip'
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
  inpname = mklname(1:i-1)//'.aip'
 end if

 fchname = ' '
 call find_specified_suffix(mklname, '.mkl', i)
 call get_a_random_int(k)
 write(fchname,'(A,I0,A)') mklname(1:i-1)//'_',k,'.fch'

 call mkl2fch_wrap(mklname, fchname, .false.)

 i = SYSTEM('fch2amo '//TRIM(fchname))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in program mkl2amo: failed to call utility fch2amo.'
  write(6,'(A)') 'mkl2amo is a wrapper of mkl2fch and fch2amo, so fch2amo is re&
                 &quired to be'
  write(6,'(A)') 'called successfully.'
  stop
 end if

 call delete_file(fchname)
 call find_specified_suffix(fchname, '.fch', i)
 inpname0 = fchname(1:i-1)//'.aip'
 orbfile0 = fchname(1:i-1)//'.amo'

 call find_specified_suffix(inpname, '.aip', i)
 orbfile = inpname(1:i-1)//'.amo'

 i = RENAME(TRIM(inpname0), TRIM(inpname))
 i = RENAME(TRIM(orbfile0), TRIM(orbfile))
end program mkl2amo

