! written by jxzou at 20230728: a wrapper of mkl2fch and fch2inp for ORCA ->
! GAMESS

program mkl2inp
 use util_wrapper, only: mkl2fch_wrap
 implicit none
 integer :: i, k, SYSTEM, RENAME
 character(len=240) :: mklname, fchname, inpname0, inpname

 i = iargc()
 if(.not. (i==1 .or. i==2)) then
  write(6,'(/,A)') 'ERROR in program mkl2inp: wrong command line argument!'
  write(6,'(A)')   'Example 1: mkl2inp a.mkl'
  write(6,'(A,/)') 'Example 2: mkl2inp a.mkl b.inp'
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
  inpname = mklname(1:i-1)//'.inp'
 end if

 fchname = ' '
 call find_specified_suffix(mklname, '.mkl', i)
 call get_a_random_int(k)
 write(fchname,'(A,I0,A)') mklname(1:i-1)//'_',k,'.fch'

 call mkl2fch_wrap(mklname, fchname, .false.)

 i = SYSTEM('fch2inp '//TRIM(fchname))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in program mkl2inp: failed to call utility fch2inp.'
  write(6,'(A)') 'mkl2inp is a wrapper of mkl2fch and fch2inp, so fch2inp is re&
                 &quired to be'
  write(6,'(A)') 'called successfully.'
  stop
 end if

 call delete_file(fchname)
 call find_specified_suffix(fchname, '.fch', i)
 inpname0 = fchname(1:i-1)//'.inp'

 i = RENAME(TRIM(inpname0), TRIM(inpname))
end program mkl2inp

