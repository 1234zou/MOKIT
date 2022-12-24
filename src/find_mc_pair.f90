! written by jxzou at 2018: print 
! updated by jxzou at 20190119: print zero information when no pair is valid
! updated by jxzou at 20191224: add hint information

program main
 implicit none
 integer :: i
 character(len=240) gmsname

 gmsname = ' '
 i = iargc()
 if(i /= 1) then
  write(6,'(A)') 'ERROR in subroutine find_mc_pair: wrong command line parameter!'
  write(6,'(A)') 'Example: ./find_mc_pair a.gms'
  stop
 end if

 call getarg(1,gmsname)
 call require_file_exist(gmsname)
 call find_mc_pair(gmsname)
end program main

subroutine find_mc_pair(gmsname)
 implicit none
 integer :: i, j, k, npair, tmp_pair(2)
 integer, parameter :: fid = 11
 integer, allocatable :: pair(:,:)
 real(kind=8) tmp_coeff(2), tmp_coeff1(2), tmp_coeff2(2)
 real(kind=8), parameter :: threshold = -1.0d-5
 real(kind=8), allocatable :: pair_coeff(:,:)
 character(len=240) buffer
 character(len=240), intent(in) :: gmsname

 buffer =  ' '
 open(unit=fid,file=TRIM(gmsname),status='old',position='rewind')

 ! find the number of pairs
 do while(.true.)
  read(fid,'(A)') buffer
  i = index(buffer,'ROHF-GVB INPUT')
  if(i /= 0) exit
 end do
 if(i == 0) then
  write(6,'(A)') "ERROR in subroutine find_mc_pair: no 'ROHF-GVB INPUT' string found!"
  write(6,'(A)') 'The .gms file '//TRIM(gmsname)//' is not complete!'
  close(fid)
  return
 end if

 ! skip 3 lines
 do i = 1, 3, 1
  read(fid,'(A)') buffer
 end do

 read(fid,'(A)') buffer
 i = index(buffer,'=')
 buffer(1:i) = ' '
 buffer = ADJUSTL(buffer)
 read(buffer,*) npair
 ! npair found

 allocate(pair_coeff(2,npair), pair(2,npair))
 pair_coeff = 0.0d0
 pair = 0

 ! find the pair coefficients
 do while(.true.)
  read(fid,'(A)') buffer
  i = index(buffer,'PAIR INFORMATION')
  if(i /= 0) exit
 end do
 if(i == 0) then
  write(6,'(A)') 'ERROR in subroutine find_mc_pair: no PAIR INFORMATION found!'
  write(6,'(A)') 'The .gms file '//TRIM(gmsname)//' is not complete!'
  close(fid)
  return
 end if

 ! skip 3 lines
 do i = 1, 3, 1
  read(fid,'(A)') buffer
 end do

 do i = 1, npair, 1
  read(fid,*) j, pair(1,i), pair(2,i), pair_coeff(1,i), pair_coeff(2,i)
 end do
 ! pair coefficients read done
 close(fid)

 ! check the pair coefficients
 do i = 1, npair, 1
  tmp_coeff = pair_coeff(:,i)
  if(DABS(tmp_coeff(1)) < DABS(tmp_coeff(2))) then
   pair_coeff(1,i) = DABS(tmp_coeff(2))
   pair_coeff(2,i) = -DABS(tmp_coeff(1))
   k = pair(2,i)
   pair(2,i) = pair(1,i)
   pair(1,i) = k
  end if
 end do
 ! check done

 ! sort the pair coefficients
 tmp_coeff = 0.0d0
 tmp_coeff1 = 0.0d0
 tmp_coeff2 = 0.0d0
 do i = 1, npair-1, 1
  tmp_coeff1 = pair_coeff(:,i)
  do j = i+1, npair, 1
   tmp_coeff2 = pair_coeff(:,j)
   if(tmp_coeff2(2) < tmp_coeff1(2)) then
    tmp_coeff = pair_coeff(:,j)
    pair_coeff(:,j) = pair_coeff(:,i)
    pair_coeff(:,i) = tmp_coeff
    tmp_coeff1 = pair_coeff(:,i) ! update tmp_coeff1
    tmp_pair = pair(:,i)
    pair(:,i) = pair(:,j)
    pair(:,j) = tmp_pair
   end if
  end do
 end do
 ! sort done

 ! output the number of required pair coefficients
 k = count(pair_coeff(2,:) < threshold)
 write(6,'(/)',advance='no')
 if(k > 0) then
  write(6,'(6X,A)') 'ORBITAL   CI COEFFICIENTS'
  write(6,'(A)')    ' PAIR 1   2     ORB 1     ORB 2'
  do i = 1, k, 1
   write(6,'(I3,2I4,2X,2F10.6)') i, pair(1:2,i), pair_coeff(1:2,i)
  end do
 else
  write(6,'(A)') 'No pair '
 end if
 ! output done

 deallocate(pair_coeff, pair)
end subroutine find_mc_pair

