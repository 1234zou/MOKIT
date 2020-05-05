! written by jxzou at 20191021: sort (part of) MOs in descending order of the pair coefficients
!                               of the 1st natural orbital in each pair

! Note: the input file must be in GAMESS format, i.e., .dat or .inp

program main
 implicit none
 integer :: i
 integer :: nbf, nif, nocc, nopen, npair
 integer, parameter :: iout = 6
 character(len=10) :: buf
 character(len=240) :: datname

 nocc = 0; nopen = 0; npair = 0
 buf = ' '
 datname = ' '
 i = iargc()
 if(i /= 6) then
  write(iout,'(/,A)') 'ERROR in subroutine gvb_sort_pairs: wrong command line arguments!'
  write(iout,'(/,A)') 'Format: gvb_sort_pairs a.dat nbf nif nocc nopen npair'
  write(iout,'(/,A)') 'Example: gvb_sort_pairs a.dat 548 548 40 1 71'
  stop
 end if

 call getarg(1, datname)
 call getarg(2, buf)
 read(buf,*) nbf
 call getarg(3, buf)
 read(buf,*) nif
 call getarg(4, buf)
 read(buf,*) nocc
 call getarg(5, buf)
 read(buf,*) nopen
 call getarg(6, buf)
 read(buf,*) npair

 if(nif > nbf) then
  write(iout,'(/,A)') 'ERROR in subroutine gvb_sort_pairs: nif>nbf. Impossible!'
  stop
 end if
 if(npair == 0) then
  write(iout,'(/,A)') 'ERROR in subroutine gvb_sort_pairs: npair=0. Not allowed!'
  stop
 end if

 call gvb_sort_pairs(datname, nbf, nif, nocc, nopen, npair)
 stop
end program main

subroutine gvb_sort_pairs(datname, nbf, nif, nocc, nopen, npair)
 implicit none
 integer :: i, j, k, m, nleft, nline
 integer :: datid, fid
 integer, intent(in) :: nbf, nif, nocc, nopen, npair
 integer, parameter :: iout = 6
 real(kind=8), allocatable :: pair_coeff(:,:), mo_coeff(:,:)
 real(kind=8), allocatable :: tmp_coeff1(:), tmp_coeff2(:)
 character(len=5) :: str1
 character(len=30) :: str2
 character(len=240) :: buf, fname
 character(len=240), intent(in) :: datname

 str1 = ' '
 str2 = ' '
 buf = ' '
 fname = ' '

 if(nocc+nopen+2*npair > nif) then
  write(iout,'(A)') 'ERROR in subroutine gvb_sort_pairs: (nocc+nopen+2*npair)>nif!'
  write(iout,'(A,4I5)') 'nocc, nopen, npair, nif=', nocc, nopen, npair, nif
  stop
 end if

 open(newunit=datid,file=TRIM(datname),status='old',position='rewind')

 ! find pair coefficients
 do while(.true.)
  read(datid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(index(buf,'CICOEF(') /= 0) exit
 end do
 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine gvb_sort_pairs: no 'CICOEF(' found!"
  write(iout,'(A)') 'The input file '//TRIM(datname)//' is not complete!'
  close(datid)
  stop
 end if

 ! read pair coefficients
 BACKSPACE(datid)
 allocate(pair_coeff(2,npair), source=0.0d0)
 do i = 1, npair, 1
  read(datid,'(A)') buf
  k = index(buf,'=')
  read(buf(k+1:),*) pair_coeff(1,i)
  k = index(buf,',')
  read(buf(k+1:),*) pair_coeff(2,i)
 end do
 ! pair coefficients read done

 ! read MO coefficients
 do while(.true.)
  read(datid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:5) == '$VEC') exit
 end do
 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine gvb_sort_pairs: no '$VEC' found in&
                  & file "//TRIM(datname)//'!'
  close(datid)
  stop
 end if

 allocate(mo_coeff(nbf,nif), source=0.0d0)
 nline = nbf/5
 nleft = nbf - nline*5
 do i = 1, nif, 1
  k = 1
  do j = 1, nline, 1
   read(datid,'(A)') buf
   buf = buf(6:)
   read(buf,'(5ES15.8)') mo_coeff(k:k+4,i)
   k = k + 5
  end do ! for j
  if(nleft > 0) then
   read(datid,'(A)') buf
   buf = buf(6:)
   str1 = ' '
   str2 = ' '
   write(str1,'(I5)') nleft
   str2 = '('//TRIM(ADJUSTL(str1))//'ES15.8)'
   read(buf,TRIM(str2)) mo_coeff(k:nbf,i)
  end if
 end do ! for i

 ! read MOs done
 close(datid)

 ! check the pair coefficients: if ABS(2nd NO) > ABS(1st NO), swap them
 do i = 1, npair, 1
  tmp_coeff1 = DABS(pair_coeff(:,i)) ! auto-allocation
  if(tmp_coeff1(1) < tmp_coeff1(2)) then
   pair_coeff(1,i) = tmp_coeff1(2)
   pair_coeff(2,i) = -tmp_coeff1(1)
   k = nocc + nopen + 2*i - 1
   tmp_coeff2 = mo_coeff(:,k) ! auto-allocation
   mo_coeff(:,k) = mo_coeff(:,k+1)
   mo_coeff(:,k+1) = tmp_coeff2
  end if
 end do
 ! check done

 ! sort the pair coefficients (MOs will be permuted accordingly)
 do i = 1, npair-1, 1
  tmp_coeff1 = pair_coeff(:,i)
  k = nocc + nopen + 2*i - 1
  do j = i+1, npair, 1
   if(pair_coeff(1,j) > tmp_coeff1(1)) then
    tmp_coeff1 = pair_coeff(:,j)
    pair_coeff(:,j) = pair_coeff(:,i)
    pair_coeff(:,i) = tmp_coeff1
    m = nocc + nopen + 2*j - 1
    tmp_coeff1 = mo_coeff(:,k)   ! auto-reallocation
    tmp_coeff2 = mo_coeff(:,k+1) ! auto-reallocation
    mo_coeff(:,k:k+1) = mo_coeff(:,m:m+1)
    mo_coeff(:,m) = tmp_coeff1
    mo_coeff(:,m+1) = tmp_coeff2
    tmp_coeff1 = pair_coeff(:,i) ! remember to update tmp_coeff1
   end if
  end do ! for j
 end do ! for i
 ! sort done

 if(allocated(tmp_coeff1)) deallocate(tmp_coeff1)
 if(allocated(tmp_coeff2)) deallocate(tmp_coeff2)

 open(newunit=datid,file=TRIM(datname),status='old',position='rewind')
 i = index(datname,'.dat')
 if(i == 0) i = index(datname,'.inp')
 fname = datname(1:i-1)//'_sort.dat'
 open(newunit=fid,file=TRIM(fname),status='replace')
 ! copy information before the 'CICOEF('
 do while(.true.)
  read(datid,'(A)') buf
  if(index(buf, 'CICOEF(') /= 0) exit
  write(fid,'(A)') TRIM(buf)
 end do

 ! print sorted pair coefficients into a new file
 if(buf(2:5) == '$SCF') then
  write(fid,'(1X,2(A,F12.8))') '$SCF    CICOEF(  1)=', pair_coeff(1,1), ',', pair_coeff(2,1)
 else
  write(fid,'(9X,2(A,F12.8))') 'CICOEF(  1)=', pair_coeff(1,1), ',', pair_coeff(2,1)
 end if
 do i = 2, npair, 1
  write(fid,'(9X,A,I3,2(A,F12.8))') 'CICOEF(',2*i-1,')=',pair_coeff(1,i),',',pair_coeff(2,i)
 end do
 write(fid,'(1X,A)') '$END'
 ! print done

 ! skip the pair coefficients in .dat file
 do while(.true.)
  read(datid,'(A)') buf
  if(index(buf, '$END') /= 0) exit
 end do

 ! copy any information between '$END' and 'VEC'
 do while(.true.)
  read(datid,'(A)') buf
  if(index(buf, '$VEC') /= 0) exit
  write(fid,'(A)') TRIM(buf)
 end do
 ! copy done

 ! print the sorted MOs into a new file
 write(fid,'(1X,A)') '$VEC'
 nline = nbf/5
 nleft = nbf - 5*nline
 do i = 1, nif, 1
  k = MOD(i,100)
  do j = 1, nline, 1
   write(fid,'(I2,I3,5ES15.8)') k, MOD(j,1000), mo_coeff(5*j-4:5*j,i)
  end do
  if(nleft > 0) then
   write(fid,'(I2,I3,5ES15.8)') k, MOD(j,1000), mo_coeff(5*j-4:nbf,i)
  end if
 end do
 write(fid,'(1X,A)') '$END'
 ! print done

 ! skip the '$VEC' section in .dat file
 do while(.true.)
  read(datid,'(A)') buf
  if(index(buf, '$END') /= 0) exit
 end do

 ! copy the rest in .dat file
 do while(.true.)
  read(datid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid,'(A)') TRIM(buf)
 end do
 close(datid)
 close(fid)

 deallocate(pair_coeff, mo_coeff)
 return
end subroutine gvb_sort_pairs

