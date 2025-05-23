! written by jxzou at 20221112: input file converted in to Gaussian .gjf file

program inp2gjf
 implicit none
 integer :: i, j
 character(len=240) :: inpname, gjfname

 i = iargc()
 if(i<1 .or. i>2) then
  write(6,'(/,A)') 'ERROR in program inp2gjf: wrong command line arguments!'
  write(6,'(A)') 'Example 1: inp2gjf h2o.inp'
  write(6,'(A,/)') 'Example 2: inp2gjf h2o.inp water.gjf'
  stop
 end if

 call getarg(1, inpname)
 j = INDEX(inpname, '.', back=.true.)
 if(j == 0) then
  write(6,'(/,A)') 'ERROR in program inp2gjf: no suffix found in filename '//&
                    TRIM(inpname)
  stop
 end if

 if(i == 2) then
  call getarg(2, gjfname)
 else
  gjfname = inpname(1:j-1)//'.gjf'
 end if

 call inp2gjf_cp2k(inpname, gjfname)
end program inp2gjf

! convert CP2K .inp/.restart file into Gaussian .gjf file
subroutine inp2gjf_cp2k(inpname, gjfname)
 implicit none
 integer :: i, j, k, charge, mult, fid1, fid2
 character(len=3) :: str3
 character(len=5) :: str5
 character(len=8) :: str8
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname, gjfname

 str3 = ' '; str5 = ' '; str8 = ' '
 call read_charge_and_mult_from_cp2k_inp(inpname, charge, mult)

 ! find Cartesian coordinates first
 open(newunit=fid1,file=TRIM(inpname),status='old',position='rewind')

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  j = INDEX(buf, "&")
  if(j == 0) cycle
  k = INDEX(buf(j+1:), " ")
  str5 = buf(j+1:j+k-1)
  call upper(str5)
  if(str5 == 'COORD') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine inp2gjf_cp2k: failed to find '&COORD' in&
                 & file "//TRIM(inpname)
  close(fid1)
  stop
 end if

 open(newunit=fid2,file=TRIM(gjfname),status='replace')
 write(fid2,'(A,//,A,/)') '#','Generated by inp2gjf of MOKIT'
 write(fid2,'(I0,1X,I0)') charge, mult

 do while(.true.)
  read(fid1,'(A)') buf
  i = INDEX(buf, "&")
  if(i > 0) then
   j = INDEX(buf(i+1:), " ")
   str3 = buf(i+1:i+j-1)
   call upper(str3)
   if(str3 == 'END') exit
  end if
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 rewind(fid1)
 ! now find cell parameters a,b,c

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  j = INDEX(buf, "&")
  if(j == 0) cycle
  k = INDEX(buf(j+1:), " ")
  str8 = buf(j+1:j+k-1)
  call upper(str8)
  if(str8(1:4)=='CELL' .and. str8(5:8)/='_OPT') exit
 end do ! for while

 do i = 1, 3
  read(fid1,'(A)') buf
  buf = ADJUSTL(buf)
  j = INDEX(buf, ' ')
  write(fid2,'(A)') 'Tv '//TRIM(buf(j+1:))
 end do ! for i

 close(fid1)
 write(fid2,'(/)')
 close(fid2)
end subroutine inp2gjf_cp2k

subroutine read_charge_and_mult_from_cp2k_inp(inpname, charge, mult)
 implicit none
 integer :: i, fid
 integer, intent(out) :: charge, mult
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname

 charge = 0; mult = 1 ! default value
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  call upper(buf)

  i = INDEX(buf, 'CHARGE')
  if(i > 0) read(buf(i+6:),*) charge

  i = INDEX(buf, 'MULTIPLICITY')
  if(i > 0) read(buf(i+12:),*) mult
 end do ! for while

 close(fid)
end subroutine read_charge_and_mult_from_cp2k_inp

