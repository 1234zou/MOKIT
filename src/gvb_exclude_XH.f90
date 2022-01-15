! written by jxzou at 20180911
! updated by jxzou at 20190109: print X-H pairs
! moved into MOKIT by jxzou at 20210914

! This is a program/subroutine used to exclude all X-H pairs out of the GVB active space,
!  except those multiconfigurational pairs.

program main
 implicit none
 integer :: i
 integer, parameter :: iout = 6
 character(len=240) :: datname, gmsname

 i = IARGC()
 if(i /= 2) then
  write(iout,'(/,A)') ' ERROR in subroutine gvb_exclude_XH: wrong command line prameter!'
  write(iout,'(A,/)') ' Example: gvb_exclude_XH a.dat a.gms'
  stop
 end if

 call GETARG(1,datname)
 call GETARG(2,gmsname)
 call gvb_exclude_XH(datname, gmsname)
 stop
end program main

! exclude all X-H pairs out of the GVB active space (except those pairs with
! multiconfigurational character)
subroutine gvb_exclude_XH(datname, gmsname)
 implicit none
 integer :: i, j, k, m, itmp, itmp1(1), itmp2(1)
 integer :: natom, ncore, nopen, npair, npair2, nbf, nif
 integer, parameter :: iout = 6
 integer, allocatable :: bf2atom(:), order(:)
 real(kind=8), allocatable :: coeff(:,:), coeff2(:,:)
 real(kind=8), allocatable :: tmp_coeff(:)
 real(kind=8), allocatable :: ci_coeff(:,:), ci_coeff2(:,:)
 real(kind=8) :: rtmp1, rtmp2, tmp_ci_coeff(2)
 real(kind=8), parameter :: Pt1 = 0.1d0
 character(len=2), allocatable :: elem(:)
 character(len=240), intent(in) :: datname, gmsname
 character(len=240) :: newdat, inpname
 logical, allocatable :: xhbond(:)
 logical :: alive

 ! get ncore, nopen and npair
 call read_npair_from_gms(gmsname, ncore, nopen, npair)

 ! get nbf and nif
 call read_nbf_and_nif_from_gms(gmsname, nbf, nif)

 ! get GVB CI coefficients from .dat file
 allocate(ci_coeff(2,npair))
 ci_coeff = 0d0
 call chk_ci_coeff_in_dat(datname, alive)
 if(alive) then
  call read_ci_coeff_from_dat(datname, npair, ci_coeff)
 else
  call read_ci_coeff_from_gms(gmsname, npair, ci_coeff)
 end if

 ! read MO from .dat file
 allocate(coeff(nbf,nif))
 call read_mo_from_dat(datname, nbf, nif, coeff)

 ! reverse C2 with C1 if C2 > C1
 allocate(order(2*npair))
 forall(i = 1:2*npair) order(i) = i

 do i = 1, npair, 1
  rtmp1 = DABS(ci_coeff(1,i))
  rtmp2 = DABS(ci_coeff(2,i))
  if(rtmp1 < rtmp2) then
   ci_coeff(1,i) = rtmp2
   ci_coeff(2,i) = -rtmp1
   itmp = order(2*i-1)
   order(2*i-1) = order(2*i)
   order(2*i) = itmp
  else
   if(ci_coeff(1,i) < 0d0) then
    ci_coeff(1,i) = rtmp1
    ci_coeff(2,i) = -rtmp2
   end if
  end if
 end do ! for i

 ! sort pairs by C2
 do i = 1, npair-1, 1
  rtmp1 = ci_coeff(2,i)

  do j = i+1, npair, 1
   rtmp2 = ci_coeff(2,j)

   if(rtmp1 > rtmp2) then
    tmp_ci_coeff = ci_coeff(:,i)
    ci_coeff(:,i) = ci_coeff(:,j)
    ci_coeff(:,j) = tmp_ci_coeff
    itmp = order(2*i-1)
    order(2*i-1) = order(2*j-1)
    order(2*j-1) = itmp
    itmp = order(2*i)
    order(2*i) = order(2*j)
    order(2*j) = itmp
    rtmp1 = rtmp2 ! update rtmp1
   end if
  end do ! for j
 end do ! for i

 ! update coeff to coeff2
 allocate(coeff2(nbf,2*npair))
 coeff2 = 0d0
 m = ncore + nopen
 forall(i = 1:2*npair) coeff2(:,i) = coeff(:,m+order(i))

 ! get natom
 call read_natom_from_gms(gmsname, natom)

 ! get elem and bf2atom
 allocate(elem(natom), bf2atom(nbf))
 elem = ' '
 bf2atom = 0
 call read_elem_and_bf2atom_from_gms(gmsname, natom, elem, nbf, bf2atom)

 ! excluding all X-H pairs from the GVB active space
 allocate(xhbond(npair), tmp_coeff(nbf))
 npair2 = 0
 xhbond = .false.
 tmp_coeff = 0d0

 do i = 1, npair, 1
  if(ci_coeff(2,i) < -Pt1) cycle ! skip multiconfigurational pair(s)
  tmp_coeff = coeff2(:,2*i-1) ! only consider the bonding orbital
  forall(j = 1:nbf) tmp_coeff(j) = DABS(tmp_coeff(j))
  itmp1 = MAXLOC(tmp_coeff)
  j = itmp1(1)
  tmp_coeff(j) = 0d0
  itmp2 = MAXLOC(tmp_coeff)
  k = itmp2(1)
  if(elem(bf2atom(j))=='H ' .or. elem(bf2atom(k))=='H ') then
   xhbond(i) = .true.
   npair2 = npair2 + 1
  end if
 end do ! for i

 ! now npair2 is the number of X-H pairs
 deallocate(elem, bf2atom, tmp_coeff)

 write(iout,'(A)') 'Report from subroutine gvb_exclude_XH:'
 if(npair2 == 0) then
  write(iout,'(A,/)') ' No X-H bond found in file '//TRIM(gmsname)
  return
 end if
 write(iout,'(I0,A,/)') npair2, ' X-H bond(s) found in file '//TRIM(gmsname)

 ! CI coefficients of non X-H pairs will be reserved in ci_coeff2
 allocate(ci_coeff2(2,npair-npair2), source=0d0)
 k = 0
 do i = 1, npair, 1
  if(xhbond(i)) cycle
  k = k + 1
  ci_coeff2(:,k) = ci_coeff(:,i)
 end do ! for i
 deallocate(ci_coeff)

 ! reset order
 forall(i = 1:2*npair) order(i) = i

 k = 0
 ! for X-H bonds
 do i = 1, npair, 1
  if(.not. xhbond(i)) cycle
  k = k + 1
  order(k) = 2*i - 1
  order(2*npair-k+1) = 2*i
 end do ! for i

 ! for non X-H bonds
 do i = 1, npair, 1
  if(xhbond(i)) cycle
  k = k + 1
  order(k) = 2*i - 1
  order(k+1) = 2*i
  k = k + 1
 end do ! for i

 deallocate(xhbond)
 ! update coeff2 to coeff
 forall(i = 1:2*npair) coeff(:,m+i) = coeff2(:,order(i))
 ! now is doubly_occupied -> singly_occupied -> X-H_doubly_occupied -> pair
 deallocate(coeff2, order)

 ! move the X-H doubly occupied orbitals into the doubly occupied group
 if(nopen > 0) then
  allocate(coeff2(nbf,nopen))
  coeff2 = coeff(:,ncore+1:m)

  do i = 1, npair2, 1
   coeff(:,ncore+i) = coeff(:,m+i)
  end do ! for i

  k = ncore + npair2
  coeff(:,k+1:k+nopen) = coeff2
  deallocate(coeff2)
 end if
 ! now is doubly_occupied -> X-H_doubly_occupied -> singly_occupied -> pair

 ! reverse the order of pair orbitals, since important pairs are often located
 ! close to core orbitals, we usually wish them to be located near HOMO, LUMO
 ! labels
 j = ncore + npair2 + nopen
 k = npair - npair2
 allocate(order(k), source=0)
 allocate(coeff2(nbf,k), source=0d0)

 forall(i = 1:k) order(i) = j + 2*(k-i) + 1
 forall(i = 1:k) coeff2(:,i) = coeff(:,order(i))
 forall(i = 1:k) coeff(:,j+2*i-1) = coeff2(:,i)

 forall(i = 1:k) order(i) = j + 2*(k-i) + 2
 forall(i = 1:k) coeff2(:,i) = coeff(:,order(i))
 forall(i = 1:k) coeff(:,j+2*i) = coeff2(:,i)
 deallocate(order, coeff2)

 ! print MO into .dat file
 call print_mo_into_dat(datname, nbf, nif, coeff, .false.)
 deallocate(coeff, ci_coeff2)
 ! Note that a new .dat file will be generated

 i = index(datname,'.dat')
 if(i == 0) i = index(datname,'.inp')
 newdat = datname(1:i-1)//'_new.dat'
 write(inpname,'(A)') datname(1:i-1)//'_XH.inp'
 call create_gvb_inp_from_dat_and_gms(newdat, gmsname, inpname, ncore+npair2,&
                                      npair-npair2, nopen)
 return
end subroutine gvb_exclude_XH

! read ncore, nopen and npair from a GAMESS .gms file
subroutine read_npair_from_gms(gmsname, ncore, nopen, npair)
 implicit none
 integer :: i, fid
 integer, intent(out) :: ncore, nopen, npair
 integer, parameter :: iout = 6
 character(len=240) :: buf
 character(len=240), intent(in) :: gmsname

 ncore = 0; nopen = 0; npair = 0
 buf = ' '

 open(newunit=fid,file=TRIM(gmsname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)', iostat=i) buf
  if(i /= 0) exit
  if(buf(34:41) == 'NCO    =') exit
 end do

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine read_npair_from_gms: no 'NCO    ='&
                   & found"
  write(iout,'(A)') ' in file '//TRIM(gmsname)//'.'
  close(fid)
  stop
 end if

 read(buf(42:),*) ncore
 read(fid,'(A)') buf
 read(buf(19:),*) npair
 read(buf(42:),*) nopen
 close(fid)
 return
end subroutine read_npair_from_gms

! read nbf and nif from a GAMESS output file
subroutine read_nbf_and_nif_from_gms(gmsname, nbf, nif)
 implicit none
 integer :: i, fid
 integer, intent(out) :: nbf, nif
 integer, parameter :: iout = 6
 character(len=240) :: buf
 character(len=240), intent(in) :: gmsname

 nbf = 0; nif = 0
 open(newunit=fid,file=TRIM(gmsname),status='old',position='rewind')

 ! find nbf
 do while(.true.)
  read(fid,'(A)', iostat=i) buf
  if(i /= 0) exit
  if(buf(2:22) == 'NUMBER OF CARTESIAN G') exit
 end do

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine read_nbf_and_nif_from_gms: no 'NUMBER&
                   & OF CARTESIAN G'"
  write(iout,'(A)') ' found in file '//TRIM(gmsname)//'.'
  close(fid)
  stop
 end if

 i = index(buf, '=')
 read(buf(i+1:),*) nbf
 ! This is actually the number of Cartesian basis functions, but we just need
 ! this number, because MOs in .dat file are expressed in Cartesian basis
 ! functions, no matter you use spherical harmonic or Cartesian basis functions

 ! find nif
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(11:23) == 'GUESS OPTIONS') exit
 end do

 read(fid,'(A)') buf   ! skip one line
 read(fid,'(A)') buf
 i = index(buf, 'NORB  =')
 read(buf(i+7:),*) nif
 close(fid)
 return
end subroutine read_nbf_and_nif_from_gms

! check if there are any GVB CI coefficients in a GAMESS .dat or .inp file
subroutine chk_ci_coeff_in_dat(fname, alive)
 implicit none
 integer :: i, j, fid
 character(len=240) :: buf
 character(len=240), intent(in) :: fname
 logical, intent(out) :: alive

 buf = ' '
 alive = .true.

 open(newunit=fid,file=TRIM(fname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  j = index(buf,'CICOEF(')
  if(j == 0) j = index(buf,'cicoef(')
  if(j /= 0) exit
 end do
 close(fid)

 if(i /= 0) alive = .false.
 return
end subroutine chk_ci_coeff_in_dat

! read CI coefficients from a GAMESS .dat or .inp file
subroutine read_ci_coeff_from_dat(fname, npair, coeff)
 implicit none
 integer :: i, j, k, fid
 integer, intent(in) :: npair
 integer, parameter :: iout = 6
 character(len=240) :: buf
 character(len=240), intent(in) :: fname
 real(kind=8), intent(out) :: coeff(2,npair)

 buf = ' '; coeff = 0d0
 open(newunit=fid,file=TRIM(fname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  j = index(buf,'CICOEF(')
  if(j == 0) j = index(buf,'cicoef(')
  if(j /= 0) exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine read_ci_coeff_from_dat: no GVB CI&
                    & coefficients'
  write(iout,'(A)') ' found in file '//TRIM(fname)//'!'
  close(fid)
  stop
 end if

 BACKSPACE(fid)
 do i = 1, npair, 1
  read(fid,'(A)') buf
  j = index(buf,'=')
  k = index(buf,',')
  read(buf(j+1:k-1),*) coeff(1,i)
  read(buf(k+1:),*) coeff(2,i)
 end do ! for i

 close(fid)
 return
end subroutine read_ci_coeff_from_dat

! read CI coefficients from a GAMESS .gms file
! Note: if there exist multiple sets of CI coefficients in the file,
!       only the last set will be read
subroutine read_ci_coeff_from_gms(fname, npair, coeff)
 implicit none
 integer :: i, j, fid
 integer, intent(in) :: npair
 integer, parameter :: iout = 6
 character(len=240) :: buf
 character(len=240), intent(in) :: fname
 real(kind=8), intent(out) :: coeff(2,npair)

 buf = ' '
 coeff = 0d0

 open(newunit=fid,file=TRIM(fname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid,iostat=i)
  BACKSPACE(fid,iostat=i)
  read(fid,'(A)') buf
  if(i/=0 .or. buf(7:18)=='ORBITAL   CI') exit
 end do

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine read_ci_coeff_from_gms: no GVB CI&
                   & coefficients'
  write(iout,'(A)') ' found in file '//TRIM(fname)//'!'
  close(fid)
  stop
 end if

 read(fid,'(A)') buf   ! skip one line

 do i = 1, npair, 1
  read(fid,'(A)') buf
  read(buf,*) j, j, j, coeff(1,i), coeff(2,i)
 end do ! for i

 close(fid)
 return
end subroutine read_ci_coeff_from_gms

! read natom from GAMESS .gms file
subroutine read_natom_from_gms(gmsfile, natom)
 implicit none
 integer :: i, fid
 integer, intent(out) :: natom
 integer, parameter :: iout = 6
 character(len=240) :: buf
 character(len=240), intent(in) :: gmsfile

 buf = ' '
 natom = 0
 open(newunit=fid,file=TRIM(gmsfile),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:22) == 'TOTAL NUMBER OF ATOMS') exit
 end do
 close(fid)

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine read_natom_from_gms: no 'TOTAL NUMBER&
                   & OF ATOMS'"
  write(iout,'(A)') 'found in file '//TRIM(gmsfile)//'.'
  stop
 end if

 i = index(buf,'=')
 read(buf(i+1:),*) natom
 return
end subroutine read_natom_from_gms

! read array elem and bf2atom from .gms file
! bf2atom: the map from basis functions to atoms
subroutine read_elem_and_bf2atom_from_gms(gmsfile, natom, elem, nbf, bf2atom)
 implicit none
 integer :: i, j, nelem, fid
 integer, intent(in) :: natom, nbf
 integer, intent(out) :: bf2atom(nbf)
 integer, parameter :: iout = 6
 character(len=2) :: tmp_elem
 character(len=240) :: buf
 character(len=2), intent(out) :: elem(natom)
 character(len=240), intent(in) :: gmsfile

 elem = ' '; buf = ' '; bf2atom = 0
 open(newunit=fid,file=TRIM(gmsfile),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(11:22) == 'EIGENVECTORS') exit
 end do

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine read_bf2atom_from_gms: no &
                    &'EIGENVECTORS' found"
  write(iout,'(A)') 'file '//TRIM(gmsfile)//'.'
  close(fid)
  stop
 end if

 do i = 1, 5, 1   ! skip 5 lines
  read(fid,'(A)') buf
 end do

 read(fid,'(A)') buf
 read(buf,*) j, elem(1), bf2atom(1)
 tmp_elem = ' '
 nelem = 1

 do i = 2, nbf, 1
  read(fid,'(A)') buf
  read(buf,*) j, tmp_elem, bf2atom(i)
  if(bf2atom(i) /= bf2atom(i-1)) then
   nelem = nelem + 1
   elem(nelem) = tmp_elem
  end if
 end do ! for i

 close(fid)
 forall(i = 1:natom) elem(i) = ADJUSTL(elem(i))
 return
end subroutine read_elem_and_bf2atom_from_gms

! print MOs into .dat file
! if replace is .true., the MOs in the original file will be replaced
! if replace is .false., a new *_new.dat file will be generated
subroutine print_mo_into_dat(datname, nbf, nif, coeff, replace)
 implicit none
 integer :: i, j, k, nline, nleft, fid1, fid2, RENAME
 integer, intent(in) :: nbf, nif
 real(kind=8), intent(in) :: coeff(nbf,nif)
 character(len=240), intent(in) :: datname
 character(len=240) :: newdat, buf
 logical, intent(in) :: replace

 buf = ' '; newdat = ' '
 i = index(datname,'.dat',.true.)
 if(i == 0) i = index(datname,'.inp',.true.)
 newdat = datname(1:i-1)//'_new.dat'

 open(newunit=fid1,file=TRIM(datname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(newdat),status='replace')
 do while(.true.)
  read(fid1,'(A)') buf
  write(fid2,'(A)') TRIM(buf)
  call upper(buf(3:5))
  if(buf(2:5)=='$VEC') exit
 end do ! for while

 ! print MOs
 nline = nbf/5
 nleft = nbf - 5*nline

 do i = 1, nif, 1
  k = MOD(i,100)
  do j = 1, nline, 1
   write(fid2,'(I2,I3,5ES15.8)') k, MOD(j,1000), coeff(5*j-4:5*j,i)
  end do ! for j
  if(nleft > 0) then
   write(fid2,'(I2,I3,5ES15.8)') k, MOD(j,1000), coeff(5*j-4:nbf,i)
  end if
 end do ! for i
 write(fid2,'(A)') ' $END'
 ! print MOs done

 ! skip the MOs in datname
 do while(.true.)
  read(fid1,'(A)') buf
  call upper(buf(3:5))
  if(buf(2:5) == '$END') exit
 end do

 ! copy remaining contens
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while
 close(fid2)
 ! copy done

 if(replace) then
  close(fid1,status='delete')
  i = RENAME(TRIM(newdat), TRIM(datname))
 else
  close(fid1)
 end if

 return
end subroutine print_mo_into_dat

! print GVB CI coefficients into .dat file
! if replace is .true., the CI coefficients in the original file will be
! replaced
! if replace is .false., a new *_new.dat file will be generated
subroutine print_ci_coeff_into_dat(datname, npair, coeff, replace)
 implicit none
 integer :: i, fid1, fid2, RENAME
 integer, intent(in) :: npair
 real(kind=8), intent(in) :: coeff(2,npair)
 character(len=240), intent(in) :: datname
 character(len=240) :: newdat, buf
 logical, intent(in) :: replace
 logical :: alive

 buf = ' '; newdat = ' '
 call chk_ci_coeff_in_dat(datname, alive)

 i = index(datname,'.dat',.true.)
 if(i == 0) i = index(datname,'.inp',.true.)
 newdat = datname(1:i-1)//'.tmp'

 open(newunit=fid1,file=TRIM(datname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(newdat),status='replace')

 if(alive) then
  do while(.true.)
   read(fid1,'(A)') buf
   i = index(buf,'CICOEF')
   if(i == 0) i = index(buf,'cicoef')
   if(i /= 0) exit
   write(fid2,'(A)') TRIM(buf)
  end do ! for while

  if(buf(2:5) == '$SCF') write(fid2,'(A)') ' $SCF'
 else
  do while(.true.)
   read(fid1,'(A)') buf
   if(buf(2:5) == '$VEC') exit
   write(fid2,'(A)') TRIM(buf)
  end do ! for while

  write(fid2,'(A)') ' $SCF'
 end if

 ! print GVB CI coefficients
 do i = 1, npair, 1
  write(fid2,'(2X,A7,I3,A2,E18.11,A1,E18.11)') 'CICOEF(', 2*i-1, ')=',coeff(1,i),&
                                            & ',', coeff(2,i)
 end do ! for i

 write(fid2,'(A)') ' $END'
 ! print CI coefficients done

 ! skip the GVB CI coefficients in datname
 if(alive) then
  do while(.true.)
   read(fid1,'(A)') buf
   if(index(buf,'$END')>0 .or. index(buf,'$End')>0) exit
  end do ! for while
 else
  BACKSPACE(fid1)
 end if

 ! copy remaining contents
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while
 close(fid2)
 ! copy done

 if(replace) then
  close(fid1,status='delete')
  i = RENAME(TRIM(newdat), TRIM(datname))
 else
  close(fid1)
 end if

 return
end subroutine print_ci_coeff_into_dat

! create a GAMESS GVB .inp file (inpname) according to given newdat and gmsname
subroutine create_gvb_inp_from_dat_and_gms(newdat, gmsname, inpname, ncore,&
                                           npair, nopen)
 implicit none
 integer :: i, fid1, fid2
 integer, parameter :: iout = 6
 integer, intent(in) :: ncore, npair, nopen
 character(len=240) :: buf
 character(len=240), intent(in) :: newdat, gmsname, inpname

 ! copy $CONTRL section from gmsname
 open(newunit=fid1,file=TRIM(gmsname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname),status='replace')

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(14:20) == '$CONTRL') exit
 end do ! for while

 BACKSPACE(fid1)
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(14:17) == '$SCF') exit
  write(fid2,'(A)') TRIM(buf(13:))
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine create_gvb_inp_from_dat_and_gms: no&
                   & '$SCF' found in file "//TRIM(gmsname)
  close(fid1)
  close(fid2,status='delete')
  stop
 end if

 ! write $SCF and $GUESS sections into inpname
 write(fid2,'(2(A,I0))',advance='no') ' $SCF NCO=', ncore, ' NPAIR=', npair
 if(nopen > 0) then
  write(fid2,'(A,I0,A)',advance='no') ' NSETO=',nopen,' NO(1)=1'
  do i = 2, nopen, 1
   write(fid2,'(A)',advance='no') ',1'
  end do ! for i
 end if
 write(fid2,'(A)') ' DIRSCF=.T. $END'

 ! skip the $SCF section in gmsname
 BACKSPACE(fid1)
 do while(.true.)
  read(fid1,'(A)') buf
  if(index(buf, '$END') > 0) exit
 end do ! for while

 ! copy $GUESS (and maybe other sections) from gmsname into inpname
 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(14:18) == '$DATA') exit
  write(fid2,'(A)') TRIM(buf(13:))
 end do ! for while
 close(fid1)

 ! copy the $DATA section from newdat into inpname
 open(newunit=fid1,file=TRIM(newdat),status='old',position='rewind')
 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(2:6) == '$DATA') exit
 end do ! for while
 write(fid2,'(A)') ' $DATA'

 do while(.true.)
  read(fid1,'(A)') buf
  write(fid2,'(A)') TRIM(buf)
  if(buf(2:5) == '$END') exit
 end do ! for while

 ! copy the remaining contents (except $SCF) from newdat into inpname
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit

  if(buf(2:5) == '$SCF') then
   do while(.true.)
    read(fid1,'(A)') buf
    if(index(buf,'$END') > 0) exit
   end do ! for while
   read(fid1,'(A)') buf
  end if

  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 close(fid1,status='delete')
 close(fid2)
 return
end subroutine create_gvb_inp_from_dat_and_gms

