! written by jxzou at 20220401: generate the SA-CASSCF input file

subroutine do_sa_cas()
 use mol, only: nacto, nacte
 use mr_keyword, only: hf_fch, casscf, casscf_prog, dmrgscf_prog, dryrun, excited
 implicit none
 integer :: i, system, RENAME
 character(len=10) :: cas_prog = ' '
 character(len=24) :: data_string = ' '
 character(len=240) :: inpname, pyname

 if(.not. excited) return
 if(casscf) then
  data_string = 'SA-CASSCF'
  cas_prog = casscf_prog
 else
  data_string = 'DMRG-SA-CASSCF'
  cas_prog = dmrgscf_prog
 end if
 write(6,'(A)',advance='no') TRIM(data_string)
 write(6,'(A,2(I0,A))') '(',nacte,'e,',nacto,'o) using program '//&
                           TRIM(cas_prog)

 i = index(hf_fch, '.fch', back=.true.)
 select case(TRIM(cas_prog))
 case('pyscf')
  inpname = hf_fch(1:i-1)//'.py'
  pyname = hf_fch(1:i-1)//'_SA.py'
  i = system('bas_fch2py '//TRIM(hf_fch))
  i = RENAME(TRIM(inpname), TRIM(pyname))
  call prt_sacas_script_into_py(pyname, hf_fch)
 case('gaussian')
  inpname = hf_fch(1:i-1)//'_SA.gjf'
  call prt_sacas_gjf(inpname, hf_fch)
 case('orca')
  inpname = hf_fch(1:i-1)//'_SA.inp'
  call prt_sacas_orca_inp(inpname, hf_fch)
 case('gamess')
  inpname = hf_fch(1:i-1)//'_SA.inp'
  call prt_sacas_gms_inp(inpname, hf_fch)
 case('dalton')
  inpname = hf_fch(1:i-1)//'_SA.dal'
  call prt_sacas_dalton_inp(inpname, hf_fch)
 case default
  write(6,'(A)') 'Other CASSCF programs are not supported currently.'
  stop
 end select
end subroutine do_sa_cas

! print (DMRG-)SA-CASSCF script into a given .py file
subroutine prt_sacas_script_into_py(pyname, gvb_fch)
 use mol, only: nacto, nacta, nactb
 use mr_keyword, only: mem, nproc, casscf, dmrgscf, maxM, hardwfn, crazywfn, &
  dkh2_or_x2c, RI, RIJK_bas, hf_fch, nstate
 implicit none
 integer :: i, fid1, fid2, RENAME
 character(len=21) :: RIJK_bas1
 character(len=240) :: buf, pyname1, cmofch
 character(len=240), intent(in) :: pyname, gvb_fch

 if(RI) call auxbas_convert(RIJK_bas, RIJK_bas1, 1)
 pyname1 = TRIM(pyname)//'.tmp'
 i = index(pyname, '.py', back=.true.)

 open(newunit=fid1,file=TRIM(pyname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(pyname1),status='replace')

 do while(.true.)
  read(fid1,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do

 if(dmrgscf) then
  write(fid2,'(A)') 'from pyscf import mcscf, dmrgscf, lib'
 else
  write(fid2,'(A)') 'from pyscf import mcscf, lib'
 end if
 write(fid2,'(A)') 'from py2fch import py2fch'
 write(fid2,'(A)') 'from shutil import copyfile'
 write(fid2,'(A,/)') 'import numpy as np'
 if(dmrgscf) then
  write(fid2,'(A,I0,A)') "dmrgscf.settings.MPIPREFIX ='mpirun -n ",nproc,"'"
 end if
 write(fid2,'(A,I0,A1,/)') 'lib.num_threads(',nproc,')'

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(1:9) == 'mf.kernel') exit
  write(fid2,'(A)') TRIM(buf)
 end do
 write(fid2,'(A,I0,A)') 'mf.max_memory = ', mem*1000, ' # MB'
 write(fid2,'(A)') 'mf.kernel()'

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do
 close(fid1,status='delete')

 ! mem*500 is in fact mem*1000/2. The mc.max_memory and fcisolver.max_memory seem
 ! not to share common memory, they used two memory, so I have to make them half
 if(casscf) then ! SA-CASSCF
  write(fid2,'(3(A,I0),A)',advance='no') 'mc = mcscf.CASSCF(mf,', nacto,',(',nacta,',',nactb,')'
  if(dkh2_or_x2c) then
   write(fid2,'(A)') ').x2c1e()'
  else
   if(RI) then
    write(fid2,'(A)') ").density_fit(auxbasis='"//TRIM(RIJK_bas1)//"')"
   else
    write(fid2,'(A)') ')'
   end if
  end if
  write(fid2,'(A,I0,A)') 'mc.fcisolver.max_memory = ',mem*200,' # MB'
 else ! DMRG-SA-CASSCF
  write(fid2,'(A)',advance='no') 'mc = dmrgscf.DMRGSCF(mf,'
  if(dkh2_or_x2c) then
   write(fid2,'(3(I0,A))') nacto,',(',nacta,',',nactb,')).x2c1e()'
  else
   write(fid2,'(3(I0,A))') nacto,',(',nacta,',',nactb,'))'
  end if
  write(fid2,'(A,I0)') 'mc.fcisolver.maxM = ', maxM
  write(fid2,'(A,I0,A)') 'mc.fcisolver.memory = ',CEILING(DBLE(mem)/DBLE((5*nproc))),' # GB'
 end if

 write(fid2,'(A)',advance='no') 'mc = mc.state_average_(['
 write(fid2,'(5(A,I0,A))',advance='no') ('1e0/',nstate,'e0,',i=1,nstate,1)
 write(fid2,'(A)') '0.0,0.0,0.0])'
 write(fid2,'(A,I0,A)') 'mc.max_memory = ', mem*800, ' # MB'
 write(fid2,'(A)') 'mc.max_cycle = 200'
 write(fid2,'(A,I0)') 'mc.fcisolver.spin = ', nacta-nactb
 write(fid2,'(A,I0,A)') 'mc.fix_spin_(ss=',nacta-nactb,')'

 if(hardwfn) then
  write(fid2,'(A)') 'mc.fcisolver.max_cycle = 200'
 else if(crazywfn) then
  write(fid2,'(A)') 'mc.fcisolver.level_shift = 0.2'
  write(fid2,'(A)') 'mc.fcisolver.pspace_size = 1200'
  write(fid2,'(A)') 'mc.fcisolver.max_space = 100'
  write(fid2,'(A)') 'mc.fcisolver.max_cycle = 300'
 else
  write(fid2,'(A)') 'mc.fcisolver.max_cycle = 100'
 end if
 write(fid2,'(A)') 'mc.verbose = 5'
 write(fid2,'(A)') 'mc.kernel()'

 i = index(hf_fch, '.fch', back=.true.)
 cmofch = hf_fch(1:i-1)//'_SA.fch'
 write(fid2,'(/,A)') '# save MOs into .fch file'
 write(fid2,'(A)') "copyfile('"//TRIM(gvb_fch)//"', '"//TRIM(cmofch)//"')"
 write(fid2,'(A)') 'noon = np.zeros(nif)'
 write(fid2,'(A)') "py2fch('"//TRIM(cmofch)//"',nbf,nif,mc.mo_coeff,'a',noon,False)"
 ! mc.mo_occ only exists for PySCF >= 1.7.4
 close(fid2)

 i = RENAME(TRIM(pyname1), TRIM(pyname))
 return
end subroutine prt_sacas_script_into_py

! print SA-CASSCF input file of Gaussian
subroutine prt_sacas_gjf(gjfname, hf_fch)
 use mol, only: nacto, nacte
 use mr_keyword, only: mem, nproc, nstate
 use util_wrapper, only: unfchk
 implicit none
 integer :: i, fid
 real(kind=8), allocatable :: weight(:)
 character(len=240) :: chkname
 character(len=240), intent(in) :: gjfname, hf_fch

 i = index(gjfname, '.gjf', back=.true.)
 chkname = gjfname(1:i-1)//'.chk'
 call unfchk(hf_fch, chkname)

 open(newunit=fid,file=TRIM(gjfname),status='replace')
 write(fid,'(A)') '%chk='//TRIM(chkname)
 write(fid,'(A,I0,A)') '%mem=',mem,'GB'
 write(fid,'(A,I0)') '%nprocshared=',nproc
 write(fid,'(2(A,I0),A)') '#p CASSCF(',nacte,',',nacto,',StateAverage) chkbasis&
  & nosymm int=nobasistransform'
 write(fid,'(A,/)') 'scf(maxcycle=500) guess=read geom=allcheck'
 allocate(weight(nstate))
 weight = 1d0/DBLE(nstate)
 do i = 1, nstate, 1
  write(fid,'(F10.8)') weight(i)
 end do ! for i
 deallocate(weight)

 write(fid,'(A)',advance='no')
 close(fid)
end subroutine prt_sacas_gjf

! print SA-CASSCF input file of ORCA
subroutine prt_sacas_orca_inp(inpname, hf_fch)
 use util_wrapper, only: fch2gbw
 use mol, only: nacto, nacte, mult
 use mr_keyword, only: mem, nproc, dkh2_or_x2c, RI, RIJK_bas, nstate
 implicit none
 integer :: i, fid, fid1, system, RENAME
 character(len=240) :: buf, gbwname, inpname1
 character(len=240), intent(in) :: inpname, hf_fch

 i = index(inpname, '.inp', back=.true.)
 inpname1 = inpname(1:i-1)//'.in'
 gbwname = inpname(1:i-1)//'.gbw'
 call fch2gbw(hf_fch, gbwname)

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')
 write(fid1,'(A,I0,A)') '%pal nprocs ', nproc, ' end'
 write(fid1,'(A,I0)') '%maxcore ', CEILING(1000d0*DBLE(mem)/DBLE(nproc))
 write(fid1,'(A)',advance='no') '! CASSCF'
 ! RIJK in CASSCF must be combined with CONVentional
 if(RI) write(fid1,'(A)',advance='no') ' RIJK conv '//TRIM(RIJK_bas)
 write(fid1,'(A)') ' TightSCF'
 if(dkh2_or_x2c) then
  write(fid1,'(A)') '%rel'
  write(fid1,'(A)') ' method DKH'
  write(fid1,'(A)') ' order 2'
  write(fid1,'(A)') 'end'
 end if
 write(fid1,'(A)') '%casscf'
 write(fid1,'(A,I0)') ' nel ', nacte
 write(fid1,'(A,I0)') ' norb ', nacto
 write(fid1,'(A,I0)') ' mult ', mult
 write(fid1,'(A,I0)') ' nroots ', nstate
 write(fid1,'(A)') ' maxiter 200'
 if(RI) write(fid1,'(A)') ' TrafoStep RI'
 write(fid1,'(A)') 'end'

 do i = 1, 3   ! skip 3 lines
  read(fid,'(A)') buf
 end do ! for i

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine prt_sacas_orca_inp

! print SA-CASSCF input file of GAMESS
subroutine prt_sacas_gms_inp(inpname, hf_fch)
 implicit none
 character(len=240), intent(in) :: inpname, hf_fch
end subroutine prt_sacas_gms_inp

! print SA-CASSCF input file of Dalton
subroutine prt_sacas_dalton_inp(inpname, hf_fch)
 implicit none
 character(len=240), intent(in) :: inpname, hf_fch
end subroutine prt_sacas_dalton_inp

