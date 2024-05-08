! written by jxzou at 20231109: partition a large molecule into monomers
!  and generate initial MOs as well as amplitudes from fragement calculations

program frag_qc
 implicit none
 integer :: wfn_type
 integer, parameter :: n = 2
 integer :: wfn_type0(n)
 character(len=240) :: fchname0(n), fchname
 character(len=24) :: data_string
 logical :: pos(n)

 wfn_type0 = 1; wfn_type=1; pos = .true.
 fchname0 = ['Cu_def2TZVP.fch','Cu_2.54_def2TZVP.fch']
 fchname = 'Cu_dimer_def2TZVP.fch'

 call fdate(data_string)
 write(6,'(A)') TRIM(data_string)

! call direct_sum_frag_fock_in_fch(n, fchname0, wfn_type0, pos, fchname, wfn_type)
! call direct_sum_frag_mo_in_fch(n, fchname0, wfn_type0, pos, fchname, wfn_type)
 call sum_frag_density_and_prt_into_fch(n, fchname0, pos, fchname)

 call fdate(data_string)
 write(6,'(/,A)') TRIM(data_string)
end program frag_qc

