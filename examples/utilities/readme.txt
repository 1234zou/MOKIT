---------------------------------------------------
Example 1: transfer RHF MOs from Gaussian to GAMESS
---------------------------------------------------
Run 'fch2inp 01-water_cc-pVTZ_6D10F.fch' in Shell, a file 01-water_cc-pVTZ_6D10F.inp would be generated.

(Usually the keywords in .inp file are already appropriate, but do check if advanced options are desired)

---------------------------------------------------
Example 2: transfer UHF MOs from Gaussian to GAMESS
---------------------------------------------------
Run 'fch2inp 02-bis_6D10F.fchk' in Shell, a file 02-bis_6D10F.inp would be generated.

----------------------------------------------------------------
Example 2.1: transfer UHF MOs from Gaussian to MOLCAS/OpenMolcas
----------------------------------------------------------------
Run 'fch2inporb 02-bis_6D10F.fchk' in Shell, two files(.input and .INPORB) would be generated.

-------------------------------------------------------
Example 3: transfer paired LMOs from Gaussian to GAMESS
-------------------------------------------------------
Run 'fch2inp 03-ben_cc-pVDZ_6D10F_uno_assoc_rot.fch -gvb 3' in Shell, the .inp would be generated.
You can run
'gms 03-ben_cc-pVDZ_6D10F_uno_assoc_rot.inp 01 8 >& 03-ben_cc-pVDZ_6D10F_uno_assoc_rot.gms'
to submit this job. This is a GVB(3) job of GAMESS with paired LMOs as initial guess.

-----------------------------------------------------
Example 3.1: transfer GVB NOs from GAMESS to Gaussian
-----------------------------------------------------
After Example 3 finished, you can run the following 2 steps
(1) cp 03-ben_cc-pVDZ_6D10F_uno_assoc_rot.fch 03-ben_cc-pVDZ_6D10F_gvb3.fch
(2) dat2fch 03-ben_cc-pVDZ_6D10F_uno_assoc_rot.dat 03-ben_cc-pVDZ_6D10F_gvb3.fch -gvb 3
The GVB(3) orbitals are in 03-ben_cc-pVDZ_6D10F_gvb3.fch file, which can be visualizied by GaussView.

-------------------------------------------------
Example 4: transfer RHF MOs from Gaussian to ORCA
-------------------------------------------------
Run 'fch2mkl 04-PhCOOH_def2TZVP_rhf.fchk' in Shell, two files(.mkl and .inp) would be generated.
Then run 'orca_2mkl 04-PhCOOH_def2TZVP_rhf -gbw' in Shell, you can obtain 04-PhCOOH_def2TZVP_rhf.gbw file.
Finally run 'orca 04-PhCOOH_def2TZVP_rhf.inp >& 04-PhCOOH_def2TZVP_rhf.out', you will find that the RHF energy
can be exactly reproduced using the utility fch2mkl.

--------------------------------------------------
Example 5: transfer GVB NOs from Gaussian to PySCF
--------------------------------------------------
In Example 3.1, you already got the file 03-ben_cc-pVDZ_6D10F_gvb3.fch, now run
'python 05-ben_cc-pVDZ_6D10F_gvb3_2CAS66.py >& 05-ben_cc-pVDZ_6D10F_gvb3_2CAS66.out'
The GVB NOs will be read and used as the initial guess of CAS(6,6). After CASSCF finished,
CASSCF NOs will be exported.

(Note: the basis sets in .py file are obtained by running 'bas_fch2py 03-ben_cc-pVDZ_6D10F_uno_assoc_rot.fch'
in Shell. Do not use the built-in basis sets of PySCF.)

--------------------------------------------------
Example 6: transfer UHF MOs from Gaussian to PySCF
--------------------------------------------------
Run 'python 06-bis_6D10F.py' in Shell. The UHF canonical MOs in 02-bis_6D10F.fchk will be read in.
One will find that the UHF energy immediately converges after little change.

--------------------------------------------------
Example 7: create fch from PySCF
--------------------------------------------------
Run 'python 07-OH_writefch.py'. It will create `07-OH.fch` containing basis information and MOs.
To validate the fch, run 'unfchk 07-OH.fch && g16 07-OH_restart.gjf' and the Gaussian calculation will converge in 1 cycle.
Of course, the main purpose of generating this fch is to visualize it with GaussView, Multiwfn, etc.

----------------------------------------------------------
Example 8: transfer complex GHF MOs from Gaussian to PySCF
----------------------------------------------------------
Run 'bas_fch2py 08-test1198.fch'. It will create `08-test1198.py` which contains
Cartesian coordinates, basis information and a few of necessary keywords.  
Uncomment the last three lines in `08-test1198.py`, and run 'python 08-test1198.py',
PySCF will converge in 2 cycles.  
If you further append `mf.stability()` into the end of the file `08-test1198.py`,
you will find this complex GHF solution has an internal instability.

----------------------------------------------------------
Example 9: transfer complex GHF MOs from PySCF to Gaussian
----------------------------------------------------------
Please read file `09-test1198.py`, which is nothing but several lines are added
into `08-test1198.py`. The density matrix of the instable wave function in `08-test1198.py`
is slightly perturbed in order to obtain a stable wave function. And then MOs of
the stable wave function are exported into a .fch file.  
To validate this .fch, run 'unfchk 09-test1198.fch && g16 09-test1198_restart.gjf'
and the Gaussian calculation will converge in ~19 cycle. But the electronic energy
changes only with 1.6e-9 a.u.

