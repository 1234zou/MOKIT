# Molecular Orbital KIT (MOKIT)
MOKIT offers various utilities and modules to transfer MOs among various quantum
chemistry software packages. Besides, the automr program in MOKIT can set up and
run common multi-reference calculations in a block-box way.

With MOKIT, one can perform multi-reference calculations in a quite simple way,
and utilize the best modules of each program. E.g.

  UHF(UNO) -> CASSCF -> CASPT2  
  Gaussian&emsp;&emsp;PySCF&emsp;&emsp;OpenMolcas  
or  
  UHF(UNO) -> GVB   -> CASSCF -> NEVPT2  
  Gaussian&emsp;&emsp;GAMESS&emsp;&emsp;PySCF&emsp;&emsp;PySCF

Negligible energy loss(usually<1e-6 a.u., for the same wave function method in two
programs) are ensured during transferring MOs, since the basis order of angular
momentum up to H(i.e. l=5) are considered.

Note that although MOKIT aims to make the multi-reference calculations block-box,
the users are still required to have practical experiences of quantum chemistry
computations (e.g. familiar with routine calculations in Gaussian). You are encouraged
to learn how to use Gaussian first if you are a fresh hand.

2020-06-23

Installation
------------

* Prerequisites
    - ifort (>=2015, recommended) or gfortran (>=4.8.5)
    - Intel MKL (installing Intel compilers recommended)
    - f2py (installing Anaconda Python3 recommended)

* Compile core module

        cd src
        make all

* Compile individual utility or module
  E.g.
        make fch2inp

* After 'make all', you need to set environment variables `MOKIT_ROOT`, `PATH` and `PYTHONPATH`.
  E.g. if MOKIT is installed in /home/$USER/mokit, the following should be set in ~/.bashrc:

        export MOKIT_ROOT=/home/$USER/mokit
        export PATH=$MOKIT_ROOT/bin:$PATH
        export PYTHONPATH=$MOKIT_ROOT/lib:$PYTHONPATH

Quick Start
-----------
* Each utility is self-explanatory. For example, run `fch2inp` in Shell,
  you will find

   ERROR in subroutine fch2inp: wrong command line arguments!  
   Example 1 (R(O)HF, CAS): fch2inp a.fch  
   Example 2 (UHF)        : fch2inp a.fch -uhf  
   Example 3 (GVB)        : fch2inp a.fch -gvb [npair]  
   Example 4 (ROGVB)      : fch2inp a.fch -gvb [npair] -open [nopen]

* For usages of modules in lib/, see examples/

* The input syntax of the automr program is like Gaussian gjf. For example, the input
  file 'N2_cc-pVQZ_6D10F_4.0.gjf' of the N2 molecule at d(N-N)=4.0A is shown below
  (assuming a stable UHF wave function is hold in the Gaussian .fch(k) file 'N2_cc-pVQZ_6D10F_4.0_uhf.fchk'):

        %mem=16GB
        %nprocshared=16
        #p NEVPT2/cc-pVQZ

        mokit{readuhf='N2_cc-pVQZ_6D10F_4.0_uhf.fchk',ist=1} 

  There is no need to write Cartesian coordinates since they are already provided in .fchk file. Run
  'automr N2_cc-pVQZ_6D10F_4.0.gjf >& N2_cc-pVQZ_6D10F_4.0.out' in SHELL. The automr
  program will successively perform GVB, CASCI, CASSCF, NEVPT2 computations. See examples/
  for more examples.

Some Tips
---------
* To avoid unnecessary errors, it is strongly recommended to specify keyword
  'nosymm int=nobasistransform' in Gaussian.

* Since only spherical harmonic functions are supported in ORCA, one should specify
  '5D 7F' keyword in Gaussian if he/she wants to call ORCA.

* Since only Cartesian functions are supported in GAMESS, one should specify '6D 10F'
  keyword in Gaussian if he/she wants to call GAMESS.

* Both spherical harmonic/Cartesian functions are supported in OpenMolcas, thus
  one can use either of them by specifying '5D 7F'/'6D 10F' in Gaussian (do not
  use the default 6-31G(d,p) keyword in Gaussian since it corresponds to '6D 7F')

Bug Report
----------
* If you find any bug frequently occurs, please contact the author jxzou via njumath[at]sina.cn,
  with your .fch(k) file and output files attached. The author may not answer or update code
  frequently since he is being postponed due to his PhD program.

TODO
----
* MOs trasferring among BAGEL, PSI4, NWCHEM

* Support more multireference methods like ic-MRCC, MS-CASPT2

* Develope/Implement black-box strategies of excited state calculations

Citation
--------
* If you use MOKIT in your work, please cite the gitlab page of MOKIT.

* If you use MOKIT to peform calculations involving GVB, please also cite DOI: 10.1021/acs.jctc.8b00854.

Disclaimer
----------
Copyright (c) 2020 jxzou

All rights reserved.

Redistribution and use in source and binary forms are permitted provided that the above copyright notice and this paragraph are duplicated in all such forms and that any documentation, advertising materials, and other materials related to such distribution and use acknowledge that the software was developed by the author. The name of the authors may not be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

