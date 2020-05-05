# Molecular Orbital KIT (MOKIT)
MOKIT offers various utilities and modules to transfer MOs between various quantum
chemistry software packages. Currently supported: Gaussian, GAMESS, ORCA, PySCF
and OpenMolcas.

With MOKIT, one can easily perform multi-reference calculations in a block-box
way, and utilize the best modules of each program. E.g.

  UHF(UNO) -> CASSCF -> CASPT2
  Gaussian    PySCF     OpenMolcas
or
  UHF(UNO) -> GVB   -> CASSCF -> NEVPT2
  Gaussian    GAMESS   PySCF     PySCF

Negligible energy loss(usually<1e-6 a.u., for the same wave function method in two
programs) are ensured during transferring MOs, since the basis order of angular
momentum up to H(i.e. l=5) are considered.

Note that although MOKIT aims to make the multi-reference calculations block-box,
the users are still required to have practical experiences of quantum chemistry
computations (e.g. familiar with routine calculations in Gaussian). You are encouraged
to learn how to use Gaussian first if you are a fresh hand.

2020-05-05

Installation
------------

* Prerequisites
    - ifort (>=2015, recommended) or gfortran (>=4.8.5)
    - Intel MKL (installing Intel compilers recommended)
    - f2py (installing Anaconda Python3 recommended)

* Compile core module

        cd src
        make all
        make install

* Compile individual utility or module
  E.g.
        make fch2inp

  In this case you have to manually move fch2inp to ../bin.

* After 'make install', you need to set environment variable `PATH` and `PYTHONPATH`.
  E.g. if MOKIT is installed in /home/$USER/mokit, the following should be set in ~/.bashrc:

        export PATH=/home/$USER/mokit/bin:$PATH
        export PYTHONPATH=/home/$USER/mokit/lib:$PYTHONPATH

Quick Start
-----------
* Each utility is self-explanatory. For example, run `fch2inp` in Shell,
  you will find

   ERROR in subroutine fch2inp: wrong command line arguments!
   Example 1 (R(O)HF, CAS): fch2inp a.fch
   Example 2 (UHF)        : fch2inp a.fch -uhf
   Example 3 (GVB)        : fch2inp a.fch -gvb [npair]
   Example 4 (ROGVB)      : fch2inp a.fch -gvb [npair] -open [nopen0]

* For usages of modules in lib/, see examples/

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
* If you find any bug frequently occurs, please contact the author jxzou at njumath[at]sina.cn,
  with your .fch(k) file and output files attached. The author may not answer or update code
  frequently since he is being postponed due to his PhD program.

* Better not eat apples everyday, since the old saying goes "An apple a day, keeps the doctor away".

TODO
----
* MOs trasferring among BAGEL, PSI4, NWCHEM

* A whole framework of doing multi-reference calculations in a black-box way

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

