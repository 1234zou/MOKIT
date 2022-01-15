*[English](README.md)*  
*[下载程序](https://gitlab.com/jxzou/mokit/-/archive/master/mokit-master.zip)*
# Molecular Orbital KIT (MOKIT)
MOKIT提供各种小程序和模块，用以实现在常见量子化学软件间传递分子轨道。除此之外。
MOKIT中的automr程序可以进行多参考（态）方法的自动化、黑箱式计算。

MOKIT中重要的小程序及其功能请见下图  
![MOKIT utilities with their functions](doc/orbital_transfer_CN.png)

利用MOKIT中的automr程序（结合上图中的各个小程序），您可以十分简单地（即像DFT计算
一样）进行多参考计算，并可充分利用每个量子化学软件最强的功能。例如，如下组合

  UHF(UNO) -> CASSCF -> CASPT2  
  Gaussian&emsp;&emsp;PySCF&emsp;&emsp;OpenMolcas  
or  
  UHF(UNO) -> GVB   -> CASSCF -> NEVPT2  
  Gaussian&emsp;&emsp;GAMESS&emsp;&emsp;PySCF&emsp;&emsp;PySCF  
or   
  RHF      -> GVB   -> CASSCF -> ic-MRCISD+Q  
  Gaussian&emsp;&emsp;GAMESS&emsp;&emsp;PySCF&emsp;&emsp;OpenMolcas

整个过程都是自动的。MOKIT在不同量化程序间传轨道时，考虑了基函数角动量的顺序问题（最
高支持H角动量，相当于C原子用cc-pV5Z基组，Zn原子用cc-pVQZ基组），因此同一种理论方法
（例如CASSCF）在不同量化程序中的电子能量可以很好地复现（误差通常小于10^-6 a.u.），且
几乎1-2圈收敛。

开发者还提供`Windows* 系统`下预编译好的约20个小程序，点击[下载](https://gitlab.com/jxzou/mokit/-/releases).
但请注意这些小程序的版本会滞后于master主分支代码，因此仍建议下载源码到Linux下自行编译。

请注意，尽管MOKIT程序的目标是使多参考计算实现自动化和黑箱式，无需人为干预。但用户
仍需具备使用常见量子化学软件的基本技能（例如熟悉Gaussian软件的常规DFT计算）。若
您是一名量化新手，强烈建议先学习并熟练使用Gaussian软件做常规计算，否则很可能难以
正确理解MOKIT的输出内容，或做出错误解读。

2022年1月15号

安装
----------

* 前提（编译器和库要求）
    - Fortran编译器: ifort(推荐>=2015) 或 gfortran(>=4.8.5 且 <8.0)
    - Intel MKL (推荐安装Intel编译器，内含ifort和MKL库)
    - f2py (推荐安装Anaconda Python3，内含f2py)

* 编译全部模块， 执行

        cd src
        make all

* 若仅想使用单个子程序，可仅编译单个子程序或模块  
  例如执行

        make fch2inp

* 在执行'make all'之后, 你需要设置三个环境变量`MOKIT_ROOT`, `PATH` 和 `PYTHONPATH`.  
  例如，假定你的MOKIT安装在/home/$USER/software/mokit目录, 你需要在~/.bashrc文件中设定
  以下变量（ORCA和GAMESS可执行文件的路径请按照您机器上的实际情况修改）:

        export MOKIT_ROOT=/home/$USER/software/mokit
        export PATH=$MOKIT_ROOT/bin:$PATH
        export PYTHONPATH=$MOKIT_ROOT/lib:$PYTHONPATH
        export ORCA=/opt/orca-5.0.1/orca
        export GMS=/home/$USER/software/gamess/rungms

* 原始GAMESS程序只能处理少于13对的GVB计算，但借助MOKIT现今可以实现上百对的GVB计算。
  因此请阅读[手册](doc/)4.4.10部分使用提供的脚本自动修改GAMESS代码。

快速开始
----------
* 每个小程序的使用十分简单，直接运行即可在屏幕上打印出使用说明。例如在Shell中运行小程
  序`fch2inp`，输出如下

   ERROR in subroutine fch2inp: wrong command line arguments!  
   Example 1 (R(O)HF, UHF, CAS): fch2inp a.fch  
   Example 2 (GVB)             : fch2inp a.fch -gvb [npair]  
   Example 3 (ROGVB)           : fch2inp a.fch -gvb [npair] -open [nopen]

* 对于lib/目录下Python动态库文件的使用方法，请阅读examples/utilities/目录下的readme.txt

* 自动做多参考计算的核心程序automr的输入文件采用的是Gaussian gjf文件的格式。例如，一个O-H
  键长为1.5 A的水分子输入文件'00-h2o_cc-pVDZ_1.5.gjf'示例如下

        %mem=4GB
        %nprocshared=4
        #p CASSCF/cc-pVDZ
        
        mokit{}
        
        0 1
        O      -0.23497692    0.90193619   -0.068688
        H       1.26502308    0.90193619   -0.068688
        H      -0.73568721    2.31589843   -0.068688

  只需在Shell中执行

        automr 00-h2o_cc-pVDZ_1.5.gjf >& 00-h2o_cc-pVDZ_1.5.out

  命令，automr程序会相继执行HF，GVB和CASSCF等计算,自动确定活性空间为CAS(4,4)。更多例子
  请见[examples](examples/)。

AutoMR支持调用的量子化学程序
----------
* [Gaussian](http://gaussian.com/)
* [PySCF](https://github.com/pyscf/pyscf)
* [GAMESS](https://www.msg.chem.iastate.edu/gamess/index.html)
* [OpenMolcas](https://gitlab.com/Molcas/OpenMolcas)
* [Molpro](https://www.molpro.net/)
* [ORCA](https://orcaforum.kofo.mpg.de)
* [BDF](http://182.92.69.169:7226/Introduction)
* [PSI4](https://github.com/psi4/psi4/)
* [Dalton](https://gitlab.com/dalton/dalton)

温馨提示
----------
* 若您想提供.fch(k)文件给automr程序读入，请务必在计算前在Gaussian的输入文件中加上
  关键词nosymm int=nobasistransform，以避免后续产生不必要的、不可预见的错误。

汇报Bug
----------
* 若您发现MOKIT的程序错误或bug，或有任何使用建议，可通过电子邮件njumath[at]sina.cn
  联系作者jxzou。在邮件中请将您的相关文件（例如.gjf, .fchk, .out文件等）添加为附件
  发送。由于作者正在博士延期中，不保证及时回复。

* 您也可以在此页面[Issues](https://gitlab.com/jxzou/mokit/-/issues)上新建一个问题。

* 还可加入MOKIT用户交流QQ群，群号：470745084

下一步计划（可能）
----------
* 支持BAGEL, CFOUR, NWCHEM等软件间传轨道

* 支持更多的多参考方法（如ic-MRCC）

* 开发和实现多参考的激发态计算

如何引用
----------
* 若您在您的研究中使用了MOKIT的任何一个子程序，请务必按如下引用

   Jingxiang Zou, Molecular Orbital Kit (MOKIT), https://gitlab.com/jxzou/mokit (accessed month day, year)

* 若您使用MOKIT进行了GVB方法的相关计算，请务必引用以下2篇文献

   DOI: 10.1021/acs.jctc.8b00854; DOI: 10.1021/acs.jpca.0c05216.

* 更详细的引用说明请参见文件[doc/](doc/)。您的规范引用是对开发者的极大鼓励。

Disclaimer
----------
Copyright (c) 2022 jxzou

All rights reserved.

Redistribution and use in source and binary forms are permitted provided that the above copyright notice and this paragraph are duplicated in all such forms and that any documentation, advertising materials, and other materials related to such distribution and use acknowledge that the software was developed by the author. The name of the authors may not be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

