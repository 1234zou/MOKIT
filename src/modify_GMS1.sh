#! /bin/sh
# Written by Jingxiang Zou at 20171106. This Shell script can modify the GAMESS source code
#  (only 2017 version tested) to do GVB up to 499 pairs.
# Note that this script must be run in the gamess/source directory; and your GAMESS is required to
#  be compiled before, so that objective files *.o are already in the gamess/object directory.
#  Only modified source code files will be recompiled.

echo 'Warning! This Shell script will modify your current GAMESS source code.'
echo 'You may want to save a copy of your source code(if it is important).'
read -p 'Continue? [y/n]:' select

if [ $select = "n" ]; then
 exit
fi

echo 'detecting if file modify_GMS2.f90 exists...'
if [ -e "modify_GMS2.f90" ]; then
 echo 'File modify_GMS2.f90 exists. Continue...'
else
 echo 'File modify_GMS2.f90 does not exist! This script must be used combined with modify_GMS2.f90'
 exit
fi

echo 'Modifying the source code...'

filelist='comp cphf cprohf fmoh2c fmohss grd1 grd2a guess gvb hess hss1c hss2a hss2b hss2c locpol mexing parley prppop qmfm scflib statpt vector vvos'

for fname in $filelist
do
 sed -i 's/\<CICOEF(2,12)/CICOEF(2,499)/g' $fname.src
 sed -i 's/\<F(25)/F(999)/g'               $fname.src
 sed -i 's/\<FGVB(25)/FGVB(999)/g'         $fname.src
 sed -i 's/\<ALPHA(325)/ALPHA(499500)/g'   $fname.src
 sed -i 's/\<BETA(325)/BETA(499500)/g'     $fname.src

 sed -i 's/\<TKIN(25)/TKIN(999)/g'         $fname.src
 sed -i 's/\<TPLUSV(25)/TPLUSV(999)/g'     $fname.src
 sed -i 's/\<CILOW(12)/CILOW(499)/g'       $fname.src
 sed -i 's/\<KCORB(2,12)/KCORB(2,499)/g'   $fname.src
 sed -i 's/\<CIHAM(91)/CIHAM(125250)/g'    $fname.src
done
sed -i 's/9148 FORMAT(1X,I2/9148 FORMAT(I3/g' gvb.src
sed -i "s/CICOEF(',I2,')=',F12.8,',',F12.8/CICOEF(',I3,')=',E17.10,',',E17.10/" gvb.src
sed -i 's/NHAMX\ =\ 25/NHAMX\ =\ 999/g'    scflib.src
sed -i 's/NPAIRX\ =\ 12/NPAIRX\ =\ 499/g'  scflib.src
sed -i 's/200) T/500) T/g' inputa.src

gfortran modify_GMS2.f90 -o modify_GMS2.exe
./modify_GMS2.exe
echo 'Modification finished.'

echo 'Recompile the modified source code...'

cd ../object/
for fname in $filelist
do
 rm -f $fname.o
done

cd ..
for fname in $filelist
do
 ./comp $fname
done
./comp inputa
echo 'Finish recompling.'

echo 'Try to link all object files...'
./lked gamess 01
echo 'The new GAMESS version will be gamess.01.x'

