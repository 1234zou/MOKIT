#! /bin/bash
# Written by Jingxiang Zou at 20171222.
# This Shell script can modify the Gaussian source code (only G09D.01 and G16A.03 tested)
#  to do GVB up to 100 pairs.

echo 'Warning! This Shell script will modify your current Gaussian source code.'
echo 'You may want to save a copy of your source code(if it is important).'
read -p 'Continue? [y/n]:' select

if [ $select = "n" ]; then
 exit
fi

echo 'Modifying the source code...'
sed -i 's/IOrb(110)/IOrb(220)/g'     l506.F
sed -i 's/IDIAG(110)/IDIAG(220)/ig'  l506.F
sed -i 's/\<DIAG(110)/DIAG(220)/g'   l506.F
sed -i 's/NConf(110)/NConf(220)/g'   l506.F
sed -i 's/MSymPr(110)/MSymPr(220)/g' l506.F
sed -i 's/Delta(110)/Delta(220)/ig'  l506.F
sed -i 's/\<F(101)/F(201)/g'         l506.F
sed -i 's/AJ(101,101)/AJ(201,201)/g' l506.F
sed -i 's/AK(101,101)/AK(201,201)/g' l506.F
sed -i 's/AJ(10201)/AJ(40401)/g'     l506.F
sed -i 's/AK(10201)/AK(40401)/g'     l506.F
sed -i 's/NHOrb(2,101)/NHOrb(2,201)/g'   l506.F
sed -i 's/IHOrb(2,101)/IHOrb(2,201)/g'   l506.F
sed -i 's/IShell(2,101)/IShell(2,201)/g' l506.F
sed -i 's/IFreez(101)/IFreez(201)/g'     l506.F
sed -i 's/MaxOrb\/110\//MaxOrb\/220\//g' l506.F
sed -i 's/MaxHam\/101\//MaxHam\/201\//g' l506.F

sed -i 's/IDel(110)/IDel(220)/g'   l506.F
sed -i 's/AClear(110,Delta)/AClear(220,Delta)/g' l506.F
sed -i "s/'WILL USE',I5,1X,'MATRIX ELEMENTS OUT OF',I5/'WILL USE',I7,1X,'MATRIX ELEMENTS OUT OF',I7/" l506.F

sed -i 's/Coef(100)/Coef(200)/g'   l506.F l9999.F
sed -i 's/IPr(100)/IPr(200)/g'     l506.F l9999.F
sed -i 's/IPOrb(100)/IPOrb(200)/g' l506.F l9999.F
sed -i 's/NRt(100)/NRt(200)/g'     l506.F l9999.F
sed -i 's/INX(100)/INX(200)/g'     l506.F l9999.F
sed -i 's/MxPair\/100\//MxPair\/200\//g' l506.F l9999.F
echo 'Modification finished.'

#echo 'Recompile all the source code files...'
#echo 'This will take about 25min...'
#csh
#bsd/bldg09 |tee make.log
#exit
#echo 'Finish recompling. Please check whether new Gaussian has been compiled successfully.'

