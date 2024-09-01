set -x -e
cd src
echo `which f2py`
make all -f Makefile.gnu_openblas_conda
cd ..
pip install -v --prefix=$PREFIX .
#mkdir conda_build
mv bin $PREFIX
