set -x -e
cd src
make all -f Makefile.gnu_mkl_conda
cd ..
pip install -v --prefix=$PREFIX .
#mkdir conda_build
mv bin $PREFIX
