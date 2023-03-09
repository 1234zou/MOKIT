#!/usr/bin/env python
# Copyright 2020-2023 The MOKIT Developers. All Rights Reserved.

import os, sys, subprocess
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
from distutils.command.build_scripts import build_scripts
from distutils import sysconfig

CLASSIFIERS = [
'Development Status :: 5 - Production/Stable',
'Intended Audience :: Science/Research',
'Intended Audience :: Developers',
'License :: OSI Approved :: Apache Software License',
'Programming Language :: Fortran',
'Programming Language :: Python :: 3',
'Topic :: Software Development',
'Topic :: Scientific/Engineering',
'Operating System :: POSIX',
]

NAME             = 'mokit'
MAINTAINER       = 'Jingxiang Zou'
MAINTAINER_EMAIL = 'njumath@sina.cn'
DESCRIPTION      = 'MOKIT: Molecular Orbital KIT'
LONG_DESCRIPTION = '''
MOKIT offers various utilities and modules to transfer MOs among various quantum chemistry software packages. 

Besides, the automr program in MOKIT can set up and run common multi-reference calculations in a black-box way.
'''
URL              = 'http://gitlab.com/jxzou/mokit'
DOWNLOAD_URL     = 'http://gitlab.com/jxzou/mokit'
LICENSE          = 'Apache License 2.0'
AUTHOR           = 'MOKIT developers'
AUTHOR_EMAIL     = 'njumath@sina.cn'
PLATFORMS        = ['Linux']
def get_version():
    topdir = os.path.abspath(os.path.join(__file__, '..'))
    with open(os.path.join(topdir, 'mokit', '__init__.py'), 'r') as f:
        for line in f.readlines():
            if line.startswith('__version__'):
                delim = '"' if '"' in line else "'"
                return line.split(delim)[1]
    raise ValueError("Version string not found")
VERSION = get_version()

#class MakeExt(Extension):
#    def __init__(self, name, cmdir="."):
#        Extension.__init__(self, name, [])
#        #self.cmake_lists_dir = os.path.abspath(cmdir)

from distutils.dep_util import newer
#from distutils import log
from stat import ST_MODE
from distutils.util import convert_path
class BinBuild(build_scripts):
    def initialize_options(self):
        build_scripts.initialize_options(self)
        self.build_temp = None

    def finalize_options(self):
        build_scripts.finalize_options(self)
        self.set_undefined_options("build", ("build_temp", "build_temp"))

    def copy_scripts(self):
        #build_scripts.copy_scripts(self)
        outfiles, updated_files = [], []
        subprocess.check_call(['mkdir', '-p', self.build_dir])
        for script in self.scripts:
            script = convert_path(script)
            outfile = os.path.join(self.build_dir, os.path.basename(script))
            #print(self.build_dir, outfile)
            outfiles.append(outfile)

            if not self.force and not newer(script, outfile):
                print("not copying %s (up-to-date) "% script)
                return
    
            updated_files.append(outfile)
            #print("copying and adjusting %s -> %s"%( script, self.build_dir))
            self.copy_file(script, outfile)
        if os.name == "posix":
            for file in outfiles:
                if self.dry_run:
                    print("changing mode of %s"% file)
                else:
                    oldmode = os.stat(file)[ST_MODE] & 0o7777
                    newmode = (oldmode | 0o555) & 0o7777
                    if newmode != oldmode:
                        print(
                            "changing mode of %s from %o to %o"%( file, oldmode, newmode)
                        )
                        os.chmod(file, newmode)


class MakeBuildExt(build_ext):
    def run(self):
        extension = self.extensions[0]
        assert extension.name == 'mokit_lib_placeholder'
        self.build_extensions()

    def build_extensions(self):
        self.announce('Building binaries', level=3)
        # Do not use high level parallel compilation
        print("Python3: ", sys.executable)
        print("Build Dir: ", os.getcwd())
        #extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmd = ['make', 'all', '-f', 'Makefile.gnu_mkl_conda']
        #self.spawn(cmd)
        subprocess.check_call(cmd, cwd='./src')

    # To remove the infix string like cpython-37m-x86_64-linux-gnu.so
    # Python ABI updates since 3.5
    # https://www.python.org/dev/peps/pep-3149/
#    def get_ext_filename(self, ext_name):
#        ext_path = ext_name.split('.')
#        filename = build_ext.get_ext_filename(self, ext_name)
#        name, ext_suffix = os.path.splitext(filename)
#        return os.path.join(*ext_path) + ext_suffix
    def get_ext_filename(self, ext_name):
        ext_path = ext_name.split('.')
        filename = build_ext.get_ext_filename(self, ext_name)
        name, ext_suffix = os.path.splitext(filename)
        return os.path.join(*ext_path) + ext_suffix

from distutils.command.build import build

build.sub_commands = [c for c in build.sub_commands if c[0] == "build_ext"] + [
    c for c in build.sub_commands if c[0] != "build_ext"
]

def find_exes():
    topdir = os.path.abspath(os.path.join(__file__, '..'))
    with open(os.path.join(topdir, 'src', 'Makefile.main'), 'r') as f:
        read = False
        all_exes = []
        while(True): 
            line = f.readline()
            if line.startswith('exe:'):
                read = True
            if read:
                if len(line.strip())==0:
                    break
                exes = line.split()
                for exe in exes:
                    if exe != 'exe:' and exe != '\\': all_exes.append('bin/'+exe)
    print('all exes', all_exes)
    return all_exes + ['bin/mo_svd', 'bin/get_mokit_loc.py']


setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=LONG_DESCRIPTION,
    url=URL,
    download_url=DOWNLOAD_URL,
    project_urls={
        "Documentation": "http://jeanwsr.gitlab.io/mokit-doc-mdbook",
        "Source Code": "http://gitlab.com/jxzou/mokit/src",
        "Bug Tracker": "http://gitlab.com/jxzou/mokit/-/issues",
    },
    license=LICENSE,
    classifiers=CLASSIFIERS,
    keywords='quantum-chemistry computational-chemistry',
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    platforms=PLATFORMS,
    include_package_data=True,
    packages=find_packages(),
    ext_modules = [Extension('mokit_lib_placeholder', [])],
    cmdclass={'build_ext': MakeBuildExt, 'build_scripts': BinBuild},
    install_requires=[
        'numpy>=1.18',
        'mkl>=2020',
        'libgfortran5'],
    #extras_require={'h5py':['h5py>=2.7']}
    scripts = find_exes()
)
