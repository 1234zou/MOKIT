#!/usr/bin/env python
# Copyright 2020-2023 The MOKIT Developers. All Rights Reserved.

import os
from setuptools import setup
from setuptools.command.build_ext import build_ext

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
'Operating System :: Unix',
'Operating System :: MacOS',
]

NAME             = 'mokit'
MAINTAINER       = 'Jingxiang Zou'
MAINTAINER_EMAIL = 'njumath@sina.cn'
DESCRIPTION      = 'MOKIT: Molecular Orbital KIT'
URL              = 'http://gitlab.com/jxzou/mokit'
DOWNLOAD_URL     = 'http://gitlab.com/jxzou/mokit'
LICENSE          = 'Apache License 2.0'
AUTHOR           = 'MOKIT developers'
AUTHOR_EMAIL     = 'njumath@sina.cn'
PLATFORMS        = ['Linux', 'Mac OS-X', 'Unix']
def get_version():
    topdir = os.path.abspath(os.path.join(__file__, '..'))
    with open(os.path.join(topdir, 'lib', '__init__.py'), 'r') as f:
        for line in f.readlines():
            if line.startswith('__version__'):
                delim = '"' if '"' in line else "'"
                return line.split(delim)[1]
    raise ValueError("Version string not found")
VERSION = get_version()

class MakeBuildExt(build_ext):
    def build_make(self, extension):
        self.announce('Building binaries', level=3)
        # Do not use high level parallel compilation
        cmd = ['make', 'all']
        self.spawn(cmd)

    # To remove the infix string like cpython-37m-x86_64-linux-gnu.so
    # Python ABI updates since 3.5
    # https://www.python.org/dev/peps/pep-3149/
#    def get_ext_filename(self, ext_name):
#        ext_path = ext_name.split('.')
#        filename = build_ext.get_ext_filename(self, ext_name)
#        name, ext_suffix = os.path.splitext(filename)
#        return os.path.join(*ext_path) + ext_suffix

setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=DESCRIPTION,
    url=URL,
    download_url=DOWNLOAD_URL,
    project_urls={
        "Documentation": "http://gitlab.com/jxzou/mokit/doc",
        "Source Code": "http://gitlab.com/jxzou/mokit/src",
        "Bug Tracker": "http://gitlab.com/jxzou/mokit/-/issues",
    },
    license=LICENSE,
    classifiers=CLASSIFIERS,
    keywords='quantum-chemistry computational-chemistry',
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    platforms=PLATFORMS,
    packages=['src'],
    cmdclass={'build_ext': MakeBuildExt},
    install_requires=['numpy>=1.13,!=1.16,!=1.17'],
    extras_require={'h5py':['h5py>=2.7']}
)
