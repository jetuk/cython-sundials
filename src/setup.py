# -*- coding: utf-8 -*-

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os
import numpy as np

print np.get_include()

#sundials_lib = '/home/james/.local/lib/'
#sundials_inc = '/home/james/.local/include/'

#sundials_static_libs = []

#for f in os.listdir(sundials_lib):
#    if 'libsundials' in f and os.path.splitext(f)[1] == '.a':
#        sundials_static_libs.append( os.path.join(sundials_lib,f) )
        
#print sundials_static_libs

#os.utime("sundials.pyx", None)
#os.utime("cvode.pyx", None)
#os.utime("kinsol.pyx", None

debug = False
if debug:
    compile_args= []
else:
    compile_args = ["-O3",]

setup(
    name='pySundials',
    version='0.1',
    description='Python bindings for Sundials (CVODE, KINSOL)',
    author='James E Tomlinson',
    author_email='tomo.bbe@gmail.com',
    packages=['pySundials',],
    cmdclass = {'build_ext': build_ext},
    package_data = {'pySundials' : ['sundials.pxd', 'cvode.pxd', 'kinsol.pxd',
                                    'libsundials.pxd', 'libcvode.pxd', 'libkinsol.pxd']},
    ext_modules = [
    
    Extension("pySundials.sundials", ["pySundials/sundials.pyx",],
                    libraries=['sundials_nvecserial','lapack','blas'],
                    include_dirs=[np.get_include(),],
                    #library_dirs=[sundials_lib,],
                    extra_compile_args=compile_args,
                    pyrex_gdb=debug,
            ),
    Extension("pySundials.cvode", ["pySundials/cvode.pyx",],
                    libraries=['sundials_cvode','lapack','blas'],
                    include_dirs=[np.get_include(),],
                    #library_dirs=[sundials_lib,], 
                    extra_compile_args=compile_args,
                    pyrex_gdb=debug,
            ),
    Extension("pySundials.kinsol", ["pySundials/kinsol.pyx",],
                    libraries=['sundials_kinsol','lapack','blas'],
                    include_dirs=[np.get_include(),],
                    #library_dirs=[sundials_lib,], 
                    extra_compile_args=compile_args,
                    pyrex_gdb=debug,
            ),
            ],
)

                                                 #'sundials_cvode',
                                                 #'sundials_kinsol',