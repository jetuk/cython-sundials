# -*- coding: utf-8 -*-

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from cython_generators import createKinsolProperties
import os
import numpy as np

print np.get_include()

sundials_lib = '.'#~/.local/lib/'
sundials_inc = '.'
#sundials_lib = r'D:\apps\sundials-2.5.0\lib'
#sundials_inc = r'D:\apps\sundials-2.5.0\include'

#sundials_static_libs = []

#for f in os.listdir(sundials_lib):
#    if 'libsundials' in f and os.path.splitext(f)[1] == '.a':
#        sundials_static_libs.append( os.path.join(sundials_lib,f) )
        
#print sundials_static_libs

#os.utime("sundials.pyx", None)
#os.utime("cvode.pyx", None)
#os.utime("kinsol.pyx", None

createKinsolProperties(os.path.join('pySundials','kinsol_properties.pxi'))

debug = False
if debug:
    compile_args= []
else:
    compile_args = ["-O2",]
    
    
extra_libraries = ['lapack','blas']

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
                    libraries=['sundials_nvecserial',]+extra_libraries,
                    include_dirs=[np.get_include(),sundials_inc],
                    library_dirs=[sundials_lib,],
                    extra_compile_args=compile_args,
                    
            ),
    Extension("pySundials.cvode", ["pySundials/cvode.pyx",],
                    libraries=['sundials_cvode',]+extra_libraries,
                    include_dirs=[np.get_include(),sundials_inc],
                    library_dirs=[sundials_lib,], 
                    extra_compile_args=compile_args,
                    
            ),
    Extension("pySundials.kinsol", ["pySundials/kinsol.pyx",],
                    libraries=['sundials_kinsol',]+extra_libraries,
                    include_dirs=[np.get_include(),sundials_inc],
                    library_dirs=[sundials_lib,], 
                    extra_compile_args=compile_args,
                    
            ),
            ],
)

                                                 #'sundials_cvode',
                                                 #'sundials_kinsol',