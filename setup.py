from setuptools import find_packages
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy 
cmdclass = {}
ext = [Extension(
    "openmmlib.polymerCython",                 # name of extension
    ["openmmlib/polymerCython.pyx","openmmlib/qcpSource/qcprot.c"],           # filename of our Pyrex/Cython source
    language="c",              # this causes Pyrex/Cython to create C++ source
    include_dirs=[".", numpy.get_include()],
    library_dirs=["."],
    )]


ext.append(Extension(
    "openmmlib.fastContacts",                 # name of extension
    ["openmmlib/fastContacts.pyx"],           # filename of our Pyrex/Cython source
    language="c++",              # this causes Pyrex/Cython to create C++ source    
    include_dirs=[".", numpy.get_include()],
    library_dirs=["."],
    ))


cmdclass.update({'build_ext': build_ext} )

setup(
    name='openmmlib',
    url='http://mirnylab.bitbucket.org/hiclib/index.html',
    description=('Hi-C data analysis library.'),
      ext_modules=ext,
      include_package_data=True,
       package_data = {
                   '': ['openmmlib/getCpu*']},
      cmdclass = cmdclass,
       packages=['openmmlib'],

)
