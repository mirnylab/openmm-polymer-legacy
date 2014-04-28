from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize


ext = Extension(
    "polymerCython",                 # name of extension
    ["polymerCython.pyx","qcprot.c"],           # filename of our Pyrex/Cython source
    language="c",              # this causes Pyrex/Cython to create C++ source    

    include_dirs=["."],
    library_dirs=["."],
    cmdclass = {'build_ext': build_ext}
    )


setup(
    name = "polymerCython",
    cmdclass={"build_ext":build_ext},
    ext_modules=[ext]
)
