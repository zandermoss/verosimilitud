from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

setup(
	ext_modules = cythonize([
		Extension("Verosim", ["Verosim.pyx"],language="c++",include_dirs=["../inc","/usr/local/hdf5/include","/home/pinkpig/src/dlib-18.18"],library_dirs=["../"],runtime_library_dirs=["../"],libraries=["verosimilitud"])
		])
	
)
