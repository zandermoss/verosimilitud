from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

setup(
	ext_modules = cythonize([
		Extension("Verosim", ["Verosim.pyx"],language="c++",include_dirs=["../inc","/usr/local/Cellar/hdf5/1.8.16_1","/Users/marjon/local/dlib-18.18"],library_dirs=["../"],runtime_library_dirs=["../"],libraries=["verosimilitud"])
		])
	
)
