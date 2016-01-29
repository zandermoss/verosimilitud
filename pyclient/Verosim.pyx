cimport cVerosimilitud
from libcpp.vector cimport vector
import numpy as np


cdef class Verosim:
	cdef cVerosimilitud.Verosimilitud* _c_verosimilitud
	def __cinit__(self):
		self._c_verosimilitud = <cVerosimilitud.Verosimilitud *>new cVerosimilitud.Verosimilitud()
		#self._c_verosimilitud = new cVerosimilitud.Verosimilitud()
		if self._c_verosimilitud is NULL:
			raise MemoryError()

	def getsheep(self):
		cdef vector[double] v
		v= self._c_verosimilitud.SendSheep()
		return np.asarray(v)
	
	def printsheep(self,vector[double] sheep):
		self._c_verosimilitud.PrintSheep(sheep)
		return

	def __dealloc__(self):
		cdef cVerosimilitud.Verosimilitud *temp
		if self._c_verosimilitud is not NULL:
			temp=<cVerosimilitud.Verosimilitud *> self._c_verosimilitud
			del temp
