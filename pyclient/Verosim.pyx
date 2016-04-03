cimport cVerosimilitud
from libcpp.vector cimport vector
import numpy as np



cdef double callback(vector[double] argument, void *f):
	return (<object>f)(<object>argument)

cdef class Verosim:
	cdef cVerosimilitud.Verosimilitud* _c_verosimilitud
	def __cinit__(self,unsigned int numneu):
		self._c_verosimilitud = <cVerosimilitud.Verosimilitud *>new cVerosimilitud.Verosimilitud(numneu)
		#self._c_verosimilitud = new cVerosimilitud.Verosimilitud()
		if self._c_verosimilitud is NULL:
			raise MemoryError()

	def SetDecayStructure(self, vector[vector[double]] dcy_lambda):
		self._c_verosimilitud.SetDecayStructure(dcy_lambda)
		return

	def SetMassStructure(self, vector[vector[double]] pmns_lambda):
		self._c_verosimilitud.SetMassStructure(pmns_lambda)
		return

	def SetDecayEigenvalues(self, vector[double] dcy_eig):
		self._c_verosimilitud.SetDecayEigenvalues(dcy_eig)
		return

	def PrintThing(self):
		self._c_verosimilitud.PrintThing()
		return

	def GetExpDims(self):
		return self._c_verosimilitud.GetExpDims()

	def GetDataDims(self):
		return self._c_verosimilitud.GetExpDims()

	def GetEproxEdges(self):
		return self._c_verosimilitud.GetEproxEdges()

	def GetCosZenithEdges(self):
		return self._c_verosimilitud.GetCosZenithEdges()


	def CalculateExpectation(self):
		self._c_verosimilitud.CalculateExpectation()

	def Chi2MinNuisance(self,vector[double] nuisance):
		return self._c_verosimilitud.Chi2MinNuisance(nuisance)

	def GetDataVec(self,scale):
		return self._c_verosimilitud.GetDataVec(scale)
	#def GetDataVec(self):
	#	return self._c_verosimilitud.GetDataVec()

	def GetExpectationVec(self,vector[double] nuisance):
		return self._c_verosimilitud.GetExpectationVec(nuisance)
	#def GetExpectationVec(self):
	#	return self._c_verosimilitud.GetExpectationVec()

	def SetEproxCuts(self,vector[double] cuts):
		self._c_verosimilitud.SetEproxCuts(cuts)


	"""
	def SetDeSolver(self,object obj):
		self._c_verosimilitud.SetDeSolver(obj)
		return
	"""	

	def SetDeSolver(self,f):
		self._c_verosimilitud.SetDeSolver(callback,<void*>f)
		return

	def OscillationProbability(self,energy,zenith):
		self._c_verosimilitud.OscillationProbability(energy,zenith)
		return

		
	def __dealloc__(self):
		cdef cVerosimilitud.Verosimilitud *temp
		if self._c_verosimilitud is not NULL:
			temp=<cVerosimilitud.Verosimilitud *> self._c_verosimilitud
			del temp
