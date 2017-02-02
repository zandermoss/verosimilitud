cimport cVerosimilitud
from libcpp.vector cimport vector
import numpy as np


cdef class Verosim:
	cdef cVerosimilitud.Verosimilitud* _c_verosimilitud
	def __cinit__(self,unsigned int numneu, data_path, flux_path,effective_area_path, nu_vec, antinu_vec):

		self._c_verosimilitud = <cVerosimilitud.Verosimilitud *>new cVerosimilitud.Verosimilitud(numneu,data_path,flux_path,effective_area_path,nu_vec, antinu_vec)
		if self._c_verosimilitud is NULL:
			raise MemoryError()


	def MinLLH(self, param, low_bound, high_bound, param_to_minimize):
		return self._c_verosimilitud.MinLLH(param, low_bound, high_bound, param_to_minimize)

	def LLH(self, param):
		return self._c_verosimilitud.LLH(param)

	def LLHGrad(self, param, index):
		return self._c_verosimilitud.LLHGrad(param,index)

	def GetExpDims(self):
		return self._c_verosimilitud.GetExpDims()

	def GetDataDims(self):
		return self._c_verosimilitud.GetExpDims()

	def GetEproxEdges(self):
		return self._c_verosimilitud.GetEproxEdges()
	def GetEnergyEdges(self):
		return self._c_verosimilitud.GetEnergyEdges()

	def GetCosZenithEdges(self):
		return self._c_verosimilitud.GetCosZenithEdges()


	def CalculateExpectation(self):
		self._c_verosimilitud.CalculateExpectation()

#	def Chi2MinNuisance(self,vector[double] nuisance):
#		return self._c_verosimilitud.Chi2MinNuisance(nuisance)

	#def GetDataVec(self,scale):
	#	return self._c_verosimilitud.GetDataVec(scale)
	def GetDataVec(self):
		return self._c_verosimilitud.GetDataVec()

	def GetPertExpectationVec(self,vector[double] nuisance):
		return self._c_verosimilitud.GetPertExpectationVec(nuisance)
#	def GetExpectationVec(self):
#		return self._c_verosimilitud.GetExpectationVec()
	def GetAreaVec(self):
		return self._c_verosimilitud.GetAreaVec()
#	def GetFluxVec(self):
#		return self._c_verosimilitud.GetFluxVec()

	def SetEproxCuts(self,vector[double] cuts):
		self._c_verosimilitud.SetEproxCuts(cuts)

	def SetSimpsNIntervals(self,nintervals):
		self._c_verosimilitud.SetSimpsNIntervals(nintervals)

	def SimpsAvg(self,coszmin, coszmax, emin, emax, anti, nintervals):
		return self._c_verosimilitud.SimpsAvg(coszmin, coszmax, emin, emax, anti, nintervals)


		
	def __dealloc__(self):
		cdef cVerosimilitud.Verosimilitud *temp
		if self._c_verosimilitud is not NULL:
			temp=<cVerosimilitud.Verosimilitud *> self._c_verosimilitud
			del temp
