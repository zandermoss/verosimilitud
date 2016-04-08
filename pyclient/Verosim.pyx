cimport cVerosimilitud
from libcpp.vector cimport vector
import numpy as np



cdef double callback(vector[double] argument, void *f):
	return (<object>f)(<object>argument)

cdef class Verosim:
	cdef cVerosimilitud.Verosimilitud* _c_verosimilitud
	def __cinit__(self,unsigned int numneu, year):
		if year=="2010":
			loyear=0
			hiyear=1
		elif year=="2011":
			loyear=1
			hiyear=2
		elif year=="both":
			loyear=0
			hiyear=2
		else:
			raise ValueError("Bad datayears argument: choices are '2010', '2011', or 'both'")

		self._c_verosimilitud = <cVerosimilitud.Verosimilitud *>new cVerosimilitud.Verosimilitud(numneu,loyear,hiyear)
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
	def GetEnergyEdges(self):
		return self._c_verosimilitud.GetEnergyEdges()

	def GetCosZenithEdges(self):
		return self._c_verosimilitud.GetCosZenithEdges()


	def CalculateExpectation(self):
		self._c_verosimilitud.CalculateExpectation()

	def Chi2MinNuisance(self,vector[double] nuisance):
		return self._c_verosimilitud.Chi2MinNuisance(nuisance)

	#def GetDataVec(self,scale):
	#	return self._c_verosimilitud.GetDataVec(scale)
	def GetDataVec(self):
		return self._c_verosimilitud.GetDataVec()

	def GetPertExpectationVec(self,vector[double] nuisance):
		return self._c_verosimilitud.GetPertExpectationVec(nuisance)
	def GetExpectationVec(self):
		return self._c_verosimilitud.GetExpectationVec()
	def GetAreaVec(self):
		return self._c_verosimilitud.GetAreaVec()
	def GetFluxVec(self):
		return self._c_verosimilitud.GetFluxVec()

	def SetEproxCuts(self,vector[double] cuts):
		self._c_verosimilitud.SetEproxCuts(cuts)

	def SetSimpsNIntervals(self,nintervals):
		self._c_verosimilitud.SetSimpsNIntervals(nintervals)

	def SimpsAvg(self,coszmin, coszmax, emin, emax, anti, nintervals):
		return self._c_verosimilitud.SimpsAvg(coszmin, coszmax, emin, emax, anti, nintervals)


	"""
	def SetDeSolver(self,object obj):
		self._c_verosimilitud.SetDeSolver(obj)
		return
	"""	

	def SetDeSolver(self,f):
		self._c_verosimilitud.SetDeSolver(callback,<void*>f)
		return

	def OscillationProbability(self,energy,zenith, anti):
		self._c_verosimilitud.OscillationProbability(energy,zenith,anti)
		return

		
	def __dealloc__(self):
		cdef cVerosimilitud.Verosimilitud *temp
		if self._c_verosimilitud is not NULL:
			temp=<cVerosimilitud.Verosimilitud *> self._c_verosimilitud
			del temp
