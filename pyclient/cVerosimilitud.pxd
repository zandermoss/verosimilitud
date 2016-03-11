from libcpp.vector cimport vector

cdef extern from "Verosimilitud.h":

	ctypedef double (*pyoscfunc)(vector[double] argument, void* user_data)

	cdef cppclass Verosimilitud:
		Verosimilitud(unsigned int numneu)
		double LogLikelihood(vector[double] pp,vector[double] np)
		void SetDecayStructure(vector[vector[double]] dcy_lambda)
		vector[double] Likelihood()
		vector[double] GetEproxEdges()
		vector[double] GetCosZenithEdges()
		vector[unsigned int] GetExpDims()
		void SetMassStructure(vector[vector[double]] pmns_lambda)
		void SetDecayEigenvalues(vector[double] dcy_eig)
		void PrintThing()
		double OscillationProbability(double energy,double zenith)
		void SetDeSolver(pyoscfunc de_solver, void* user_data)
		
