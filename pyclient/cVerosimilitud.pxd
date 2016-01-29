from libcpp.vector cimport vector

cdef extern from "Verosimilitud.h":
	cdef cppclass Verosimilitud:
		Verosimilitud()
		double LogLikelihood(vector[double] pp,vector[double] np)
		double OscillationProbability(double energy,double zenith, vector[double] pp)
		vector[double] SendSheep()
		void PrintSheep(vector[double])
