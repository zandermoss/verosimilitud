from libcpp.vector cimport vector

cdef extern from "Verosimilitud.h":

	ctypedef double (*pyoscfunc)(vector[double] argument, void* user_data)

	cdef cppclass Verosimilitud:
		#Verosimilitud(unsigned int numneu,unsigned int loyear, unsigned int hiyear, str flux_path, str effective_area_path, str detector_correction_path, osc_func, callback)
		Verosimilitud(unsigned int numneu,unsigned int loyear, unsigned int hiyear, str flux_path, str effective_area_path, str detector_correction_path)

		void SetEproxCuts(vector[double] cuts)
		void CalculateExpectation()
		vector[double] Chi2MinNuisance(vector[double] nuisance)
		#vector[double] GetDataVec(double scale)
		vector[double] GetDataVec()
		vector[double] GetPertExpectationVec(vector[double] nuisance)
#		vector[double] GetExpectationVec()
		vector[double] GetAreaVec()
		vector[double] GetFluxVec()
		#vector[double] GetExpectationVec()
		double SimpsAvg(double coszmin, double coszmax, double emin, double emax, double anti, int nintervals)
		void SetSimpsNIntervals(int nintervals)

		vector[double] GetEproxEdges()
		vector[double] GetEnergyEdges()
		vector[double] GetCosZenithEdges()
		vector[unsigned int] GetExpDims()
		vector[unsigned int] GetDataDims()

		double OscillationProbability(double energy,double zenith, double anti)
		void SetDeSolver(pyoscfunc de_solver, void* user_data)
		
