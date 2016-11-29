from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "Verosimilitud.h":

	cdef cppclass Verosimilitud:
		#Verosimilitud(unsigned int numneu,unsigned int loyear, unsigned int hiyear, str flux_path, str effective_area_path, str detector_correction_path, osc_func, callback)
		Verosimilitud(unsigned int numneu,unsigned int loyear, unsigned int hiyear, char* data_path, char* flux_path, char* effective_area_path, char* detector_correction_path, vector[double] nu_vec, vector[double] antinu_vec)

		vector[double] MinLLH(vector[double] param, vector[double] low_bound, vector[double] high_bound, vector[bool] param_to_minimize)

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

