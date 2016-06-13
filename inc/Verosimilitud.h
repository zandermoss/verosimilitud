#ifndef __VEROSIMILITUD_H_INCLUDED__
#define __VEROSIMILITUD_H_INCLUDED__

#include "Tools.h"
#include <vector>
#include <math.h>
#include <Python.h>
#include "Tensor.h"
#include "EffectiveArea.h"
#include "ICData.h"
#include "ConventionalFlux.h"
#include <functional>

#include <dlib/optimization.h>

//! The log-likelihood/Chi-squared calculation class.
/*!
	The Verosimilitud class has three essential functions:
	-1
		Calculate a histogram of expected rate over reconstructed energy and angle. This calculation uses survival probabilities from NuSHEEP as well as effective area and flux information from the data release. The specifics of the given oscillation model are stipulated in the calling script using the interface to NuSHEEP.
	-2
		Load IceCube data from the IceCube release.
	-3
		Calculate a Chi^2 corresponding to the prediction of the specified oscillation model and the IC (IceCube) data. Saturated and unsaturated log-likelihood values are calculated, and then used to approximate the corresponding chi^2. This chi^2 is then minimized over the relevant nuisance parameters with penalties distributed according to prior nuisance distributions.

FIXME add links to the relevant functions!
*/

class Verosimilitud {

	public:
		//--------------------------------------------------------//
		//! The Constructor.
		/*!
			The constructor for Verosimilitud initializes data structures and loads in IC data as well as effective detector area and atmospheric muon neutrino fluxes. Data is loaded and stored in an ICData object. Effective area is loaded and stored in an EffectiveArea object, and flux is loaded and stored in a ConventionalFlux object. Default energy and zenith angle cuts are applied (no cut). Finally, bin centers are calculated from "::"NeutrinoEnergyEdges", "::"CosZenithEdges, and "::"EnergyProxyEdges for use in simulation.
			\param my_numneu An integer specifying the number of neutrinos, or, more precisely, the dimension of the neutrino state space.
			\param loyear An integer providing a lower bound for the loops over data years.
			\param hiyear An integer providing an upper bound for the loops over data years. If [loyear,hiyear]=[0,2], both the 2010 and 2011 datasets are included. If [0,1], only 2010, and if [1,2], only 2011.
		*/


		Verosimilitud(unsigned int my_numneu,unsigned int loyear, unsigned int hiyear);


		

		//--------------------------------------------------------//
		//! The Destructor
		/*!
			Just cleaning up dynamic memory. No special operations.
		*/

		~Verosimilitud();


	
		//--------------------------------------------------------//
		//! Typedef for python function for NuSHEEP survival probability calculation.
		/*!
			Some arcane C trickery for passing relevant function data using the void pointer userdata. I never actually came to understand how this works. I just swiped the code from stackexchange:
			\param userdata void pointer used for sordid data transmission.
			\param argument	the actual argument, whose nature will be clarified in the documentation surrounding "::"OscillationProbability
		*/

		typedef double (*pyoscfunc)(std::vector<double> argument, void *userdata);

		//! Function wrapping a call to the NuSHEEP survival probability calculation.
			

	    double OscillationProbability(double energy,double zenith, double anti);

		void SetDecayStructure(std::vector<std::vector<double> > my_dcy_lambda);
		void SetDeSolver(pyoscfunc de_solve, void* user_data);

		std::vector<unsigned int> GetExpDims(void);
		std::vector<unsigned int> GetDataDims(void);
		std::vector<double> GetEproxEdges(void);
		std::vector<double> GetEnergyEdges(void);
		std::vector<double> GetCosZenithEdges(void);
		void CalculateExpectation(void);


		void SetSimpsNIntervals(int nintervals);


		std::vector<double> GetFluxVec(void);
		std::vector<double> GetAreaVec(void);

		std::vector<double> GetExpectationVec(void);
		std::vector<double> GetPertExpectationVec(std::vector<double> nuisance);
		std::vector<double> GetDataVec(void);
		//std::vector<double> GetDataVec(double scale);

		

		std::vector<double> Chi2MinNuisance(std::vector<double> nuisance);
		double Chi2(const dlib::matrix<double,0,1>& nuisance);
		dlib::matrix<double,0,1> Chi2Gradient(const dlib::matrix<double,0,1>& nuisance);

		double TestFunction(double e, double z, double a);


		double SimpsAvg(double coszmin, double coszmax, double emin, double emax, double anti, int nintervals);


		void SetEproxCuts(std::vector<double> cuts);
		void SetCoszCuts(std::vector<double> cuts);




		const unsigned int NeutrinoEnergyBins=280;
    	const unsigned int CosZenithBins=11;
    	const unsigned int EnergyProxyBins=50;

		

	protected:


		unsigned int data_years[2];

		int simps_nintervals;

    	std::vector<double>  pp;
    	std::vector<double>  np;
   		
		unsigned int numneu;
		std::vector<std::vector<double> > decay_lambda;
		std::vector<std::vector<double> > pmns_lambda;
		std::vector<double> decay_eig;

		pyoscfunc de_solver;
		void* user_data;	
 
//    	marray data: reco energy, reco zenith // histogram
 //   	marray simulation: reco energy, true energy, reco zenith, true zenith // histogram

		//initialize bins for energy and zenith angles. Here I will use lower bounds.
		std::vector<double> ebins;
		std::vector<double> zbins;
		//unsigned int esize=ebins.size(); 
		//unsigned int zsize=zbins.size(); 

		std::vector<double> LogNuisancePriors(std::vector<double> np);
		double Flux(double energy, double zenith);
		double X_sec(double energy);


		std::vector<unsigned int> exp_dimvec;
		std::vector<unsigned int> dat_dimvec;
		std::vector<double>* eprox_edges;
		std::vector<double>* energy_edges;
		std::vector<double>* coszenith_edges;

		std::vector<double>* eprox_centers;
		std::vector<double>* coszenith_centers;


		Tensor* expectation;




	    EffectiveArea* eff_area;
	    Tensor* area;
	    Tensor* NeutrinoEnergyEdges;
	    Tensor* CosZenithEdges;
	    Tensor* EnergyProxyEdges;

		Tensor* data;
		
		ICData* icd[2];
	
	    ConventionalFlux* conv_flux;
	    Tensor* flux;
	    Tensor* detcorr;
	



		std::vector<double>* eprox_cuts;
		std::vector<double>* cosz_cuts;
	
//Nuisance parameters:

	double norm_mean=1;
	double norm_sigma=0.4;
	double gamma_mean=0;
	double gamma_sigma=0.03;



   
// 	   Set_data()
    
//	   Set_simulation()
};

class Chi2_caller
{
	private:
		Verosimilitud* my_verosim;
	public:
		Chi2_caller(Verosimilitud * verosim)
		{
			my_verosim = verosim;
		}


		double operator() (const dlib::matrix<double,0,1>& nuisance) const
		{
			return my_verosim->Chi2(nuisance);
		}

};




class Chi2grad_caller
{
	private:
		Verosimilitud* my_verosim;
	public:

		Chi2grad_caller(Verosimilitud * verosim)
		{
			my_verosim = verosim;
		}


		dlib::matrix<double,0,1> operator() (const dlib::matrix<double,0,1>& nuisance) const
		{
			return my_verosim->Chi2Gradient(nuisance);
		}

};



#endif // __VEROSIMILILTUD_H_INCLUDED__
