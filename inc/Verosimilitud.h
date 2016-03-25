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

class Verosimilitud {

	public:

		Verosimilitud(unsigned int my_numneu);
		~Verosimilitud();
	    double LogLikelihood(std::vector<double> pp,std::vector<double> np);
		typedef double (*pyoscfunc)(std::vector<double> argument, void *userdata);
	    double OscillationProbability(double energy,double zenith);
		void SetDecayStructure(std::vector<std::vector<double> > my_dcy_lambda);
		void SetDeSolver(pyoscfunc de_solve, void* user_data);
		void SetMassStructure(std::vector<std::vector<double> > my_pmns_lambda);
		void SetDecayEigenvalues(std::vector<double> my_dcy_eig);
        void Print1D(std::vector<double> array);
        void Print2D(std::vector<std::vector<double> > matrix);
        void PrintThing(void);
		void testhdf(void);
		std::vector<unsigned int> GetExpDims(void);
		std::vector<double> GetEproxEdges(void);
		std::vector<double> GetCosZenithEdges(void);
		std::vector<double> Likelihood(void);
		std::vector<double> CalculateExpectation(void);

		

	protected:
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
		std::vector<double>* eprox_edges;
		std::vector<double>* coszenith_edges;

		Tensor* expectation;



		const unsigned int NeutrinoEnergyBins=280;
    	const unsigned int CosZenithBins=11;
    	const unsigned int EnergyProxyBins=50;

	    EffectiveArea* eff_area;
	    Tensor* area;
	    Tensor* NeutrinoEnergyEdges;
	    Tensor* CosZenithEdges;
	    Tensor* EnergyProxyEdges;

		Tensor* data;
		
		ICData* icd;
	
	    ConventionalFlux* conv_flux;
	    Tensor* flux;
	    Tensor* detcorr;
	






   
// 	   Set_data()
    
//	   Set_simulation()
};

#endif // __VEROSIMILILTUD_H_INCLUDED__
