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
		std::vector<unsigned int> GetDataDims(void);
		std::vector<double> GetEproxEdges(void);
		std::vector<double> GetCosZenithEdges(void);
		std::vector<double> CalculateExpectation(void);


		std::vector<double> GetExpectationVec(void);
		std::vector<double> GetExpectationVec(std::vector<double> nuisance);
		std::vector<double> GetDataVec(void);
		std::vector<double> GetDataVec(double scale);

		

		double Chi2MinNuisance(std::vector<double>* param);
		double Chi2(const dlib::matrix<double,0,1>& nuisance);
		dlib::matrix<double,0,1> Chi2Gradient(const dlib::matrix<double,0,1>& nuisance);


		void SetEproxCuts(std::vector<double> cuts);
		void SetCoszCuts(std::vector<double> cuts);

		const unsigned int NeutrinoEnergyBins=280;
    	const unsigned int CosZenithBins=11;
    	const unsigned int EnergyProxyBins=50;

		

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
		std::vector<unsigned int> dat_dimvec;
		std::vector<double>* eprox_edges;
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
		
		ICData* icd;
	
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





#endif // __VEROSIMILILTUD_H_INCLUDED__
