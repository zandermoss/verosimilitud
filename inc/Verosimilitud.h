#ifndef __VEROSIMILITUD_H_INCLUDED__
#define __VEROSIMILITUD_H_INCLUDED__

#include "Tools.h"
#include <vector>
#include <math.h>
#include <Python.h>

class Verosimilitud {

	public:

		Verosimilitud(unsigned int my_numneu);
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
   
// 	   Set_data()
    
//	   Set_simulation()
};

#endif // __VEROSIMILILTUD_H_INCLUDED__
