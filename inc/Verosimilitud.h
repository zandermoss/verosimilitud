#ifndef __VEROSIMILITUD_H_INCLUDED__
#define __VEROSIMILITUD_H_INCLUDED__

#include "Tools.h"
#include <vector>
#include <math.h>

class Verosimilitud {

	public:

		Verosimilitud(void);
	    double LogLikelihood(std::vector<double> pp,std::vector<double> np);
	    double OscillationProbability(double energy,double zenith, std::vector<double> pp);
		std::vector<double> SendSheep(void);
		void PrintSheep(std::vector<double> sheep);

	protected:
    	std::vector<double>  pp;
    	std::vector<double>  np;
    
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
