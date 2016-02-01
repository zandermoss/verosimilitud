#include "Verosimilitud.h"
#include <iostream>
#include <stdexcept>
#include <Python.h>



		Verosimilitud::Verosimilitud(unsigned int my_numneu)
		{
			numneu=my_numneu;	
			std::cout << "BLEAT" << std::endl;
		}
		/*
		std::vector<double> Verosimilitud::SendSheep(void)
		{
			std::vector<double> sheep;
			double lamb = 0;
			for (unsigned int i=0; i<3; i++)
			{
				lamb+=0.1;
				sheep.push_back(lamb);
			}
			return sheep;
		}
		*/

        void Verosimilitud::SetDecayStructure(std::vector<std::vector<double> > my_dcy_lambda)
		{
	        if((my_dcy_lambda.size() != my_dcy_lambda[0].size()) 
				|| (my_dcy_lambda.size() != numneu))
			{
	        	throw std::runtime_error("Decay lambda matrix has the wrong dimensions.");
			}

			for (unsigned int i=0; i<my_dcy_lambda.size(); i++)
			{
				decay_lambda.push_back(my_dcy_lambda[i]);
			}	
		}

        void Verosimilitud::SetMassStructure(std::vector<std::vector<double> > my_pmns_lambda)
		{
	        if((my_pmns_lambda.size() != my_pmns_lambda[0].size()) 
				|| (my_pmns_lambda.size() != numneu))
			{
	        	throw std::runtime_error("Decay lambda matrix has the wrong dimensions.");
			}

			for (unsigned int i=0; i<my_pmns_lambda.size(); i++)
			{
				pmns_lambda.push_back(my_pmns_lambda[i]);
			}	
		}

        void Verosimilitud::SetDecayEigenvalues(std::vector<double> my_dcy_eig)
		{
	        if(my_dcy_eig.size() != numneu)
			{
	        	throw std::runtime_error("Decay eigenvalue vector has the wrong dimension.");
			}

			decay_eig=std::vector<double>(my_dcy_eig);
		}

		void Verosimilitud::Print1D(std::vector<double> array)
		{
			PrintArray(array);
		}
		void Verosimilitud::Print2D(std::vector<std::vector<double> > matrix)
		{
			PrintMatrix(matrix);
		}
		void Verosimilitud::PrintThing(void)
		{
			PrintMatrix(decay_lambda);
		}
/*
	    double LogLikelihood(std::vector<double> pp,std::vector <double> np)
		{
			double mu;
			double loglike=0;
			for (unsigned int ei=0; ei<esize; ei++)
			{
				for (unsigned int zi=0; zi<zsize; zi++)
				{
					mu = Measurement_Interval*Flux(reco_energy[ei],reco_zenith[zi])*
						X_sec(reco_energy[ei]) * OscProb(reco_energy[ei],reco_zenith[zi]); 
					
					loglike+=LogPoissonProbability(data[ei][zi],mu);
				}
			}
	
			//loglike+=LogNuisancePriors(np); Implement this bayesian behavior in 
			//a child class: verosimilitud_bayesian
	
	        mu[recoenergy][recozenith] = mu(pp,np) = sum_true_energy_true_zenith simulation[recoenergy][recozenith][trueenergy][truezenith];
	    // sum over each of the bins in marray::data
	    // here this loop/sum is over recoenergy and recozenith bins
	        return -sum(LogPoissonProbability(xdata,mu[recoenergy][recozenith])) + sum(LogNuisancePriors);
	    }
*/

	    double Verosimilitud::OscillationProbability(double energy,double zenith)
		{
	   	 	//Here we call nusheep. Although, would it be better to precalculate an array?
			//I think for now I will implement direct calls.
			std::vector<double> argument;
			argument.push_back(energy);
			argument.push_back(zenith);
			double osc_prob = de_solver(argument,user_data);	
			std::cout << "  ZENITH: " << zenith << "  ENERGY: " << energy << "  PROB: " << osc_prob << std::endl;				
			return osc_prob;
   		}

		void Verosimilitud::SetDeSolver(pyoscfunc my_de_solver, void* my_user_data)
		{
			de_solver=my_de_solver;
			user_data=my_user_data;		
		}



    
   
		/* 
		std::vector<double> LogNuisancePriors(np)
		{
			//

        	return std::vector of nuisance parameter penalization
        	return {LogGaussianProbability(np[0],np[0]_prior), ... }
		}
		*/	


		/*
		double Flux(double energy, double zenith)
		{
			//Return from data set?
			//How do flux uncertanties enter into this calculation?
			//For now, restrict to muon flavor only
		}
		*/


		/*
		double X_sec(double energy)
		{
			//Again, from data-set?
			//Restrict to muon only
		} 
		*/ 
  
 //	   Set_data()
    
//	   Set_simulation()
