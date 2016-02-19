#include "Verosimilitud.h"
#include <iostream>
#include <stdexcept>
#include <Python.h>
#include "Flux.h"
#include "Tensor.h"
#include "H5Cpp.h"
#include "EffectiveArea.h"

/*
//From effective_area_demo.c
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
//--------
*/


/*

//Read a dataset into a buffer, asusming that the allocated size is correct.
//If anything goes wrong just bail out.
void readDataSet(hid_t container_id, const char* path, double* targetBuffer){
	hid_t dataset_id = H5Dopen(container_id, path, H5P_DEFAULT);
	herr_t status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, targetBuffer);
	if(status<0){
		fprintf(stderr,"Failed to read dataset '%s'\n",path);
		exit(1);
	}
	H5Dclose(dataset_id);
}

void readDoubleAttr(hid_t container_id, const char* path, const char* name, double* targetBuffer){
	hid_t attr_id = H5Aopen_by_name(container_id, path, name, H5P_DEFAULT, H5P_DEFAULT);
	herr_t status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, targetBuffer);
	if(status<0){
		fprintf(stderr,"Failed to read attribute '%s::%s'\n",path,name);
		exit(1);
	}
	H5Aclose(attr_id);
}

int main(void){
	
	const unsigned int neutrinoEnergyBins=280;
	const unsigned int cosZenithBins=11;
	const unsigned int energyProxyBins=50;
	//the edges used for the bins of the effective area histograms
	//these are the same for all effective areas
	double* trueEnergyEdges=(double*)malloc((neutrinoEnergyBins+1)*sizeof(double));
	double* cosZenithEdges=(double*)malloc((cosZenithBins+1)*sizeof(double));
	double* energyProxyEdges=(double*)malloc((energyProxyBins+1)*sizeof(double));
	
	const unsigned int histogramDims[3]={neutrinoEnergyBins,cosZenithBins,energyProxyBins};
	
	multidim effArea2010NuMu=alloc_multi(3,histogramDims);
	multidim effArea2010NuMuBar=alloc_multi(3,histogramDims);
	multidim effArea2010NuTau=alloc_multi(3,histogramDims);
	multidim effArea2010NuTauBar=alloc_multi(3,histogramDims);
	multidim effArea2011NuMu=alloc_multi(3,histogramDims);
	multidim effArea2011NuMuBar=alloc_multi(3,histogramDims);
	multidim effArea2011NuTau=alloc_multi(3,histogramDims);
	multidim effArea2011NuTauBar=alloc_multi(3,histogramDims);
	
	multidim effArea2010NuMu_Err=alloc_multi(3,histogramDims);
	multidim effArea2010NuMuBar_Err=alloc_multi(3,histogramDims);
	multidim effArea2010NuTau_Err=alloc_multi(3,histogramDims);
	multidim effArea2010NuTauBar_Err=alloc_multi(3,histogramDims);
	multidim effArea2011NuMu_Err=alloc_multi(3,histogramDims);
	multidim effArea2011NuMuBar_Err=alloc_multi(3,histogramDims);
	multidim effArea2011NuTau_Err=alloc_multi(3,histogramDims);
	multidim effArea2011NuTauBar_Err=alloc_multi(3,histogramDims);
	
	double livetime2010;
	double livetime2011;
	
	//make some arrays of aliases to make it easy to loop over these things
	multidim* effectiveAreas[2][4]={{&effArea2010NuMu,&effArea2010NuMuBar,
	                                 &effArea2010NuTau,&effArea2010NuTauBar},
	                                {&effArea2011NuMu,&effArea2011NuMuBar,
	                                 &effArea2011NuTau,&effArea2011NuTauBar}};
	multidim* effectiveAreaErrs[2][4]={{&effArea2010NuMu_Err,&effArea2010NuMuBar_Err,
	                                    &effArea2010NuTau_Err,&effArea2010NuTauBar_Err},
	                                   {&effArea2011NuMu_Err,&effArea2011NuMuBar_Err,
	                                    &effArea2011NuTau_Err,&effArea2011NuTauBar_Err}};
	double* livetimes[2]={&livetime2010,&livetime2011};
	
	//the expected distribution of astrophysical events
	//in the energy proxy observable
	double* expectation=(double*)malloc(energyProxyBins*sizeof(double));
	double* expectationUncertainty=(double*)malloc(energyProxyBins*sizeof(double));
	
	herr_t status=0;
	hid_t file_id = H5Fopen("effective_area.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
	
	readDoubleAttr(file_id,"/2010","experimental_livetime(seconds)",&livetime2010);
	readDoubleAttr(file_id,"/2011","experimental_livetime(seconds)",&livetime2011);
	
	readDataSet(file_id, "/2010/nu_mu/bin_edges_0", trueEnergyEdges);
	readDataSet(file_id, "/2010/nu_mu/bin_edges_1", cosZenithEdges);
	readDataSet(file_id, "/2010/nu_mu/bin_edges_2", energyProxyEdges);
	
	readDataSet(file_id, "/2010/nu_mu/area", effArea2010NuMu.data);
	readDataSet(file_id, "/2010/nu_mu_bar/area", effArea2010NuMuBar.data);
	readDataSet(file_id, "/2010/nu_tau/area", effArea2010NuTau.data);
	readDataSet(file_id, "/2010/nu_tau_bar/area", effArea2010NuTauBar.data);
	readDataSet(file_id, "/2011/nu_mu/area", effArea2011NuMu.data);
	readDataSet(file_id, "/2011/nu_mu_bar/area", effArea2011NuMuBar.data);
	readDataSet(file_id, "/2011/nu_tau/area", effArea2011NuTau.data);
	readDataSet(file_id, "/2011/nu_tau_bar/area", effArea2011NuTauBar.data);
	
	readDataSet(file_id, "/2010/nu_mu/area_uncertainty", effArea2010NuMu_Err.data);
	readDataSet(file_id, "/2010/nu_mu_bar/area_uncertainty", effArea2010NuMuBar_Err.data);
	readDataSet(file_id, "/2010/nu_tau/area_uncertainty", effArea2010NuTau_Err.data);
	readDataSet(file_id, "/2010/nu_tau_bar/area_uncertainty", effArea2010NuTauBar_Err.data);
	readDataSet(file_id, "/2011/nu_mu/area_uncertainty", effArea2011NuMu_Err.data);
	readDataSet(file_id, "/2011/nu_mu_bar/area_uncertainty", effArea2011NuMuBar_Err.data);
	readDataSet(file_id, "/2011/nu_tau/area_uncertainty", effArea2011NuTau_Err.data);
	readDataSet(file_id, "/2011/nu_tau_bar/area_uncertainty", effArea2011NuTauBar_Err.data);
	
	H5Fclose(file_id);
	
	printf("2010 livetime: %lf s\n",livetime2010);
	printf("2011 livetime: %lf s\n",livetime2011);
	
	//a test flux to evaluate
	powerlawFlux flux;
	flux.normalization=1.63e-18;
	flux.index=-2.22;
	
	//To compute the expected number of observed events we need to:
	//1a. Multiply each effective area by the average flux in each bin
	//1b. Multiply by the phase space in each bin (true energy and solid angle)
	//1c. Multiply by the livetime for which the effective area is relevant
	//2. Sum over the effective areas for all particle types and detector configurations
	//3. Sum over dimensions not of interest (true neutrino energy, possibly zenith angle)
	//In this case we will compute the expectation as a function of the energy
	//proxy only, so we will project out both true energy and zenith angle.
	
	memset(expectation, 0, energyProxyBins*sizeof(double));
	memset(expectationUncertainty, 0, energyProxyBins*sizeof(double));
	unsigned int indices[3];
	
	for(unsigned int i=0; i<energyProxyBins; i++){
		indices[2]=i;
		//for both years
		for(unsigned int y=0; y<2; y++){
			//for each particle type
			for(unsigned int p=0; p<4; p++){
				//for each true energy bin
				for(unsigned int e=0; e<neutrinoEnergyBins; e++){
					indices[0]=e;
					
					double enMin=trueEnergyEdges[e];
					double enMax=trueEnergyEdges[e+1];
					
					//for each cosine zenith angle bin
					for(unsigned int z=0; z<cosZenithBins; z++){
						indices[1]=z;
						
						double cosZenithMin=cosZenithEdges[z];
						double cosZenithMax=cosZenithEdges[z+1];
						
						//the product of the average flux in the bin and the
						//phase space in the bin is simply the integral of the
						//flux over the bin
						double fluxIntegral=
						  integratePowerlawFlux(flux,enMin,enMax) //energy intergal
						  *(cosZenithMax-cosZenithMin) //zenith integral
						  *2*M_PI; //azimuth integral
						
						double effectiveArea=*index_multi(*effectiveAreas[y][p],indices);
						effectiveArea*=1.0e4; //convert m^2 to cm^2
						expectation[i] += effectiveArea * fluxIntegral * *livetimes[y];
						
						//We can also compute the uncertainty due to limited simulation statistics.
						//This requires adding the error terms in quadrature.
						double effectiveAreaErr=*index_multi(*effectiveAreaErrs[y][p],indices);
						effectiveAreaErr*=1.0e4; //convert m^2 to cm^2
						effectiveAreaErr*=fluxIntegral * *livetimes[y];
						effectiveAreaErr*=effectiveAreaErr; //square!
						expectationUncertainty[i] += effectiveAreaErr;
					}
				}
			}
		}
	}
	
	//Print out the expected number of events in each energy proxy bin, with
	//statistical uncertainties.
	printf("Astrophysical flux expectation as a function of energy proxy\n");
	double total=0;
	for(unsigned int i=0; i<energyProxyBins; i++){
		printf(" %lf: %lf [%lf,%lf]\n",energyProxyEdges[i],expectation[i],
		       expectation[i]-sqrt(expectationUncertainty[i]),
		       expectation[i]+sqrt(expectationUncertainty[i]));
		total+=expectation[i];
	}
	printf("Total: %lf\n",total);
//-----------end effec_area_demo_insertion
*/

		Verosimilitud::Verosimilitud(unsigned int my_numneu)
		{
			numneu=my_numneu;	
			std::cout << "BLEAT" << std::endl;
		}


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
        	&return {LogGaussianProbability(np[0],np[0]_prior), ... }
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

//void Verosimilitud::testhdf(void)
int main(void)
{
	EffectiveArea ea = EffectiveArea();

	return 0;
}
