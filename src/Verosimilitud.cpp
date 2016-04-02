#include "Verosimilitud.h"
#include <iostream>
#include <stdexcept>
#include <Python.h>
#include "Flux.h"
#include "Tensor.h"
#include "H5Cpp.h"
#include "EffectiveArea.h"
#include "ConventionalFlux.h"
#include "ICData.h"
#include <cmath>
#include <cfloat>

#include <dlib/optimization.h>
#include <dlib/member_function_pointer.h>


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
	
	//the expected distribution of astrophysical events
	//in the energy proxy observable
	double* expectation=(double*)malloc(energyProxyBins*sizeof(double));
	double* expectationUncertainty=(double*)malloc(energyProxyBins*sizeof(double));
	
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



	eff_area = new EffectiveArea();

	conv_flux = new ConventionalFlux();

    unsigned int edge_indices[4] = {0,0,0,0};
	NeutrinoEnergyEdges = eff_area->GetEdge(edge_indices);
    edge_indices[3] = 1;
    CosZenithEdges = eff_area->GetEdge(edge_indices);
    edge_indices[3] = 2;
    EnergyProxyEdges = eff_area->GetEdge(edge_indices);
	unsigned int ind1 = 6;
	unsigned int ind2 = EnergyProxyBins-17;
	std::cout << "EDGE 6 " << EnergyProxyEdges->Index(&ind1) << std::endl;
	std::cout << "EDGE FINAL  " << EnergyProxyEdges->Index(&ind2) << std::endl;
	icd = new ICData(EnergyProxyEdges, CosZenithEdges);
	icd->OpenCSV("observed_events.dat");
	//std::vector<unsigned int> data_cosz(CosZenithBins,0);
	//std::vector<unsigned int> data_eprox(EnergyProxyBins,0);
	//icd.BinData(&data_cosz,&data_eprox);

	icd->ReadCSV();
	unsigned int datadims[2];
	datadims[0]=EnergyProxyBins;
	datadims[1]=CosZenithBins;
	data = new Tensor(2,datadims,0);
    icd->BinData(data);

	eprox_cuts = new std::vector<double>(2,0);
	cosz_cuts = new std::vector<double>(2,0);

	(*eprox_cuts)[0]=0;
	(*eprox_cuts)[1]=EnergyProxyBins;

	(*cosz_cuts)[0]=0;
	(*cosz_cuts)[1]=CosZenithBins;


	//calculate bin centers

	eprox_centers = new std::vector<double>(EnergyProxyBins,0);
	coszenith_centers = new std::vector<double>(CosZenithBins,0);
	unsigned int ip1;

	for (unsigned int i=0; i<EnergyProxyBins; i++)
	{
		ip1=i+1;
		(*eprox_centers)[i] = (EnergyProxyEdges->Index(&i)+EnergyProxyEdges->Index(&ip1))/2.0;
	}

	for (unsigned int i; i<CosZenithBins; i++)
	{
		ip1=i+1;
		(*coszenith_centers)[i] = (CosZenithEdges->Index(&i)+CosZenithEdges->Index(&ip1))/2.0;
	}
	



	}


	void Verosimilitud::SetCoszCuts(std::vector<double> cuts)
	{
		if(cuts.size() != (*cosz_cuts).size()) 
		{
			throw std::runtime_error("Cuts Size Incorrect.");
		}

		for (unsigned int i; i<cuts.size(); i++)
		{
			(*cosz_cuts)[i]=cuts[i];
		}
	}

	void Verosimilitud::SetEproxCuts(std::vector<double> cuts)
	{
		if(cuts.size() != (*eprox_cuts).size()) 
		{
			throw std::runtime_error("Cuts Size Incorrect.");
		}

		for (unsigned int i; i<cuts.size(); i++)
		{
			(*eprox_cuts)[i]=cuts[i];
		}
	}


		Verosimilitud::~Verosimilitud()
		{
			delete eff_area;
			delete conv_flux;
			delete data;
			delete icd;
			delete eprox_cuts;
			delete cosz_cuts;
			delete eprox_centers;
			delete coszenith_centers;
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


std::vector<unsigned int> Verosimilitud::GetExpDims()
{
	return exp_dimvec;
}

std::vector<double> Verosimilitud::GetEproxEdges()
{
	return *eprox_edges;
}

std::vector<double> Verosimilitud::GetCosZenithEdges()
{
	return *coszenith_edges;
}

std::vector<double> Verosimilitud::CalculateExpectation()
{
	unsigned int dyn_indices[3];
	unsigned int exp_indices[2];
	unsigned int area_indices[3];
	unsigned int edge_indices[4];


	unsigned int exp_dims[2]={EnergyProxyBins,CosZenithBins};
	expectation = new Tensor(2,exp_dims,0);

	unsigned int counter=0;

	double countsheep=0;

	std::cout << "AtLoops" << std::endl;
	
	//Loop over energy proxies.
	for(unsigned int ep=0; ep<EnergyProxyBins; ep++)
	{
		dyn_indices[2]=ep;
	
		//Loop over years 2010, 2011
		for(unsigned int year=0; year<2; year++)
		{
			//Do muon neutrinos only: no loop over flavor.
			unsigned int flavor=0;

			double livetime= eff_area->GetLivetime(year);			

			//Loop over matter/antimatter
			for (unsigned int anti=0; anti<2; anti++)
			{

				unsigned int area_indices[3] = {year,flavor,anti};
				area = eff_area->GetArea(area_indices);

				detcorr = conv_flux->GetDetCorr(&year);
				flux = conv_flux->GetFlux(&anti);


				//Loop over true energy bins.
				for(unsigned int e=0; e<NeutrinoEnergyBins; e++)
				{
					dyn_indices[0]=e;
					unsigned int ep1=e+1;	
					double eMin=NeutrinoEnergyEdges->Index(&e);
					double eMax=NeutrinoEnergyEdges->Index(&ep1);

					//Loop over coszenith bins.
					for(unsigned int z=0; z<CosZenithBins; z++)
					{
						counter++;
						//std::cout << counter << std::endl;
						dyn_indices[1]=z;
						unsigned int zp1=z+1;	
						double CosZenithMin=CosZenithEdges->Index(&z);
						double CosZenithMax=CosZenithEdges->Index(&zp1);

						double FluxIntegral = flux->Index(dyn_indices);
						//implement DOM correction?
						double DomCorr = detcorr->Index(dyn_indices);
						double EffAreaVal = area->Index(dyn_indices); 
						EffAreaVal*=1.0e4; //m^2 to cm^2

						double eavg = (eMin+eMax)/2.0;
						double zavg = (CosZenithMin+CosZenithMax)/2.0;
						//double oscprob = OscillationProbability(eavg,zavg);
	
						unsigned int exp_indices[2]={ep,z};

						double last=expectation->Index(exp_indices);
						double current=last+FluxIntegral*EffAreaVal*livetime*DomCorr;

						expectation->SetIndex(exp_indices,current);
						countsheep+=FluxIntegral*EffAreaVal*livetime*DomCorr;
	
					}
				}
			}
		}
	}
	std::cout<<"COUNTSHEEP: " << countsheep << std::endl;

	double* eprox_edges_array =EnergyProxyEdges->GetDataPointer();
	unsigned int eprox_size = EnergyProxyEdges->GetDataLength();
	double* coszenith_edges_array = CosZenithEdges->GetDataPointer();
	unsigned int coszenith_size= CosZenithEdges->GetDataLength();

	//eprox_edges=new std::vector<double>(eprox_edges_array, eprox_edges_array+eprox_size);
	//coszenith_edges = new std::vector<double>(coszenith_edges_array, coszenith_edges_array+coszenith_size);


	double * exparray = expectation->GetDataPointer();
	unsigned int explen = expectation->GetDataLength();
	std::vector<double> expvec(exparray, exparray + explen);
	unsigned int expdims[]={0,0};
	expectation->GetDims(expdims);
	exp_dimvec.push_back(expdims[0]);
	exp_dimvec.push_back(expdims[1]);
	//std::cout << exp_dims[0] << std::endl;
	std::cout << expdims[0] << std::endl;
	std::cout << expdims[1] << std::endl;
	std::cout << exp_dimvec[0] << std::endl;
	std::cout << exp_dimvec[1] << std::endl;
	return expvec;

}




double Verosimilitud::Chi2MinNuisance(std::vector<double>* param)
{
//Minimize over nuisance parameters. 

  dlib::matrix<double,0,1> nuisance(2);
  nuisance(0)=(*param)[0];
  nuisance(1)=(*param)[1];


  dlib::find_min(dlib::bfgs_search_strategy(), dlib::objective_delta_stop_strategy(1e-7),Chi2_caller(this),Chi2grad_caller(this),nuisance,-1);

  (*param)[0]=nuisance(0);
  (*param)[1]=nuisance(1);
	return Chi2(nuisance);
}





double Verosimilitud::Chi2(const dlib::matrix<double,0,1>& nuisance)
{
	unsigned int indices[2];
	double scalar_exp;
	double scalar_data;
	double llh=0;
	double sllh=0;
	int count=0;
	double prob;
	double pert_scalar_exp;



	const double norm = nuisance(0);
	const double gamma = nuisance(1);


	//calculate unsaturated log-likelihood with nuisance-perturbed expectation

	for(unsigned int ep=(*eprox_cuts)[0]; ep<(*eprox_cuts)[1]; ep++)
   	{
		for(unsigned int z=(*cosz_cuts)[0]; z<(*cosz_cuts)[1]; z++)
		{
			indices[0]=ep;
			indices[1]=z;
            scalar_exp=expectation->Index(indices);
			scalar_data=data->Index(indices);
			pert_scalar_exp=norm*scalar_exp*pow((*eprox_centers)[ep]/34592.0,gamma); // fixme what is 34592 number?
			prob = LogPoissonProbability(scalar_data,pert_scalar_exp);
			if (std::isnan(prob))
			{
				prob=0;
			}
			llh+=prob;
			count++;
			//std::cout <<"llh: "<< llh << std::endl;
			//std::cout << "count: " << count << " / " << ((*eprox_cuts)[1]-(*eprox_cuts)[0])*((*cosz_cuts)[1]-(*cosz_cuts)[0]) << std::endl;
		}
	}	

	llh+=LogGaussianProbability(norm,norm_mean,norm_sigma);
	llh+=LogGaussianProbability(gamma,gamma_mean,gamma_sigma);


	//Calculate saturated log-likelihood

	for(unsigned int ep=(*eprox_cuts)[0]; ep<(*eprox_cuts)[1]; ep++)
   	{
		for(unsigned int z=(*cosz_cuts)[0]; z<(*cosz_cuts)[1]; z++)
		{
			indices[0]=ep;
			indices[1]=z;
			scalar_data=data->Index(indices);
			double satprob;
			satprob=LogPoissonProbability(scalar_data,scalar_data);
			if (std::isnan(satprob))
			{
				satprob=0;
			}
			sllh+=satprob;

		}
	}	
	
				
	return 2*(sllh-llh);
}









dlib::matrix<double,0,1> Verosimilitud::Chi2Gradient(const dlib::matrix<double,0,1>& nuisance)
{
	unsigned int indices[2];
	double scalar_exp;
	double scalar_data;
	double pert_scalar_exp;

	const double norm = nuisance(0);
	const double gamma = nuisance(1);

	double center;
	
	double grad0=0;
	double grad1=0;


	double strange_constant=34592.0;

	for(unsigned int ep=(*eprox_cuts)[0]; ep<(*eprox_cuts)[1]; ep++)
   	{
		center=(*eprox_centers)[ep];
		for(unsigned int z=(*cosz_cuts)[0]; z<(*cosz_cuts)[1]; z++)
		{
			indices[0]=ep;
			indices[1]=z;
            scalar_exp=expectation->Index(indices);
			scalar_data=data->Index(indices);
			pert_scalar_exp=norm*scalar_exp*pow(center/strange_constant,gamma); // fixme what is 34592 number?
			if(!(std::isnan(LogPoissonProbability(scalar_data,pert_scalar_exp))))
			{
				grad0+=(1- scalar_data/pert_scalar_exp)*scalar_exp*pow(center/strange_constant,gamma);
				grad1+=(1- scalar_data/pert_scalar_exp)*norm*scalar_exp*log(center/strange_constant)*pow(center/strange_constant,gamma);
			}

		}
	}	

	grad0+=(norm-norm_mean)/pow(norm_sigma,2);
	grad1+=(gamma-gamma_mean)/pow(gamma_sigma,2);

	grad0*=2;
	grad1*=2;


	dlib::matrix<double,0,1> gradient(2);

	gradient(0) = grad0;
	gradient(1) = grad1;

	return gradient;
}








int main(void)
{

	Verosimilitud v(3);

	std::vector<double> my_eprox_cuts(2,0);
	std::vector<double> my_cosz_cuts(2,0);
	
	my_eprox_cuts[0]=6;
	my_eprox_cuts[1]=v.EnergyProxyBins-16;

	v.SetEproxCuts(my_eprox_cuts);

	v.CalculateExpectation();


	std::vector<double> nuisance_init(2,0);
	nuisance_init[0]=1; //set norm to one, leave gamma at zero
	nuisance_init[1]=0;
	double ret = v.Chi2MinNuisance(&nuisance_init);	


	std::cout << "Chi2: " << ret << std::endl;
	std::cout << "Norm: " << nuisance_init[0] << std::endl;
	std::cout << "Gamma: " << nuisance_init[1] << std::endl;

	return 0;
}
