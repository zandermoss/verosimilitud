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
#include <math.h>
#include <string>

#include <dlib/optimization.h>
#include <dlib/member_function_pointer.h>


Verosimilitud::Verosimilitud(unsigned int my_numneu,unsigned int loyear, unsigned int hiyear)
{
	numneu=my_numneu;	


	//Which years do we want to include?
	data_years[0]=loyear;
	data_years[1]=hiyear;

	simps_nintervals=2;

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


	//Bin data from 2010, 2011, or both!
	std::string fnames[2] = {"2010.dat","2011.dat"};

	unsigned int datadims[2];
	datadims[0]=EnergyProxyBins;
	datadims[1]=CosZenithBins;
	data = new Tensor(2,datadims,0);

	for (unsigned int y=data_years[0]; y<data_years[1]; y++)
	{
		icd[y] = new ICData(EnergyProxyEdges, CosZenithEdges);
		icd[y]->OpenCSV(fnames[y]);
		icd[y]->ReadCSV();
    	icd[y]->BinData(data);
	}



	//Apply the appropriate cuts
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


	


	void Verosimilitud::SetSimpsNIntervals(int nintervals)
	{
		simps_nintervals=nintervals;
	}

	void Verosimilitud::SetCoszCuts(std::vector<double> cuts)
	{
		/*
		if(cuts.size() != (*cosz_cuts).size()) 
		{
			throw std::runtime_error("Cuts Size Incorrect.");
		}
		*/

		for (unsigned int i; i<cuts.size(); i++)
		{
			(*cosz_cuts)[i]=cuts[i];
			std::cout <<"COSZ CUTS: " << (*cosz_cuts)[i] << std::endl;
		}
	}

	void Verosimilitud::SetEproxCuts(std::vector<double> cuts)
	{
		std::cout << "BEGINNING OF EPROXCUT" << std::endl;

		/*
		if(cuts.size() != (*eprox_cuts).size()) 
		{
			throw std::runtime_error("Cuts Size Incorrect.");
		}
		*/

		for (unsigned int i; i<cuts.size(); i++)
		{
			(*eprox_cuts)[i]=cuts[i];
			std::cout <<"EPROX CUTS: "  << (*eprox_cuts)[i] << std::endl;
		}
	}


		Verosimilitud::~Verosimilitud()
		{
			delete eff_area;
			delete conv_flux;
			delete data;
			delete icd[0];
			delete icd[1];
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

	    double Verosimilitud::OscillationProbability(double energy,double zenith, double anti)
		{
	   	 	//Here we call nusheep. Although, would it be better to precalculate an array?
			//I think for now I will implement direct calls.
			std::vector<double> argument;
			argument.push_back(energy);
			argument.push_back(zenith);
			argument.push_back(anti);
			double osc_prob = de_solver(argument,user_data);	
			std::cout << "  ZENITH: " << zenith << "  ENERGY: " << energy << " ANTI: " << anti  << "  PROB: " << osc_prob << std::endl;				
			return osc_prob;
   		}

		void Verosimilitud::SetDeSolver(pyoscfunc my_de_solver, void* my_user_data)
		{
			de_solver=my_de_solver;
			user_data=my_user_data;		
		}





std::vector<unsigned int> Verosimilitud::GetExpDims()
{
	return exp_dimvec;
}

std::vector<unsigned int> Verosimilitud::GetDataDims()
{
	return dat_dimvec;
}

std::vector<double> Verosimilitud::GetEnergyEdges()
{
	return *energy_edges;
}

std::vector<double> Verosimilitud::GetEproxEdges()
{
	return *eprox_edges;
}

std::vector<double> Verosimilitud::GetCosZenithEdges()
{
	return *coszenith_edges;
}


double Verosimilitud::TestFunction(double e, double z, double a)
{
	return log(pow(3*e,4)*pow(z,4));
}

double Verosimilitud::SimpsAvg(double coszmin,double coszmax, double emin, double emax, double anti, int nintervals_e)
{
	double width = (emax-emin)/(double)nintervals_e;
	double integral=0;
	double mean=0;
	double coszval;
	for(int j=0; j<2; j++)
	{
		coszval = coszmin+(double)j*(coszmax-coszmin);
		integral+=OscillationProbability(emin,acos(coszval),anti);
		for (int i=1; i<(nintervals_e/2); i++)
		{
			integral+=2*OscillationProbability(emin+2*(double)i*width,acos(coszval),anti);
		}
		for (int i=1; i<(nintervals_e/2+1); i++)
		{
			integral+=4*OscillationProbability(emin+(2*(double)i-1)*width,acos(coszval),anti);
		}
		integral+=OscillationProbability(emax,acos(coszval),anti);
		integral*=width/3.0;
		mean+=integral/(2*(emax-emin));
		integral=0;
	}	
	return mean;		

}
/*
double Verosimilitud::ESimpsAvg(double coszval, double emin, double emax, double anti, int nintervals_e)
{
	double width = (emax-emin)/(double)nintervals;
	double integral=0;

	integral+=OscillationProbability(emin,coszval,anti);
	for (int i=1; i<(nintervals/2 -1); i++)
	{
		integral+=2*OscillationProbability(emin+2*(double)i*width,coszval,anti);
	}
	for (int i=1; i<(nintervals/2); i++)
	{
		integral+=4*OscillationProbability(emin+(2*(double)i-1)*width,coszval,anti);
	}
	integral+=OscillationProbability(emax,coszval,anti);
	integral*=width/3.0;
	
	return integral		

}

double Verosimilitud::CosSimpsAvg(double coszmin, double coszmax, double emin, double emax, double anti, int nintervals_cosz, int nintervals_e)
{
	double width = (coszmax-coszmin)/(double)nintervals_cosz;
	double integral=0;
	double mean;

	integral+=ESimpsAvg(coszmin,emin,emax,anti,nintervals_e);
	for (int i=1; i<(nintervals_cosz/2 -1); i++)
	{
		integral+=2*ESimpsAvg(coszmin+2*(double)i*width,emin,emax,anti,nintervals_e);
	}
	for (int i=1; i<(nintervals_cosz/2); i++)
	{
		integral+=4*ESimpsAvg(coszmin+(2*(double)i-1)*width,emin,emax,anti,nintervals_e);
	}
	integral+=ESimpsAvg(coszmax,emin,emax,anti,nintervals_e);
	integral*=width/3.0;
	mean=integral/(coszmax-coszmin);
	return mean;
			
}
*/






void Verosimilitud::CalculateExpectation()
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
	
	//Loop over true energies.
	for(unsigned int e=0; e<NeutrinoEnergyBins; e++)
	{
		std::cout << "TRUE ENERGY: " << e << std::endl;
		dyn_indices[0]=e;
		unsigned int ep1=e+1;	
		double eMin=NeutrinoEnergyEdges->Index(&e);
		double eMax=NeutrinoEnergyEdges->Index(&ep1);
		double eavg = (eMin+eMax)/2.0;
	

		//Loop over matter/antimatter
		for (unsigned int anti=0; anti<2; anti++)
		{

			//Loop over coszenith bins.
			for(unsigned int z=0; z<CosZenithBins; z++)
			{
				dyn_indices[1]=z;
				unsigned int zp1=z+1;	
				double CosZenithMin=CosZenithEdges->Index(&z);
				double CosZenithMax=CosZenithEdges->Index(&zp1);
				double zavg = (CosZenithMin+CosZenithMax)/2.0;

				double oscprob = OscillationProbability(eavg,acos(zavg),(double)anti);
				//double oscprob = SimpsAvg(CosZenithMin,CosZenithMax,eMin,eMax,(double)anti, simps_nintervals);


				//Loop over years 2010, 2011
				for(unsigned int year=data_years[0]; year<data_years[1]; year++)
				{
					//Do muon neutrinos only: no loop over flavor.
					unsigned int flavor=0;

					double livetime= eff_area->GetLivetime(year);			

					unsigned int area_indices[3] = {year,flavor,anti};
					area = eff_area->GetArea(area_indices);

					detcorr = conv_flux->GetDetCorr(&year);
					flux = conv_flux->GetFlux(&anti);



					//Loop over energy proxy bins.
					for(unsigned int ep=0; ep<EnergyProxyBins; ep++)
					{
						dyn_indices[2]=ep;

						counter++;
						//std::cout << counter << std::endl;

						double FluxIntegral = flux->Index(dyn_indices);
						//implement DOM correction?
						double DomCorr = detcorr->Index(dyn_indices);
						double EffAreaVal = area->Index(dyn_indices); 
						EffAreaVal*=1.0e4; //m^2 to cm^2

	
						unsigned int exp_indices[2]={ep,z};

						double last=expectation->Index(exp_indices);
						double current=last+FluxIntegral*EffAreaVal*livetime*DomCorr*oscprob;
						//double current=last+FluxIntegral*EffAreaVal*livetime*DomCorr;

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
	double* energy_edges_array =NeutrinoEnergyEdges->GetDataPointer();
	unsigned int energy_size = NeutrinoEnergyEdges->GetDataLength();
	double* coszenith_edges_array = CosZenithEdges->GetDataPointer();
	unsigned int coszenith_size= CosZenithEdges->GetDataLength();

	eprox_edges=new std::vector<double>(eprox_edges_array, eprox_edges_array+eprox_size);
	energy_edges=new std::vector<double>(energy_edges_array, energy_edges_array+energy_size);
	coszenith_edges = new std::vector<double>(coszenith_edges_array, coszenith_edges_array+coszenith_size);

/*
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
*/

}

std::vector<double> Verosimilitud::GetExpectationVec(void)
{
	double * exparray = expectation->GetDataPointer();
	unsigned int explen = expectation->GetDataLength();
	std::vector<double> expvec(exparray, exparray + explen);
	unsigned int expdims[]={0,0};
	expectation->GetDims(expdims);
	exp_dimvec.push_back(expdims[0]);
	exp_dimvec.push_back(expdims[1]);
	return expvec;
}


std::vector<double> Verosimilitud::GetFluxVec(void)
{

	unsigned int anti=0;
    Tensor* my_flux = conv_flux->GetFlux(&anti);

	double * fluxarray = my_flux->GetDataPointer();
	unsigned int fluxlen = my_flux->GetDataLength();
	std::cout << "FLUXLEN: " << fluxlen << std::endl;
	std::vector<double> fluxvec(fluxlen,0);
	for (int i=0; i<fluxlen; i++)
	{
		fluxvec[i]=fluxarray[i];
	}

	return fluxvec;
}

std::vector<double> Verosimilitud::GetAreaVec(void)
{
	unsigned int area_indices[3]={0,0,0};
	Tensor* my_area = eff_area->GetArea(area_indices);
	double * areaarray = my_area->GetDataPointer();
	unsigned int arealen = my_area->GetDataLength();
	std::cout << "AREALEN: " << arealen << std::endl;
	std::vector<double> areavec(arealen,0);
	for (int i=0; i<arealen; i++)
	{
		areavec[i]=areaarray[i];
	}
	return areavec;
}


std::vector<double> Verosimilitud::GetPertExpectationVec(std::vector<double> nuisance)
{
	unsigned int indices[2];
	double scalar_exp;
	double pert_scalar_exp;

	double norm=nuisance[0];
	double gamma=nuisance[1];


    unsigned int exp_dims[2]={EnergyProxyBins,CosZenithBins};
    Tensor* perturbed_expectation = new Tensor(2,exp_dims,0);



   for(unsigned int ep=0; ep<EnergyProxyBins; ep++)
    {
        for(unsigned int z=0; z<CosZenithBins; z++)
        {
            indices[0]=ep;
            indices[1]=z;
            scalar_exp=expectation->Index(indices);
            pert_scalar_exp=norm*scalar_exp*pow((*eprox_centers)[ep]/34592.0,gamma); // fixme what is 34592 number?
			perturbed_expectation->SetIndex(indices,pert_scalar_exp);
        }
    }




	double * exparray = perturbed_expectation->GetDataPointer();
	unsigned int explen = perturbed_expectation->GetDataLength();
	std::vector<double> expvec(exparray, exparray + explen);
	unsigned int expdims[]={0,0};
	perturbed_expectation->GetDims(expdims);
	exp_dimvec.push_back(expdims[0]);
	exp_dimvec.push_back(expdims[1]);
	delete perturbed_expectation;
	return expvec;
}



std::vector<double> Verosimilitud::GetDataVec(void)
{
	double * datarray = data->GetDataPointer();
	unsigned int datlen = data->GetDataLength();
	std::vector<double> datvec(datarray, datarray + datlen);
	unsigned int datdims[]={0,0};
	data->GetDims(datdims);
	dat_dimvec.push_back(datdims[0]);
	dat_dimvec.push_back(datdims[1]);
	return datvec;
}

/*
std::vector<double> Verosimilitud::GetDataVec(double scale)
{
   unsigned int indices[2];
    double pert_scalar_data;
    double scalar_data;


    unsigned int exp_dims[2]={EnergyProxyBins,CosZenithBins};
    Tensor* perturbed_data = new Tensor(2,exp_dims,0);


   
   for(unsigned int ep=0; ep<EnergyProxyBins; ep++)
    {
        for(unsigned int z=0; z<CosZenithBins; z++)
        {
            indices[0]=ep;
            indices[1]=z;
            scalar_data=data->Index(indices);
            pert_scalar_data=scale*scalar_data; // fixme what is 34592 number?
            perturbed_data->SetIndex(indices,pert_scalar_data);
        }
    }



	double * datarray = perturbed_data->GetDataPointer();
	unsigned int datlen = perturbed_data->GetDataLength();
	std::vector<double> datvec(datarray, datarray + datlen);
	unsigned int datdims[]={0,0};
	perturbed_data->GetDims(datdims);
	dat_dimvec.push_back(datdims[0]);
	dat_dimvec.push_back(datdims[1]);
	delete perturbed_data;

	return datvec;

}
*/






std::vector<double> Verosimilitud::Chi2MinNuisance(std::vector<double> param)
{
//Minimize over nuisance parameters. 

  dlib::matrix<double,0,1> nuisance(2);
  nuisance(0)=(param)[0];
  nuisance(1)=(param)[1];

	
  dlib::matrix<double,0,1> lo_bounds(2);
  dlib::matrix<double,0,1> hi_bounds(2);
  lo_bounds(0)=0.01;
  lo_bounds(1)=-2.0;
  hi_bounds(0)=2.0;
  hi_bounds(1)=2.0;


  dlib::find_min_box_constrained(dlib::bfgs_search_strategy(), dlib::objective_delta_stop_strategy(1e-7),Chi2_caller(this),Chi2grad_caller(this),nuisance,lo_bounds,hi_bounds);


	std::vector<double> ret(3,0);
	ret[0]=Chi2(nuisance);
	ret[1]=nuisance(0);
	ret[2]=nuisance(1);

	return ret;
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

	std::cout << "EPROX CUTS: " <<(*eprox_cuts)[0] << "    "  << (*eprox_cuts)[1] << std::endl;

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
		}
	}	


	double norm_penalty = LogGaussianProbability(norm,norm_mean,norm_sigma);
	double gamma_penalty= LogGaussianProbability(gamma,gamma_mean,gamma_sigma);
	
	std::cout << "NORM PENALTY  " << norm_penalty << std::endl;
	std::cout << "GAMMA PENALTY  " << gamma_penalty << std::endl;
	std::cout << "GAMMA  " << gamma << std::endl;
	std::cout << "GAMMAmu  " << gamma_mean << std::endl;
	std::cout << "GAMMAsig  " << gamma_sigma << std::endl;


	llh+=norm_penalty + gamma_penalty;


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
	

	std::cout << "CHI2: " << 2*(sllh-llh) << std::endl;
	std::cout << std::endl;
				
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
/*
	double step=0.1;
	for (int x=0; x<10; x++)
	{
		std::cout << LogGaussianProbability(1-x*step,1,0.4) << std::endl;
	}
	for (int x=0; x<10; x++)
	{
		std::cout << LogGaussianProbability(1+x*step,1,0.4) << std::endl;
	}
	Verosimilitud v(3);
	double last=0;
	for (int x=1; x<10; x++)
	{
		//std::cout << "N: " << 2*x << " AVG: " <<	v.SimpsAvg(1.57,3.14,100,106,0,2*x) << std::endl;
		double val=v.SimpsAvg(1,2,1,2,0,2*x);
		std::cout << "N: " << 2*x << " AVG: " <<	fabs(val-last) << std::endl;
		std::cout << "N: " << 2*x << " AVG: " <<	val << std::endl;
		last=val;
	}
*/

/*
	std::vector<double> my_eprox_cuts(2,0);
	std::vector<double> my_cosz_cuts(2,0);
	
	my_eprox_cuts[0]=6;
	my_eprox_cuts[1]=v.EnergyProxyBins-16;

	v.SetEproxCuts(my_eprox_cuts);

	v.CalculateExpectation();


	std::vector<double> nuisance_init(2,0);
	nuisance_init[0]=1; //set norm to one, leave gamma at zero
	nuisance_init[1]=0;
	std::vector<double> ret = v.Chi2MinNuisance(nuisance_init);	


	std::cout << "Chi2: " << ret[0] << std::endl;
	std::cout << "Norm: " << ret[1] << std::endl;
	std::cout << "Gamma: " << ret[2] << std::endl;
*/
	return 0;
}
