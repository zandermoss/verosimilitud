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

	for (unsigned int y=0; y<2; y++)
	{
		icd[y] = new ICData(EnergyProxyEdges, CosZenithEdges);
	}

	for (unsigned int y=data_years[0]; y<data_years[1]; y++)
	{
		icd[y]->OpenCSV(fnames[y].c_str());
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

	


	void Verosimilitud::SetSimpsNIntervals(int nintervals)
	{
		simps_nintervals=nintervals;
	}

	void Verosimilitud::SetCoszCuts(std::vector<double> cuts)
	{
		for (unsigned int i; i<cuts.size(); i++)
		{
			(*cosz_cuts)[i]=cuts[i];
			std::cout <<"COSZ CUTS: " << (*cosz_cuts)[i] << std::endl;
		}
	}

	void Verosimilitud::SetEproxCuts(std::vector<double> cuts)
	{

		for (unsigned int i; i<cuts.size(); i++)
		{
			(*eprox_cuts)[i]=cuts[i];
			std::cout <<"EPROX CUTS: "  << (*eprox_cuts)[i] << std::endl;
		}
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


double Verosimilitud::TestFunction(double e, double z)
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
	
	unsigned int which_flux;

	unsigned int exp_dims[2]={EnergyProxyBins,CosZenithBins};

//	expectation = new Tensor(2,exp_dims,0);

	Tensor* expectation[4];	
	for (int i = 0; i<4; i++){
		expectation[i]= new Tensor(2,exp_dims,0);
		}

	unsigned int counter=0;

	std::cout << "AtLoops" << std::endl;
	
	//Loop over true energies.
	for(unsigned int e=0; e<NeutrinoEnergyBins; e++)
	{
//		std::cout << "TRUE ENERGY: " << e << std::endl;
		dyn_indices[0]=e;
		unsigned int ep1=e+1;	
		double eMin=NeutrinoEnergyEdges->Index(&e);
		double eMax=NeutrinoEnergyEdges->Index(&ep1);
		double eavg = (eMin+eMax)/2.0;
	
		//Loop over mesons
		for (unsigned int meson=0; meson<2; meson++)
		{
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

					//double oscprob = SimpsAvg(CosZenithMin,CosZenithMax,eMin,eMax,(double)anti, simps_nintervals);
					//uncomment above to let neutrinos oscillate

					//Loop over years 2010, 2011
					for(unsigned int year=data_years[0]; year<data_years[1]; year++)
					{
//						std::cout << "YEAR: " << year << std::endl;
						//Do muon neutrinos only: no loop over flavor.
						unsigned int flavor=0;

						double livetime= eff_area->GetLivetime(year);			

						unsigned int area_indices[3] = {year,flavor,anti};
						area = eff_area->GetArea(area_indices);

						detcorr = conv_flux->GetDetCorr(&year);
						
						if 		(meson == 0 && anti == 0){ which_flux = 0; }
						else if (meson == 0 && anti == 1){ which_flux = 1; }
						else if (meson == 1 && anti == 0){ which_flux = 2; }
						else if (meson == 1 && anti == 1){ which_flux = 3; }
						
						flux = conv_flux->GetFlux(&which_flux);

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

							double last=expectation[which_flux]->Index(exp_indices);
							//double current=last+FluxIntegral*EffAreaVal*livetime*DomCorr*oscprob;
							double current=last+FluxIntegral*EffAreaVal*livetime*DomCorr;

							expectation[which_flux]->SetIndex(exp_indices,current); //make two separate expectations ->make expectation an array
	
						}
					}
				}
			}
		}
	}

	double* eprox_edges_array =EnergyProxyEdges->GetDataPointer();
	unsigned int eprox_size = EnergyProxyEdges->GetDataLength();
	double* energy_edges_array =NeutrinoEnergyEdges->GetDataPointer();
	unsigned int energy_size = NeutrinoEnergyEdges->GetDataLength();
	double* coszenith_edges_array = CosZenithEdges->GetDataPointer();
	unsigned int coszenith_size= CosZenithEdges->GetDataLength();

	eprox_edges=new std::vector<double>(eprox_edges_array, eprox_edges_array+eprox_size);
	energy_edges=new std::vector<double>(energy_edges_array, energy_edges_array+energy_size);
	coszenith_edges = new std::vector<double>(coszenith_edges_array, coszenith_edges_array+coszenith_size);

}

/*
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
*/

std::vector<double> Verosimilitud::GetFluxVec(void)
// Why is this hard-coded for anti=0?
// This is not going to work with the separate Pi, K fluxes
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
	double scalar_exp[4];
	double pert_scalar_exp;

	double norm=nuisance[0];
	double gamma=nuisance[1];
	double r_kpi=nuisance[2];
	double r_nubarnu=nuisance[3];

	double E0 = 34592.0;
	double center;

    unsigned int exp_dims[2]={EnergyProxyBins,CosZenithBins};
    Tensor* perturbed_expectation = new Tensor(2,exp_dims,0);



   for(unsigned int ep=0; ep<EnergyProxyBins; ep++)
    {
    	center = (*eprox_centers)[ep];
    	
        for(unsigned int z=0; z<CosZenithBins; z++)
        {
            indices[0]=ep;
            indices[1]=z;
            
            for (unsigned int i=0; i<4; i++){
            	scalar_exp[i]=expectation[i].Index(indices);
            }
            
            pert_scalar_exp = norm*pow(center/E0,-gamma)*(scalar_exp[0]+r_kpi*scalar_exp[2]+r_nubarnu*(scalar_exp[1]+r_kpi*scalar_exp[3]));
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







std::vector<double> Verosimilitud::Chi2MinNuisance(std::vector<double> param)
{
//Minimize over nuisance parameters. 
  dlib::matrix<double,0,1> nuisance(4);
  nuisance(0)=(param)[0];
  nuisance(1)=(param)[1];
  nuisance(2)=(param)[2];
  nuisance(3)=(param)[3];

	
  dlib::matrix<double,0,1> lo_bounds(4);
  dlib::matrix<double,0,1> hi_bounds(4);
  lo_bounds(0)=0.01;
  lo_bounds(1)=-2.0;
  lo_bounds(2)=0.01;
  lo_bounds(3)=0.01;
  hi_bounds(0)=2.0;
  hi_bounds(1)=2.0;
  hi_bounds(2)=2.0;
  hi_bounds(3)=2.0;


	// Marjon: error is in the line below
  dlib::find_min_box_constrained(dlib::bfgs_search_strategy(), dlib::objective_delta_stop_strategy(1e-7),Chi2_caller(this),Chi2grad_caller(this),nuisance,lo_bounds,hi_bounds);

	std::vector<double> ret(5,0);
	ret[0]=Chi2(nuisance);
	ret[1]=nuisance(0);
	ret[2]=nuisance(1);
	ret[3]=nuisance(2);
	ret[4]=nuisance(3);

	return ret;
}





double Verosimilitud::Chi2(const dlib::matrix<double,0,1>& nuisance)
{
	unsigned int indices[2];
	double scalar_exp[4];
	double scalar_data;
	double llh=0;
	double sllh=0;
	int count=0;
	double prob;
	double tot_scalar_exp;
	double raw_scalar_data;

	const double norm = nuisance(0);
	const double gamma = nuisance(1);
	const double r_kpi = nuisance(2);		//ratio of kaons to pions
	const double r_nubarnu = nuisance(3);	//ratio of nubars to nus
	
	double center;
	double E0 = 34592.0;

	//calculate unsaturated log-likelihood with nuisance-perturbed expectation

//	std::cout << "EPROX CUTS: " <<(*eprox_cuts)[0] << "    "  << (*eprox_cuts)[1] << std::endl;

	for(unsigned int ep=(*eprox_cuts)[0]; ep<(*eprox_cuts)[1]; ep++)
   	{
   	
   	center=(*eprox_centers)[ep];
   	
		for(unsigned int z=(*cosz_cuts)[0]; z<(*cosz_cuts)[1]; z++)
		{
			indices[0]=ep;
			indices[1]=z;
  
			for (unsigned int i=0; i<4; i++){   
            	scalar_exp[i]=expectation[i].Index(indices);
			}

			scalar_data=data->Index(indices);

			//raw_scalar_data=data->Index(indices);
			//scalar_data=1.3*raw_scalar_data*pow((*eprox_centers)[ep]/34592.0,0.2); // fixme what is 34592 number?
			//tot_scalar_exp=norm*scalar_exp*pow((*eprox_centers)[ep]/34592.0,gamma); // fixme what is 34592 number?

			tot_scalar_exp = norm*pow(center/E0,-gamma)*(scalar_exp[0]+r_kpi*scalar_exp[2]+r_nubarnu*(scalar_exp[1]+r_kpi*scalar_exp[3]));
	
			prob = LogPoissonProbability(scalar_data,tot_scalar_exp);
			if (std::isnan(prob))
			{
				prob=0;
			}
			llh+=prob;
			count++;
		}
	}	


	double norm_penalty = LogGaussianProbability(norm,norm_mean,norm_sigma);
	double gamma_penalty = LogGaussianProbability(gamma,gamma_mean,gamma_sigma);
	double r_kpi_penalty = LogGaussianProbability(r_kpi,r_kpi_mean,r_kpi_sigma);
	double r_nubarnu_penalty = LogGaussianProbability(r_nubarnu,r_nubarnu_mean,r_nubarnu_sigma);
	
//	std::cout << "NORM PENALTY  " << norm_penalty << std::endl;
//	std::cout << "GAMMA PENALTY  " << gamma_penalty << std::endl;
//	std::cout << "GAMMA  " << gamma << std::endl;
//	std::cout << "GAMMAmu  " << gamma_mean << std::endl;
//	std::cout << "GAMMAsig  " << gamma_sigma << std::endl;


	llh+=norm_penalty + gamma_penalty + r_kpi_penalty + r_nubarnu_penalty;


	//Calculate saturated log-likelihood

	for(unsigned int ep=(*eprox_cuts)[0]; ep<(*eprox_cuts)[1]; ep++)
   	{
		for(unsigned int z=(*cosz_cuts)[0]; z<(*cosz_cuts)[1]; z++)
		{
			indices[0]=ep;
			indices[1]=z;
			scalar_data=data->Index(indices);
			//raw_scalar_data=data->Index(indices);
			//scalar_data=1.3*raw_scalar_data*pow((*eprox_centers)[ep]/34592.0,0.2); // fixme what is 34592 number?
			double satprob;
			satprob=LogPoissonProbability(scalar_data,scalar_data);
			if (std::isnan(satprob))
			{
				satprob=0;
			}
			sllh+=satprob;

		}
	}	
	

//	std::cout << "CHI2: " << 2*(sllh-llh) << std::endl;
//	std::cout << std::endl;
				
	return 2*(sllh-llh);
}









dlib::matrix<double,0,1> Verosimilitud::Chi2Gradient(const dlib::matrix<double,0,1>& nuisance)
{
	unsigned int indices[2];
	double scalar_exp[4];
	double scalar_data;
	double tot_scalar_exp;

	const double norm = nuisance(0);
	const double gamma = nuisance(1);
	const double r_kpi = nuisance(2);		//ratio of kaons to pions
	const double r_nubarnu = nuisance(3);	//ratio of nubars to nus	
	
	const double E0 = 34592.0;
	
	double raw_scalar_data;

	double center;
	
	double grad0=0;
	double grad1=0;
	double grad2=0;
	double grad3=0;


	double strange_constant=34592.0;

	for(unsigned int ep=(*eprox_cuts)[0]; ep<(*eprox_cuts)[1]; ep++)
   	{
		center=(*eprox_centers)[ep];
		for(unsigned int z=(*cosz_cuts)[0]; z<(*cosz_cuts)[1]; z++)
		{
			indices[0]=ep;
			indices[1]=z;
			
			for (unsigned int i = 0; i<4; i++){
            	scalar_exp[i]=expectation[i].Index(indices);
			}
			
			scalar_data=data->Index(indices);

			//raw_scalar_data=data->Index(indices);
			//scalar_data=1.3*raw_scalar_data*pow((*eprox_centers)[ep]/34592.0,0.2); // fixme what is 34592 number?

			tot_scalar_exp=norm*pow(center/E0,-gamma) * ((scalar_exp[0]+r_kpi*scalar_exp[2])+r_nubarnu*(scalar_exp[1]+r_kpi*scalar_exp[3]));
			
			if(!(std::isnan(LogPoissonProbability(scalar_data,tot_scalar_exp))))
			{
				grad0 -= tot_scalar_exp/norm * (scalar_data/tot_scalar_exp - 1);
				grad1 -= -tot_scalar_exp*log(center/E0) * (scalar_data/tot_scalar_exp - 1);
				grad2 -= norm*pow(center/E0,-gamma)*(scalar_exp[2] + r_nubarnu*scalar_exp[3]) * (scalar_data/tot_scalar_exp - 1);
				grad3 -= norm*pow(center/E0,-gamma)*(scalar_exp[1] + r_kpi*scalar_exp[3]) * (scalar_data/tot_scalar_exp - 1);
			}

		}
	}	

	grad0+=(norm-norm_mean)/pow(norm_sigma,2);
	grad1+=(gamma-gamma_mean)/pow(gamma_sigma,2);
	grad2+=(r_kpi-r_kpi_mean)/pow(r_kpi_sigma,2);
	grad3+=(r_nubarnu-r_nubarnu_mean)/pow(r_nubarnu_sigma,2);
	
	grad0*=2;
	grad1*=2;
	grad2*=2;
	grad3*=2;


	dlib::matrix<double,0,1> gradient(4);

	gradient(0) = grad0;
	gradient(1) = grad1;
	gradient(2) = grad2;
	gradient(3) = grad3;

	return gradient;
}

