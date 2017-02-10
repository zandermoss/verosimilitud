#include "Verosimilitud.h"
#include <iostream>
#include <stdexcept>
//#include <Python.h>
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

void Verosimilitud::init(unsigned int numneu,
                         char* char_data_path, char* char_flux_path, char* char_effective_area_path)
{
  // Convert char to string
  std::string data_path(char_data_path);
  std::string flux_path(char_flux_path);
  std::string effective_area_path(char_effective_area_path);

  // initialize cache
  gradient_cache = dlib::matrix<double, 0, 1>(5);

  simps_nintervals = 2;

  eff_area = new EffectiveArea(effective_area_path);
  conv_flux = new ConventionalFlux(flux_path);

	num_efficiencies = eff_area->GetNEff(); 
  NeutrinoEnergyEdges = eff_area->GetEdge(0);
  CosZenithEdges = eff_area->GetEdge(1);
  EnergyProxyEdges = eff_area->GetEdge(2);

	std::string data_fname = data_path;

  unsigned int datadims[2];
  datadims[0] = EnergyProxyBins;
  datadims[1] = CosZenithBins;
  data = new Tensor(2, datadims, 0);
  icd = new ICData(EnergyProxyEdges, CosZenithEdges);

  icd->OpenCSV(data_fname.c_str());
  icd->ReadCSV();
  icd->BinData(data);

  // Apply the appropriate cuts
  eprox_cuts = new std::vector<double>(2, 0);
  cosz_cuts = new std::vector<double>(2, 0);

  // eprox bin	#		eprox center
  //			5		357.167
  //			6		449.647
  //			22		17900.8
  //			23		22535.7

  //(*eprox_cuts)[0] = 6;
  //(*eprox_cuts)[1] = 23;
  (*eprox_cuts)[0]=0;
  (*eprox_cuts)[1]=EnergyProxyBins;

  (*cosz_cuts)[0] = 0;
  (*cosz_cuts)[1] = CosZenithBins;

  // calculate bin centers

  eprox_centers = new std::vector<double>(EnergyProxyBins, 0);
  coszenith_centers = new std::vector<double>(CosZenithBins, 0);
  unsigned int ip1;

  for (unsigned int i = 0; i < EnergyProxyBins; i++) {
    ip1 = i + 1;
    (*eprox_centers)[i] =
        (EnergyProxyEdges->Index(&i) + EnergyProxyEdges->Index(&ip1)) / 2.0;
  }

  for (unsigned int i = 0; i < CosZenithBins; i++) {
    ip1 = i + 1;
    (*coszenith_centers)[i] =
        (CosZenithEdges->Index(&i) + CosZenithEdges->Index(&ip1)) / 2.0;
  }

  double *eprox_edges_array = EnergyProxyEdges->GetDataPointer();
  unsigned int eprox_size = EnergyProxyEdges->GetDataLength();
  double *energy_edges_array = NeutrinoEnergyEdges->GetDataPointer();
  unsigned int energy_size = NeutrinoEnergyEdges->GetDataLength();
  double *coszenith_edges_array = CosZenithEdges->GetDataPointer();
  unsigned int coszenith_size = CosZenithEdges->GetDataLength();

  eprox_edges = new std::vector<double>(eprox_edges_array,
                                        eprox_edges_array + eprox_size);
  energy_edges = new std::vector<double>(energy_edges_array,
                                         energy_edges_array + energy_size);
  coszenith_edges = new std::vector<double>(
      coszenith_edges_array, coszenith_edges_array + coszenith_size);


  SetSimpsNIntervals(2);
  CalculateExpectation();
}

Verosimilitud::~Verosimilitud() {
  delete eff_area;
  delete conv_flux;
  delete data;
  delete icd;
  delete eprox_cuts;
  delete cosz_cuts;
  delete eprox_centers;
  delete coszenith_centers;
}

void Verosimilitud::SetSimpsNIntervals(int nintervals) {
  simps_nintervals = nintervals;
}

void Verosimilitud::SetCoszCuts(std::vector<double> cuts) {
  for (unsigned int i = 0; i < cuts.size(); i++) {
    (*cosz_cuts)[i] = cuts[i];
    std::cout << "COSZ CUTS: " << (*cosz_cuts)[i] << std::endl;
  }
}

void Verosimilitud::SetEproxCuts(std::vector<double> cuts) {

  for (unsigned int i = 0; i < cuts.size(); i++) {
    (*eprox_cuts)[i] = cuts[i];
    std::cout << "EPROX CUTS: " << (*eprox_cuts)[i] << std::endl;
  }
}

double Verosimilitud::OscillationProbability(size_t energy_index, size_t zenith_index, size_t anti) const {
  if(anti == 0){
    return nu_osc_prob_array[(*energy_edges).size()*(zenith_index)+(energy_index)];
  } else {
    return nubar_osc_prob_array[(*energy_edges).size()*(zenith_index)+(energy_index)];
  }
}


double Verosimilitud::LinInter(double x,double xM, double xP, double yM, double yP) const {
  double f2=(x-xM)/(xP-xM);
  double f1=1-f2;
  return f1*yM + f2*yP;
}

double Verosimilitud::OscillationProbability(double energy, double costh,
                                             double anti) const {
//  double costh=cos(zenith);
  auto cthit=std::lower_bound((*coszenith_edges).begin(),(*coszenith_edges).end(),costh);
  if(cthit==(*coszenith_edges).end())
    throw std::runtime_error("zenith not found in the array.");
  if(cthit!=(*coszenith_edges).begin())
    cthit--;
  size_t cth_M=std::distance((*coszenith_edges).begin(),cthit);

  auto eit=std::lower_bound((*energy_edges).begin(),(*energy_edges).end(),energy);
  if(eit==(*energy_edges).end())
    throw std::runtime_error("energy not found in the array.");
  if(eit!=(*energy_edges).begin())
    eit--;
  size_t e_M=std::distance((*energy_edges).begin(),eit);

  double phiMM,phiMP,phiPM,phiPP;
  size_t e_size = (*energy_edges).size();
  if(anti < 0.5){
    phiMM=nu_osc_prob_array[e_size*(cth_M)+(e_M)];
    phiMP=nu_osc_prob_array[e_size*(cth_M)+(e_M+1)];
    phiPM=nu_osc_prob_array[e_size*(cth_M+1)+(e_M)];
    phiPP=nu_osc_prob_array[e_size*(cth_M+1)+(e_M+1)];
  } else {
    phiMM=nubar_osc_prob_array[e_size*(cth_M)+(e_M)];
    phiMP=nubar_osc_prob_array[e_size*(cth_M)+(e_M+1)];
    phiPM=nubar_osc_prob_array[e_size*(cth_M+1)+(e_M)];
    phiPP=nubar_osc_prob_array[e_size*(cth_M+1)+(e_M+1)];
  }

  return(LinInter(costh,(*coszenith_edges)[cth_M],(*coszenith_edges)[cth_M+1],
          LinInter(energy,(*energy_edges)[e_M],(*energy_edges)[e_M+1],phiMM,phiMP),
          LinInter(energy,(*energy_edges)[e_M],(*energy_edges)[e_M+1],phiPM,phiPP)));
}

void Verosimilitud::SetDeSolver(pyoscfunc my_de_solver, void *my_user_data) {
  de_solver = my_de_solver;
  user_data = my_user_data;
}

std::vector<unsigned int> Verosimilitud::GetExpDims() { return exp_dimvec; }

std::vector<unsigned int> Verosimilitud::GetDataDims() { return dat_dimvec; }

std::vector<double> Verosimilitud::GetEnergyEdges() { return *energy_edges; }

std::vector<double> Verosimilitud::GetEproxEdges() { return *eprox_edges; }

std::vector<double> Verosimilitud::GetCosZenithEdges() {
  return *coszenith_edges;
}

double Verosimilitud::TestFunction(double e, double z) {
  return log(pow(3 * e, 4) * pow(z, 4));
}

double Verosimilitud::SimpsAvg(double coszmin, double coszmax, double emin,
                               double emax, double anti, int nintervals_e) {

  double width = (emax - emin) / (double)nintervals_e;
  double integral = 0;
  double mean = 0;
  double coszval;

  for (int j = 0; j < 2; j++) {
    coszval = coszmin + (double)j * (coszmax - coszmin);
    //integral += OscillationProbability(emin, acos(coszval), anti);
    integral += OscillationProbability(emin, coszval, anti);
    for (int i = 1; i < (nintervals_e / 2); i++) {
      integral += 2 * OscillationProbability(emin + 2 * (double)i * width,
                                             coszval, anti);
                                            // acos(coszval), anti);
    }
    for (int i = 1; i < (nintervals_e / 2 + 1); i++) {
      integral += 4 * OscillationProbability(emin + (2 * (double)i - 1) * width,
                                             coszval, anti);
                                             //acos(coszval), anti);
    }
    //integral += OscillationProbability(emax, acos(coszval), anti);
    integral += OscillationProbability(emax, coszval, anti);
    integral *= width / 3.0;
    mean += integral / (2 * (emax - emin));
    integral = 0;
  }
  return mean;
}

void Verosimilitud::CalculateExpectation() {
  unsigned int eff_area_indices[3];
	unsigned int flux_indices[2];
  unsigned int which_flux;
  unsigned int exp_dims[2] = {EnergyProxyBins, CosZenithBins};

	//Resize for efficiency entries
	expectation.resize(num_efficiencies);

  for (size_t effi = 0; effi < num_efficiencies; effi++) {
		//Resize for flux entries
  	(expectation[effi]).resize(4);
  	for (size_t fluxi = 0; fluxi < 4; fluxi++) {
    	(expectation[effi])[fluxi] = std::make_shared<Tensor>(2, exp_dims, 0);
		}
  }

  unsigned int counter = 0;

	double livetime = eff_area->GetLivetime();

  // Loop over true energies.
  for (unsigned int e = 0; e < NeutrinoEnergyBins; e++) {
    eff_area_indices[2] = e;
		flux_indices[1]=e;
    unsigned int ep1 = e + 1;
    double eMin = NeutrinoEnergyEdges->Index(&e);
    double eMax = NeutrinoEnergyEdges->Index(&ep1);

    // Loop over matter/antimatter
    for (unsigned int anti = 0; anti < 2; anti++) {
      // Loop over coszenith bins.
      for (unsigned int z = 0; z < CosZenithBins; z++) {
        eff_area_indices[1] = z;
				flux_indices[0] = z;
        unsigned int zp1 = z + 1;
        double CosZenithMin = CosZenithEdges->Index(&z);
        double CosZenithMax = CosZenithEdges->Index(&zp1);

        double oscprob;
        if (ioscillation)
          oscprob = SimpsAvg(CosZenithMin, CosZenithMax, eMin, eMax,
                             (double)anti, simps_nintervals);
        else
          oscprob = 1;

				// Loop over efficiencies
				for (unsigned int effi = 0; effi < num_efficiencies; effi++){
          unsigned int area_indices[2] = {effi, anti};
          area = eff_area->GetArea(area_indices);


	        // Loop over mesons
	        for (unsigned int meson = 0; meson < 2; meson++) {
	
	          if (meson == 0 && anti == 0) {
	            which_flux = 0;
	          } else if (meson == 0 && anti == 1) {
	            which_flux = 1;
	          } else if (meson == 1 && anti == 0) {
	            which_flux = 2;
	          } else if (meson == 1 && anti == 1) {
	            which_flux = 3;
	          }
	
	          flux = conv_flux->GetFlux(&which_flux);

            // Loop over energy proxy bins.
            for (unsigned int ep = 0; ep < EnergyProxyBins; ep++) {
              eff_area_indices[0] = ep;

              counter++;

              double AvgFlux = flux->Index(flux_indices);
              double EffAreaVal = area->Index(eff_area_indices);



              unsigned int exp_indices[2] = {ep, z};

              double last = (expectation[effi])[which_flux]->Index(exp_indices);
//              double current =
//                  last +
//                 	AvgFlux * EffAreaVal * livetime * oscprob;
              double current =
                  last +
                 	AvgFlux * EffAreaVal * oscprob;
              (expectation[effi])[which_flux]->SetIndex(exp_indices, current);
						}
          }
        }
      }
    }
  }
}

std::vector<double> Verosimilitud::GetAreaVec(void) {
  unsigned int area_indices[2] = {0, 0};
  Tensor *my_area = eff_area->GetArea(area_indices);
  double *areaarray = my_area->GetDataPointer();
  unsigned int arealen = my_area->GetDataLength();
  std::cout << "AREALEN: " << arealen << std::endl;
  std::vector<double> areavec(arealen, 0);
  for (unsigned int i = 0; i < arealen; i++) {
    areavec[i] = areaarray[i];
  }
  return areavec;
}



std::vector<double>
Verosimilitud::GetPertExpectationVec(std::vector<double> _nuisance) {
	unsigned int nnuisance = _nuisance.size();
  dlib::matrix<double, 0, 1> nuisance(nnuisance);
  for (unsigned int i = 0; i < nnuisance; i++) {
			nuisance(i) = _nuisance[i];
  }

  unsigned int exp_dims[2] = {EnergyProxyBins, CosZenithBins};
	Tensor* perturbed_expectation = new Tensor(2,exp_dims,0);
	PerturbExpectation(nuisance, perturbed_expectation);
  double *exparray = perturbed_expectation->GetDataPointer();
  unsigned int explen = perturbed_expectation->GetDataLength();
  std::vector<double> expvec(exparray, exparray + explen);
  unsigned int expdims[2];
  perturbed_expectation->GetDims(expdims);
  exp_dimvec.push_back(expdims[0]);
  exp_dimvec.push_back(expdims[1]);
	delete perturbed_expectation;
  return expvec;
}


std::vector<double> Verosimilitud::GetDataVec(void) {
  double *datarray = data->GetDataPointer();
  unsigned int datlen = data->GetDataLength();
  std::vector<double> datvec(datarray, datarray + datlen);
  unsigned int datdims[] = {0, 0};
  data->GetDims(datdims);
	std::cout << "DATDIM0: " << datdims[0] << std::endl;
	std::cout << "DATDIM1: " << datdims[1] << std::endl;
  dat_dimvec.push_back(datdims[0]);
  dat_dimvec.push_back(datdims[1]);
  return datvec;
}


std::vector<double> Verosimilitud::MinLLH(std::vector<double> param,
                                          std::vector<double> low_bound,
                                          std::vector<double> high_bound,
                                          std::vector<bool> param_to_minimize) {
  if (param.size() != 5) {
    std::cout << "param size not equal to 5. break." << std::endl;
    exit(1);
  }
  if (param.size() != param_to_minimize.size()) {
    std::cout << "sizes param do not match. break" << std::endl;
    exit(1);
  }

  if (low_bound.size() != param_to_minimize.size()) {
    std::cout << "sizes low do not match. break" << std::endl;
    exit(1);
  }

  if (high_bound.size() != param_to_minimize.size()) {
    std::cout << "sizes high do not match. break" << std::endl;
    exit(1);
  }


  unsigned int number_of_parameters_to_minimize =
      std::count(param_to_minimize.begin(), param_to_minimize.end(), true);

  // set initial values and boundaries
  dlib::matrix<double, 0, 1> nuisance(number_of_parameters_to_minimize);
  dlib::matrix<double, 0, 1> lo_bounds(number_of_parameters_to_minimize);
  dlib::matrix<double, 0, 1> hi_bounds(number_of_parameters_to_minimize);
  unsigned int j = 0;
  for (unsigned int i = 0; i < param.size(); i++) {
    if (param_to_minimize[i]) {
      nuisance(j) = param[i];
      lo_bounds(j) = low_bound[i];
      hi_bounds(j) = high_bound[i];
      j++;
    }
  }

  dlib::find_min_box_constrained(dlib::bfgs_search_strategy(),
                 dlib::objective_delta_stop_strategy(1e-7),
                 LLH_caller(this,param,param_to_minimize),
                 LLHgrad_caller(this,param,param_to_minimize),
                 nuisance,
                 lo_bounds,
                 hi_bounds);

//  dlib::find_min_box_constrained(
//      dlib::bfgs_search_strategy(), dlib::gradient_norm_stop_strategy(1e-6),
//      //dlib::lbfgs_search_strategy(10), dlib::gradient_norm_stop_strategy(1e-1),
//      LLH_caller(this, param, param_to_minimize),
//      LLHgrad_caller(this, param, param_to_minimize), nuisance, lo_bounds,
//      hi_bounds);

  std::vector<double> ret(param.size() + 1, 0);

  dlib::matrix<double, 0, 1> param_eval(param.size());
  unsigned int jj = 0;
  for (unsigned int i = 0; i < param.size(); i++) {
    if (param_to_minimize[i]) {
      ret[i] = nuisance(jj);
      param_eval(i) = nuisance(jj);
      jj++;
    } else {
      ret[i] = param[i];
      param_eval(i) = param[i];
    }
  }
	//FIXME: return LLH or Chi2 estimate?
  ret[param.size()] = LLH(param_eval);

  return ret;
}


unsigned int Verosimilitud::EfficiencySwitch(double eff) const{
	//Check that the efficiency is in range.
	if (eff<(eff_area->GetEff(0))){
    throw std::runtime_error("Efficiency out of range (low).");
	}
	else if (eff>(eff_area->GetEff(num_efficiencies-1))){ 
    throw std::runtime_error("Efficiency out of range (high).");
	}
	double t_index;
	//Find nearest efficiency tensor index.
	//Not very efficient, but the vector is only five elements long!
	for (unsigned int i=0; i<num_efficiencies-1; i++){
		double effedges[] = {eff_area->GetEff(i),eff_area->GetEff(i+1)};
		if ((eff>=effedges[0]) && (eff<effedges[1])){
			t_index=i;
		}
	}
	return t_index;
}

void Verosimilitud::PerturbExpectation(const dlib::matrix<double, 0, 1> &nuisance, Tensor* perturbed_expectation) const{
  if(nuisance.size()!=5)
    throw std::runtime_error("Incorrect number of nuisance parameters supplied.");

  const double norm = nuisance(0);
  const double gamma = nuisance(1);
  const double r_kpi = nuisance(2);     // ratio of kaons to pions
  const double r_nubarnu = nuisance(3); // ratio of nubars to nus
  const double efficiency = nuisance(4);

	unsigned int low_eff_index = EfficiencySwitch(efficiency);
//	std::cout << "LE Index: " << low_eff_index << std::endl;
	double low_eff = eff_area->GetEff(low_eff_index);
	double high_eff = eff_area->GetEff(low_eff_index+1);
	double eff_pert_scalars[2];


  double center;
  double E0 = 34592.0;

  unsigned int indices[2];
  double scalar_exp[4];

  for (unsigned int ep = (*eprox_cuts)[0]; ep < (*eprox_cuts)[1]; ep++) {
    center = (*eprox_centers)[ep];
    for (unsigned int z = (*cosz_cuts)[0]; z < (*cosz_cuts)[1]; z++) {
      indices[0] = ep;
      indices[1] = z;

			for (unsigned int effi=0; effi<2; effi++){
	      for (unsigned int flux = 0; flux < 4; flux++) {
	        scalar_exp[flux] = (expectation[effi+low_eff_index])[flux]->Index(indices);
	      }
	
	        eff_pert_scalars[effi]= (scalar_exp[0] + r_kpi * scalar_exp[2] +
	                        r_nubarnu * (scalar_exp[1] + r_kpi * scalar_exp[3]));
			}

	   	double fully_perturbed_scalar = norm * pow(center / E0, -gamma) * (((eff_pert_scalars[1] - eff_pert_scalars[0])/(high_eff-low_eff))*(efficiency - low_eff) + eff_pert_scalars[0]);
			perturbed_expectation->SetIndex(indices,fully_perturbed_scalar);
    }
  }
}


double Verosimilitud::LLH(const dlib::matrix<double, 0, 1> &nuisance) const {

  const double norm = nuisance(0);
  const double gamma = nuisance(1);
  const double r_kpi = nuisance(2);     // ratio of kaons to pions
  const double r_nubarnu = nuisance(3); // ratio of nubars to nus
  const double efficiency = nuisance(4); 

  unsigned int indices[2];
  double llh = 0;
  double prob;

	//Calculate the reduced perturbed expectation tensor.
  unsigned int exp_dims[2] = {EnergyProxyBins, CosZenithBins};
	Tensor* perturbed_expectation = new Tensor(2,exp_dims,0);
	PerturbExpectation(nuisance, perturbed_expectation);

  // calculate unsaturated log-likelihood with nuisance-perturbed expectation
  for (unsigned int ep = (*eprox_cuts)[0]; ep < (*eprox_cuts)[1]; ep++) {
    for (unsigned int z = (*cosz_cuts)[0]; z < (*cosz_cuts)[1]; z++) {
      indices[0] = ep;
      indices[1] = z;

      double scalar_data = data->Index(indices);
      double scalar_exp = perturbed_expectation->Index(indices); 

      prob = LogPoissonProbability(scalar_data, scalar_exp);
      if (std::isnan(prob)) {
        prob = 0;
      }
      llh += prob;
    }
  }

	delete perturbed_expectation;


  double norm_penalty = LogGaussianProbability(norm, norm_mean, norm_sigma);
  double gamma_penalty = LogGaussianProbability(gamma, gamma_mean, gamma_sigma);
  double r_kpi_penalty = LogGaussianProbability(r_kpi, r_kpi_mean, r_kpi_sigma);
  double r_nubarnu_penalty =
      LogGaussianProbability(r_nubarnu, r_nubarnu_mean, r_nubarnu_sigma);
	double efficiency_penalty = LogGaussianProbability(efficiency, efficiency_mean, efficiency_sigma);

  llh += norm_penalty + gamma_penalty + r_kpi_penalty + r_nubarnu_penalty + efficiency_penalty;

  return llh;
}

dlib::matrix<double, 0, 1> Verosimilitud::LLHGradient(const dlib::matrix<double, 0, 1> &nuisance) const {
  unsigned int indices[2];
  double scalar_exp[4];

  const double norm = nuisance(0);
  const double gamma = nuisance(1);
  const double r_kpi = nuisance(2);     // ratio of kaons to pions
  const double r_nubarnu = nuisance(3); // ratio of nubars to nus
  const double efficiency = nuisance(4); 

	unsigned int low_eff_index = EfficiencySwitch(efficiency);
	double low_eff = eff_area->GetEff(low_eff_index);
	double high_eff = eff_area->GetEff(low_eff_index+1);
	double eff_pert_scalars[2];
	double eff_scalars[2];

  double center;
  double E0 = 34592.0;

  double grad0 = 0;
  double grad1 = 0;
  double grad2 = 0;
  double grad3 = 0;
  double grad4 = 0;

  for (unsigned int ep = (*eprox_cuts)[0]; ep < (*eprox_cuts)[1]; ep++) {
    center = (*eprox_centers)[ep];
    for (unsigned int z = (*cosz_cuts)[0]; z < (*cosz_cuts)[1]; z++) {
      indices[0] = ep;
      indices[1] = z;

      double scalar_data = data->Index(indices);

			//Taking efficiency gradient.
			for (unsigned int effi=0; effi<2; effi++){
	      for (unsigned int flux = 0; flux < 4; flux++) {
	        scalar_exp[flux] = (expectation[effi+low_eff_index])[flux]->Index(indices);
	      }
	        eff_pert_scalars[effi]= (scalar_exp[0] + r_kpi * scalar_exp[2] +
	                        r_nubarnu * (scalar_exp[1] + r_kpi * scalar_exp[3]));
			}


	   	double fully_perturbed_scalar = norm * pow(center / E0, -gamma) * (((eff_pert_scalars[1] - eff_pert_scalars[0])/(high_eff-low_eff))*(efficiency - low_eff) + eff_pert_scalars[0]);


      if (!(std::isnan(LogPoissonProbability(scalar_data, fully_perturbed_scalar)))) {
	   		grad4 += norm * pow(center / E0, -gamma) * ((eff_pert_scalars[1] - eff_pert_scalars[0])/(high_eff-low_eff))*(scalar_data / fully_perturbed_scalar - 1);
			}


			//Taking first 4 gradients.
	    for (unsigned int flux = 0; flux < 4; flux++) {
				for (unsigned int effi=0; effi<2; effi++){
					eff_scalars[effi] = (expectation[effi+low_eff_index])[flux]->Index(indices);
	      }

	        scalar_exp[flux] =  (((eff_scalars[1] - eff_scalars[0])/(high_eff-low_eff))*(efficiency - low_eff) + eff_scalars[0]);
			}


      if (!(std::isnan(LogPoissonProbability(scalar_data, fully_perturbed_scalar)))) {
        grad0 += (fully_perturbed_scalar / norm) * (scalar_data / fully_perturbed_scalar - 1);
        grad1 += (-1.0)* fully_perturbed_scalar * log(center / E0) *
                 (scalar_data / fully_perturbed_scalar - 1);
        grad2 += norm * pow(center / E0, -gamma) *
                 (scalar_exp[2] + r_nubarnu * scalar_exp[3]) *
                 (scalar_data / fully_perturbed_scalar - 1);
        grad3 += norm * pow(center / E0, -gamma) *
                 (scalar_exp[1] + r_kpi * scalar_exp[3]) *
                 (scalar_data / fully_perturbed_scalar - 1);
			}


		}
	}

	grad0 += (-1.0)*(norm - norm_mean) / pow(norm_sigma, 2);
	grad1 += (-1.0)*(gamma - gamma_mean) / pow(gamma_sigma, 2);
	grad2 += (-1.0)*(r_kpi - r_kpi_mean) / pow(r_kpi_sigma, 2);
	grad3 += (-1.0)*(r_nubarnu - r_nubarnu_mean) / pow(r_nubarnu_sigma, 2);
	grad4 += (-1.0)*(efficiency - efficiency_mean) / pow(efficiency_sigma, 2);
	
  gradient_cache(0) = grad0;
  gradient_cache(1) = grad1;
  gradient_cache(2) = grad2;
  gradient_cache(3) = grad3;
  gradient_cache(4) = grad4;

  return gradient_cache;
}


double Verosimilitud::Chi2(const dlib::matrix<double, 0, 1> &nuisance) const {
  // Calculate saturated log-likelihood
  unsigned int indices[2];
  double sllh = 0;
  for (unsigned int ep = (*eprox_cuts)[0]; ep < (*eprox_cuts)[1]; ep++) {
    for (unsigned int z = (*cosz_cuts)[0]; z < (*cosz_cuts)[1]; z++) {
      indices[0] = ep;
      indices[1] = z;
      double scalar_data = data->Index(indices);
      double satprob = LogPoissonProbability(scalar_data, scalar_data);
      if (std::isnan(satprob)) {
        satprob = 0;
      }
      sllh += satprob;
    }
  }
	double llh = LLH(nuisance);
  return 2. * (sllh - llh);
}

dlib::matrix<double, 0, 1> Verosimilitud::Chi2Gradient(const dlib::matrix<double, 0, 1> &nuisance) const {
	return -2.0*LLHGradient(nuisance);
}



double Verosimilitud::GetLLH(std::vector<double> param) {
  dlib::matrix<double, 0, 1> nuisance(5);
	for (int i=0; i<5; i++){
  	nuisance(i) = (param)[i];
	}

	return LLH(nuisance);
}

double Verosimilitud::GetLLHGradient(std::vector<double> param, unsigned int index) {
  dlib::matrix<double, 0, 1> nuisance(5);
	for (int i=0; i<5; i++){
  	nuisance(i) = (param)[i];
	}
  dlib::matrix<double, 0, 1> gradient(5);
	gradient = LLHGradient(nuisance);
	return gradient(index);
}


double Verosimilitud::GetChi2(std::vector<double> param) {
  dlib::matrix<double, 0, 1> nuisance(5);
	for (int i=0; i<5; i++){
  	nuisance(i) = (param)[i];
	}

	return Chi2(nuisance);
}

double Verosimilitud::GetChi2Gradient(std::vector<double> param, unsigned int index) {
  dlib::matrix<double, 0, 1> nuisance(5);
	for (int i=0; i<5; i++){
  	nuisance(i) = (param)[i];
	}
  dlib::matrix<double, 0, 1> gradient(5);
	gradient = Chi2Gradient(nuisance);
	return gradient(index);
}


/*
void Print2Tensor(Tensor& t){
	unsigned int rank = t.GetRank();
	unsigned int* dims = new unsigned int[rank];
	t.GetDims(dims);
	unsigned int index[3];
	for (unsigned int row=0; row<dims[0]; row++){
		index[0] = row;
		for (unsigned int col=0; col<dims[1]; col++){
			index[1] = col;
			for (unsigned int dep=0; dep<dims[1]; dep++){
				index[2]=dep;	
				std::cout << t.Index(index) << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
	}
	return;
}
*/


int main(void){
	
	return 0;
}


