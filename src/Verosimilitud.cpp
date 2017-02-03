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
                         unsigned int loyear,
                         unsigned int hiyear,
                         char* char_data_path, char* char_flux_path, char* char_effective_area_path, char* char_detector_correction_path)
{
  // Convert char to string
  std::string data_path(char_data_path);
  std::string flux_path(char_flux_path);
  std::string effective_area_path(char_effective_area_path);
  std::string detector_correction_path(char_detector_correction_path);


  // initialize cache
  gradient_cache = dlib::matrix<double, 0, 1>(4);

  // Which years do we want to include?
  data_years[0] = loyear;
  data_years[1] = hiyear;

  simps_nintervals = 2;

  conv_flux = new ConventionalFlux(detector_correction_path,flux_path);
  eff_area = new EffectiveArea(effective_area_path);

  unsigned int edge_indices[4] = {0, 0, 0, 0};
  NeutrinoEnergyEdges = eff_area->GetEdge(edge_indices);
  edge_indices[3] = 1;
  CosZenithEdges = eff_area->GetEdge(edge_indices);
  edge_indices[3] = 2;
  EnergyProxyEdges = eff_area->GetEdge(edge_indices);

  // Bin data from 2010, 2011, or both!
  std::string fnames[2] = {data_path+"/2010.dat", data_path+"/2011.dat"};

  unsigned int datadims[2];
  datadims[0] = EnergyProxyBins;
  datadims[1] = CosZenithBins;
  data = new Tensor(2, datadims, 0);

  for (unsigned int y = 0; y < 2; y++) {
    icd[y] = new ICData(EnergyProxyEdges, CosZenithEdges);
  }

  for (unsigned int y = data_years[0]; y < data_years[1]; y++) {
    icd[y]->OpenCSV(fnames[y].c_str());
    icd[y]->ReadCSV();
    icd[y]->BinData(data);
  }

  // Apply the appropriate cuts
  eprox_cuts = new std::vector<double>(2, 0);
  cosz_cuts = new std::vector<double>(2, 0);

  // eprox bin	#		eprox center
  //			5		357.167
  //			6		449.647
  //			22		17900.8
  //			23		22535.7

  (*eprox_cuts)[0] = 6;
  (*eprox_cuts)[1] = 23;
  //	(*eprox_cuts)[0]=0;
  //	(*eprox_cuts)[1]=EnergyProxyBins;

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
  delete icd[0];
  delete icd[1];
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

double Verosimilitud::OscillationProbability(double energy, double zenith,
                                             double anti) const {
  double costh=cos(zenith);
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
    integral += OscillationProbability(emin, acos(coszval), anti);
    for (int i = 1; i < (nintervals_e / 2); i++) {
      integral += 2 * OscillationProbability(emin + 2 * (double)i * width,
                                             acos(coszval), anti);
    }
    for (int i = 1; i < (nintervals_e / 2 + 1); i++) {
      integral += 4 * OscillationProbability(emin + (2 * (double)i - 1) * width,
                                             acos(coszval), anti);
    }
    integral += OscillationProbability(emax, acos(coszval), anti);
    integral *= width / 3.0;
    mean += integral / (2 * (emax - emin));
    integral = 0;
  }
  return mean;
}

void Verosimilitud::CalculateExpectation() {
  unsigned int dyn_indices[3];
  unsigned int which_flux;
  unsigned int exp_dims[2] = {EnergyProxyBins, CosZenithBins};

  expectation.resize(4);
  for (size_t i = 0; i < 4; i++) {
    expectation[i] = std::make_shared<Tensor>(2, exp_dims, 0);
  }

  unsigned int counter = 0;

  // Loop over true energies.
  for (unsigned int e = 0; e < NeutrinoEnergyBins; e++) {
    dyn_indices[0] = e;
    unsigned int ep1 = e + 1;
    double eMin = NeutrinoEnergyEdges->Index(&e);
    double eMax = NeutrinoEnergyEdges->Index(&ep1);

    // Loop over matter/antimatter
    for (unsigned int anti = 0; anti < 2; anti++) {
      // Loop over coszenith bins.
      for (unsigned int z = 0; z < CosZenithBins; z++) {
        dyn_indices[1] = z;
        unsigned int zp1 = z + 1;
        double CosZenithMin = CosZenithEdges->Index(&z);
        double CosZenithMax = CosZenithEdges->Index(&zp1);

        double oscprob;
        if (ioscillation)
          oscprob = SimpsAvg(CosZenithMin, CosZenithMax, eMin, eMax,
                             (double)anti, simps_nintervals);
        else
          oscprob = 1;

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
          // Loop over years 2010, 2011
          for (unsigned int year = data_years[0]; year < data_years[1];
               year++) {
            // Do muon neutrinos only: no loop over flavor.
            unsigned int flavor = 0;

            double livetime = eff_area->GetLivetime(year);

            unsigned int area_indices[3] = {year, flavor, anti};
            area = eff_area->GetArea(area_indices);

            detcorr = conv_flux->GetDetCorr(&year);
            // Loop over energy proxy bins.
            for (unsigned int ep = 0; ep < EnergyProxyBins; ep++) {
              dyn_indices[2] = ep;

              counter++;

              double FluxIntegral = flux->Index(dyn_indices);
              // implement DOM correction?
              double DomCorr = detcorr->Index(dyn_indices);
              double EffAreaVal = area->Index(dyn_indices);
              EffAreaVal *= 1.0e4; // m^2 to cm^2

              unsigned int exp_indices[2] = {ep, z};

              double last = expectation[which_flux]->Index(exp_indices);
              double current =
                  last +
                  FluxIntegral * EffAreaVal * livetime * DomCorr * oscprob;

              expectation[which_flux]->SetIndex(exp_indices, current);
            }
          }
        }
      }
    }
  }
}

std::vector<double> Verosimilitud::GetAreaVec(void) {
  unsigned int area_indices[3] = {0, 0, 0};
  Tensor *my_area = eff_area->GetArea(area_indices);
  double *areaarray = my_area->GetDataPointer();
  unsigned int arealen = my_area->GetDataLength();
  std::cout << "AREALEN: " << arealen << std::endl;
  std::vector<double> areavec(arealen, 0);
  for (int i = 0; i < arealen; i++) {
    areavec[i] = areaarray[i];
  }
  return areavec;
}

std::vector<double>
Verosimilitud::GetPertExpectationVec(std::vector<double> nuisance) {
  unsigned int indices[2];
  double scalar_exp[4];
  double pert_scalar_exp;

  double norm = nuisance[0];
  double gamma = nuisance[1];
  double r_kpi = nuisance[2];
  double r_nubarnu = nuisance[3];

  double E0 = 34592.0;
  double center;

  unsigned int exp_dims[2] = {EnergyProxyBins, CosZenithBins};
  Tensor *perturbed_expectation = new Tensor(2, exp_dims, 0);

  for (unsigned int ep = 0; ep < EnergyProxyBins; ep++) {
    center = (*eprox_centers)[ep];

    for (unsigned int z = 0; z < CosZenithBins; z++) {
      indices[0] = ep;
      indices[1] = z;

      for (unsigned int i = 0; i < 4; i++) {
        scalar_exp[i] = expectation[i]->Index(indices);
      }

      pert_scalar_exp = norm * pow(center / E0, -gamma) *
                        (scalar_exp[0] + r_kpi * scalar_exp[2] +
                         r_nubarnu * (scalar_exp[1] + r_kpi * scalar_exp[3]));
      perturbed_expectation->SetIndex(indices, pert_scalar_exp);
    }
  }

  double *exparray = perturbed_expectation->GetDataPointer();
  unsigned int explen = perturbed_expectation->GetDataLength();
  std::vector<double> expvec(exparray, exparray + explen);
  unsigned int expdims[] = {0, 0};
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
  dat_dimvec.push_back(datdims[0]);
  dat_dimvec.push_back(datdims[1]);
  return datvec;
}

std::vector<double> Verosimilitud::MinLLH(std::vector<double> param,
                                          std::vector<double> low_bound,
                                          std::vector<double> high_bound,
                                          std::vector<bool> param_to_minimize) {
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
                 //dlib::objective_delta_stop_strategy(1e-10),
                 dlib::gradient_norm_stop_strategy(1e-8),
                 Chi2_caller(this,param,param_to_minimize),
                 Chi2grad_caller(this,param,param_to_minimize),
                 nuisance,
                 lo_bounds,
                 hi_bounds);

  /*
  dlib::find_min_box_constrained(
      dlib::bfgs_search_strategy(), dlib::gradient_norm_stop_strategy(1e-6),
      //dlib::lbfgs_search_strategy(10), dlib::gradient_norm_stop_strategy(1e-1),
      Chi2_caller(this, param, param_to_minimize),
      Chi2grad_caller(this, param, param_to_minimize), nuisance, lo_bounds,
      hi_bounds);
  */

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
  ret[param.size()] = Chi2(param_eval);

  return ret;
}

double Verosimilitud::Chi2(const dlib::matrix<double, 0, 1> &nuisance) const {
  unsigned int indices[2];
  double scalar_exp[4];
  double scalar_data;
  double llh = 0;
  double sllh = 0;
  int count = 0;
  double prob;
  double tot_scalar_exp;

  const double norm = nuisance(0);
  const double gamma = nuisance(1);
  const double r_kpi = nuisance(2);     // ratio of kaons to pions
  const double r_nubarnu = nuisance(3); // ratio of nubars to nus

  double center;
  double E0 = 34592.0;

  // calculate unsaturated log-likelihood with nuisance-perturbed expectation

  for (unsigned int ep = (*eprox_cuts)[0]; ep < (*eprox_cuts)[1]; ep++) {

    center = (*eprox_centers)[ep];

    //		std::cout<<"eprox_center: "<<center<<" ep: "<<ep<<std::endl;

    for (unsigned int z = (*cosz_cuts)[0]; z < (*cosz_cuts)[1]; z++) {
      indices[0] = ep;
      indices[1] = z;

      for (unsigned int i = 0; i < 4; i++) {

        scalar_exp[i] = expectation[i]->Index(indices);

        // std::cout<<"fixed!"<<std::endl;
      }

      scalar_data = data->Index(indices);

      tot_scalar_exp = norm * pow(center / E0, -gamma) *
                       (scalar_exp[0] + r_kpi * scalar_exp[2] +
                        r_nubarnu * (scalar_exp[1] + r_kpi * scalar_exp[3]));

      prob = LogPoissonProbability(scalar_data, tot_scalar_exp);
      if (std::isnan(prob)) {
        prob = 0;
      }
      llh += prob;
      count++;
    }
  }

  // std::cout<<"V.Chi2, stop loop"<<std::endl;

  double norm_penalty = LogGaussianProbability(norm, norm_mean, norm_sigma);
  double gamma_penalty = LogGaussianProbability(gamma, gamma_mean, gamma_sigma);
  double r_kpi_penalty = LogGaussianProbability(r_kpi, r_kpi_mean, r_kpi_sigma);
  double r_nubarnu_penalty =
      LogGaussianProbability(r_nubarnu, r_nubarnu_mean, r_nubarnu_sigma);

  //	std::cout << "NORM PENALTY  " << norm_penalty << std::endl;
  //	std::cout << "GAMMA PENALTY  " << gamma_penalty << std::endl;
  //	std::cout << "GAMMA  " << gamma << std::endl;
  //	std::cout << "GAMMAmu  " << gamma_mean << std::endl;
  //	std::cout << "GAMMAsig  " << gamma_sigma << std::endl;

  llh += norm_penalty + gamma_penalty + r_kpi_penalty + r_nubarnu_penalty;

  // Calculate saturated log-likelihood

  for (unsigned int ep = (*eprox_cuts)[0]; ep < (*eprox_cuts)[1]; ep++) {
    for (unsigned int z = (*cosz_cuts)[0]; z < (*cosz_cuts)[1]; z++) {
      indices[0] = ep;
      indices[1] = z;
      scalar_data = data->Index(indices);
      double satprob = LogPoissonProbability(scalar_data, scalar_data);
      if (std::isnan(satprob)) {
        satprob = 0;
      }
      sllh += satprob;
    }
  }

  //	std::cout << "CHI2: " << 2*(sllh-llh) << std::endl;
  //	std::cout << std::endl;

  return 2. * (sllh - llh);
}

dlib::matrix<double, 0, 1>
Verosimilitud::Chi2Gradient(const dlib::matrix<double, 0, 1> &nuisance) const {
  unsigned int indices[2];
  double scalar_exp[4];
  double scalar_data;
  double tot_scalar_exp;

  const double norm = nuisance(0);
  const double gamma = nuisance(1);
  const double r_kpi = nuisance(2);     // ratio of kaons to pions
  const double r_nubarnu = nuisance(3); // ratio of nubars to nus

  const double E0 = 34592.0;

  double center;

  double grad0 = 0;
  double grad1 = 0;
  double grad2 = 0;
  double grad3 = 0;

  //	const double strange_constant=34592.0;

  for (unsigned int ep = (*eprox_cuts)[0]; ep < (*eprox_cuts)[1]; ep++) {
    center = (*eprox_centers)[ep];

    for (unsigned int z = (*cosz_cuts)[0]; z < (*cosz_cuts)[1]; z++) {
      indices[0] = ep;
      indices[1] = z;

      for (unsigned int i = 0; i < 4; i++) {
        scalar_exp[i] = expectation[i]->Index(indices);
      }

      scalar_data = data->Index(indices);

      tot_scalar_exp = norm * pow(center / E0, -gamma) *
                       ((scalar_exp[0] + r_kpi * scalar_exp[2]) +
                        r_nubarnu * (scalar_exp[1] + r_kpi * scalar_exp[3]));

      if (!(std::isnan(LogPoissonProbability(scalar_data, tot_scalar_exp)))) {
        grad0 -= tot_scalar_exp / norm * (scalar_data / tot_scalar_exp - 1);
        grad1 -= -tot_scalar_exp * log(center / E0) *
                 (scalar_data / tot_scalar_exp - 1);
        grad2 -= norm * pow(center / E0, -gamma) *
                 (scalar_exp[2] + r_nubarnu * scalar_exp[3]) *
                 (scalar_data / tot_scalar_exp - 1);
        grad3 -= norm * pow(center / E0, -gamma) *
                 (scalar_exp[1] + r_kpi * scalar_exp[3]) *
                 (scalar_data / tot_scalar_exp - 1);
      }
    }
  }

  grad0 += 2.*(norm - norm_mean) / pow(norm_sigma, 2);
  grad1 += 2.*(gamma - gamma_mean) / pow(gamma_sigma, 2);
  grad2 += 2.*(r_kpi - r_kpi_mean) / pow(r_kpi_sigma, 2);
  grad3 += 2.*(r_nubarnu - r_nubarnu_mean) / pow(r_nubarnu_sigma, 2);

  gradient_cache(0) = grad0;
  gradient_cache(1) = grad1;
  gradient_cache(2) = grad2;
  gradient_cache(3) = grad3;

  return 2.*gradient_cache;
}

double Verosimilitud::LLH(std::vector<double> param) {

  dlib::matrix<double, 0, 1> nuisance(4);

  nuisance(0) = (param)[0];
  nuisance(1) = (param)[1];
  nuisance(2) = (param)[2];
  nuisance(3) = (param)[3];

  return Chi2(nuisance);
}
