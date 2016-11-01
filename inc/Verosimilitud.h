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

//! The log-likelihood/Chi-squared calculation class.
/*!
	The Verosimilitud class has three essential functions:
	-1
		Calculate a histogram of expected rate over reconstructed energy and angle. This calculation uses survival probabilities from NuSHEEP as well as effective area and flux information from the data release. The specifics of the given oscillation model are stipulated in the calling script using the interface to NuSHEEP.
	-2
		Load IceCube data from the IceCube release.
	-3
		Calculate a Chi^2 corresponding to the prediction of the specified oscillation model and the IC (IceCube) data. Saturated and unsaturated log-likelihood values are calculated, and then used to approximate the corresponding chi^2. This chi^2 is then minimized over the relevant nuisance parameters with penalties distributed according to prior nuisance distributions.

FIXME add links to the relevant functions!
*/

class Verosimilitud {

	public:
		//--------------------------------------------------------//
		//! The Constructor.
		/*!
			The constructor for Verosimilitud initializes data structures and loads in IC data as well as effective detector area and atmospheric muon neutrino fluxes. Data is loaded and stored in an ICData object. Effective area is loaded and stored in an EffectiveArea object, and flux is loaded and stored in a ConventionalFlux object. Default energy and zenith angle cuts are applied (no cut). Finally, bin centers are calculated from "::"NeutrinoEnergyEdges", "::"CosZenithEdges, and "::"EnergyProxyEdges for use in simulation.
			\param my_numneu An integer specifying the number of neutrinos, or, more precisely, the dimension of the neutrino state space.
			\param loyear An integer providing a lower bound for the loops over data years.
			\param hiyear An integer providing an upper bound for the loops over data years. If [loyear,hiyear]=[0,2], both the 2010 and 2011 datasets are included. If [0,1], only 2010, and if [1,2], only 2011.
		*/


		Verosimilitud(unsigned int my_numneu,unsigned int loyear, unsigned int hiyear);


		

		//--------------------------------------------------------//
		//! The Destructor
		/*!
			Just cleaning up dynamic memory. No special operations.
		*/

		~Verosimilitud();


	
		//--------------------------------------------------------//
		//! Typedef for python function for NuSHEEP survival probability calculation.
		/*!
			Some arcane C trickery for passing relevant function data using the void pointer userdata. I never actually came to understand how this works. I just swiped the code from stackexchange:
			\param userdata void pointer used for sordid data transmission.
			\param argument	the actual argument, whose nature will be clarified in the documentation surrounding "::"OscillationProbability
		*/

		typedef double (*pyoscfunc)(std::vector<double> argument, void *userdata);

		//--------------------------------------------------------//
		//! Function wrapping a call to the NuSHEEP survival probability calculation.
		/*!
			OscillationProbability is the C++ wrapper used to call the pyoscfunc function proveded by NuSHEEP.
			\param energy is the neutrino energy, in GeV.
			\param zenith is the zenith angle (angle from the south pole) corresponding to the vector along which the neutrino propagates from its atmospheric production point to the IC detector at the south pole. For the purposes of these simulations, I have omitted zenith angles from 0 to pi/2, which correspond to production in the sky above the IC detector. The neutrinos only intersect the earth at angles between pi/2 and pi. FIXME: shift this discussion to a more appropriate place? NuSHEEP dox maybe?
			
			\param is a switch between neutrino and anti-neutrino mode. 0 indicates a neutrino, and 1 an antineutrino.	
		*/

	    double OscillationProbability(double energy,double zenith, double anti);

		//--------------------------------------------------------//
		//! A function to cast the appropriate NuSHEEP python function as a pyoscfunc. 
		/*! This function is called from a python script. The relevant NuSHEEP function is passed to SetDeSolver() as a python function object. The pyoscfunc pointer is set to point to this function. Now, when OscillationProbability is called from C++, the proper python function is invoked.

		\param de_solve is the python function passed to C++
		\param user_data is the void pointer used (through witchcraft) to pass some sort of instance data between python and C++, and vice versa.
		*/	 
	
		void SetDeSolver(pyoscfunc de_solve, void* user_data);

		//--------------------------------------------------------//
		//! Returns a vector of ranks corresponding to the components of the expectation tensor filled by  "::"CalculateExpectation.
		/*!
			This tensor, of Tensor type, should be considered as a two dimensional histogram in energy proxy and zenith angle. The vector returned then provides the number of bins in each dimension.

			FIXME reference Tensor class!
		*/

		std::vector<unsigned int> GetExpDims(void);

		//--------------------------------------------------------//
		//! Returns a vector of ranks corresponding to the components of the data tensor read into (and binned) by ICData FIXME reference ICData class.
		/*!
			This tensor, of Tensor type, should be considered as a two dimensional histogram in energy proxy and zenith angle. The vector returned then provides the number of bins in each dimension.

			FIXME reference Tensor class!
		*/

		std::vector<unsigned int> GetDataDims(void);

		//--------------------------------------------------------//
		//! Returns a vector of bin edges corresponding to the energy proxy variable.
		/*!
			These bin edges are used to bin the data and construct the expectation. The edges are extracted as a tensor object from the effective_area object
			FIXME reference effective area class and Tensor class.
			This tensor object is then converted to a std::vector object. 
			\return The edges as a vector<double>
		*/

		std::vector<double> GetEproxEdges(void);

		//--------------------------------------------------------//
		//! Returns a vector of bin edges corresponding to the true energy variable.
		/*!
			These bin edges are used to bin the data and construct the expectation. The edges are extracted as a tensor object from the effective_area object
			FIXME reference effective area class and Tensor class.
			This tensor object is then converted to a std::vector object. 
			\return The edges as a vector<double>
		*/

		std::vector<double> GetEnergyEdges(void);

		//--------------------------------------------------------//
		//! Returns a vector of bin edges corresponding to the cosine of the zenith angle.
		/*!
			These bin edges are used to bin the data and construct the expectation. The edges are extracted as a tensor object from the effective_area object
			FIXME reference effective area class and Tensor class.
			This tensor object is then converted to a std::vector object. 
			\return The edges as a vector<double>
		*/

		std::vector<double> GetCosZenithEdges(void);

		//--------------------------------------------------------//
		//! Calculates expected event counts in each bin of the parameter space occupied by the IC data. 
		/*!
		Effective detector area, atmospheric muon neutrino flux, livetime, detector corrections, and "::"OscillationProbability. The calculation is a nest of summations and integrations. Summation is performed over matter and antimatter neutrinos and detector years (selectable in the argument to "::"Verosimilitud). Integration is performed over the true neutrino energy, which is used to calculate oscillation probabilities. expectation, a class member (pointer to a "::"Tensor object) represents a 2D histogram in cos(zenith) and energy proxy. This tensor is initialized and filled, bin-by-bin, in this loop. The oscillation probabilities here are not direct returns from "::"OscillationProbability, but rather averages over the bin obtained using Simpson's rule as implemented in "::"SimpsAvg.
		*/  
		void CalculateExpectation(void);




		//--------------------------------------------------------//
		//! Returns the muon neutrino flux distribution over true energy as a vector.
		/*!
			Extracts flux as an array of doubles from conv_flux, which is a "::"ConventionalFlux object. Converts the array to a vector. For a description of the utility of this function, see the documentation for "::"GetFluxVec
			\return distribution of flux over true energy.
		*/

		std::vector<double> GetFluxVec(void);


		//--------------------------------------------------------//
		//! Returns the linearized effective area tensor as a vector.
		/*!
			Extracts linearized effective area as an array of doubles from eff_area, which is an "::"EffectiveArea object. Converts the array to a vector. For reference, eff_area is defined (before linearization) over true energy, energy proxy, and cos(zenith). In concert with "::"GetEnergyEdges, "::"GetEproxEdges", and "::"GetCosZenithEdges, the original tensor can be reconstructed from the linearized output. Why linearize? If we wish to extract a tensor from a C++ function into a python script through the cython wrappers, it is easiest to extract a vector. Supplied with dimensions and edge vectors, one can loop through the linearized data and reconstruct the tensor. Although clunky, this allows full access to the data structures in Verosimilitud from the python level, which is rather useful for debugging, plotting, etc. 
			\return vector linearization of expectation tensor.
		*/
		std::vector<double> GetAreaVec(void);


		//--------------------------------------------------------//
        //! Returns the linearized expectation tensor as a vector.
        /*!
            Extracts the linearized expectation tensor as an array of doubles from expectation, which is a "::"Tensor object. expectation is calculated by "::"CalculateExpectation. Converts the array to a vector. Before linearization, the expectation tensor elements are indexed by cos(zenith) and energy proxy. For a description of the utility of
 this function, see the documentation for "::"GetFluxVec
            \return vector linearization of expectation tensor.
        */

		std::vector<double> GetExpectationVec(void);

		//--------------------------------------------------------//
        //! Returns the linearized, perturbed expectation tensor as a vector.
        /*!
			Performs the same extraction/linearization as "::"GetExpectationVec, but first applies perturbations according to input nuisance parameters. For details of these nuisance perturbations, see "::"Chi2.
			\param nuisance a vector of nuisance parameters. In its current state, nuisance[0] is a norm which multiplies every element in the expectation tensor. nuisance[1] is the spectral index, which applies a power-function tilt along the energy axis.  

            \return vector linearization of nuisance-perturbed expectation tensor.
        */

		std::vector<double> GetPertExpectationVec(std::vector<double> nuisance);


		//--------------------------------------------------------//
        //! Returns the linearized data tensor as a vector.
        /*!
            Extracts the linearized data tensor as an array of doubles from data, which is a "::"Tensor object. data is filled (in "::"Verosimilitud) by an "::"ICData object, which reads in and bins IceCube data into a 2-d histogram of the same shape as expectation. Converts the array to a vector. Before linearization, the data tensor elements are indexed by cos(zenith) and energy proxy. For a description of the utility of
 this function, see the documentation for "::"GetFluxVec
            \return vector linearization of data tensor.
        */
		std::vector<double> GetDataVec(void);

		

		//--------------------------------------------------------//
        //! Minimizes Chi2 from data/expectation log-likelihood calculation over nuisance parameters.
        /*!
			Perturbation of the expectation tensor according to a set of nuisance parameters is a method of quantifying the effects of uncertainties in physical parameters (e.g. flux normalization and spectral index) on the Chi2 and best fit points/confidence intervals. 
			Each nuisance parameter (a random variable corresponding to an uncertian physical parameter) is drawn from some prior distribution. To calculate a Chi2, the expectation is perturbed according to some vector of nuisance parameters, and the likelihood is penalized according to the pdf of this vector. If the vector elements are in the tails of their (scalar) pdfs, then, although the naive Chi2 may be small, the probability of sampling these values of the physial nuisance parameters from their underlying distributions is small, and the final likelihood (and therefore Chi2) must reflect this (for a detailed description of this calculation, consult the documentation for "::"Chi2).
			The Chi2 is then minimized, using a constrained box minimization routine from DLib (FIXME ref dlib) over the space of nuisance parameters. This minimizer makes calls to "::"Chi2 and "::"Chi2Gradient".

		\param nuisance an initial guess for the minimal-Chi2 nuisance parameters. The current configuration is a two-element vector with normalization in the first slot, and spectral index in the second.
		\return a vector containing the minimal Chi2 value (position 0), and then the 			minimal-chi2 values for each nuisance parameter (indices>0)
            
        */

		std::vector<double> Chi2MinNuisance(std::vector<double> nuisance);


		//--------------------------------------------------------//
        //! Calculates Chi2 from IC data and the calculated expectation tensors. Perturbs the expectation tensor according to input nuisance parameters.
        /*!
			In the current implementation, the total normalization of the expectation tensor is multiplied by nuisance[0]. The spectral index, nuisance[1] is used to apply a power-law tilt to the expectation tensor along the energy axis.
			The funcion then iterates over elements of data and expectation tensors, and calculates log-likelihood assuming events in each bin are poisson disributed. Total log-likelihood is then penalized according to the gaussian pdfs describing both nuisance parameters. The mean and variance of each distribution is stored as a class member variable.
			Chi2 is then approximated from calculated log-likelihood and saturated log-likelihood values.
			\param a vector of nuisance parameters. Currently: [normalization, spectral index]
			\return approximate Chi2
		*/

		double Chi2(const dlib::matrix<double,0,1>& nuisance);

		//--------------------------------------------------------//
        //! Calculates the gradient of the Chi2 returned by "::"Chi2 in the space of nuisance parameters. 
        /*!
			The form of this gradient was calculated by hand analytically in the interest of speed in the minimization, and will have to be altered when further nuisance parameters are added, although numerical differentiation is an option.
			\param the vector of nuisance parameters, in the dlib::matrix format for compatibility with the minimizer.
			\return the gradient in nuisance space (love that name), again in the dlib::matrix format.
		*/

		dlib::matrix<double,0,1> Chi2Gradient(const dlib::matrix<double,0,1>& nuisance);


		//--------------------------------------------------------//
        //! Approximates the average value of OscillationProbability over a bin in energy proxy and cos(zenith).
		/*!
			Applies the simpson method for numerical integration to approximate the average of OscillationProbability. Recall that Avg(f(x),a<=x<=b) = 1/(b-a)*Integral(f(x)*dx,a,b). This is true in more dimensions (integrating over some volume), so numerical approximation of the integral amounts to numerical approximation of the average. Computing the average over an bin in the expecatation calulation yields a more accurate estimate of the expectation. Currently, the approximation is applied in the energy proxy dimension. The default is two intervals in energy, but can be set using "::"SetSimpsNIntervals().
		This function is called from "::"CalculateExpectation(), and wraps the calls to "::"OscillationProbability().
		\param coszmin the lower bound for cosz in the bin.
		\param coszmax the upper bound for cosz in the bin.
		\param emin the lower bound for energy proxy in the bin.
		\param emax the upper bound for energy proxy in the bin.
		\param anti 0 if neutrino mode, 1 if antineutrino mode. Necessary to pass to "::"OscillationProbability()

		\return the approximate average of oscillation probability over the bin.

		*/

		double SimpsAvg(double coszmin, double coszmax, double emin, double emax, double anti, int nintervals);


		//--------------------------------------------------------//
		//! Sets the number of intervals over which to apply the simpson approximation of oscillation probabilities.
		/*!
		Sets the simps_nintervals parameter for use by "::"SimpsAvg.
		\see SimpsAvg()
		*/ 

		void SetSimpsNIntervals(int nintervals);

		//--------------------------------------------------------//
        //! A function used to test the simpson approximation. 
		/*! 
			This function was used to test the accuracy of the simpson averaging approimation corresponding to the calculation of oscillation probability bin-by-bin. I'm leaving it here because any tweaking of the simpson approximation (which might be required for more accuracy), should be tested.
			\see SimpsAvg()

			\param e a dummy representing the energy proxy.
			\param z a dummy representing the cos(zenith). 
			\return a dummy representing the oscillation probability.
		*/

		double TestFunction(double e, double z);





		//--------------------------------------------------------//
		//! Sets cuts on energy proxy range to use in the calculation of Chi2.
		/*!
			Due to poor energy reconstruction at low energies, and astrophysical contamination at high energies, the top and tail of the energy distribution should be cut.
			\param cuts a two-element vector with the lower bound first, and the upper bound second.
		*/

		void SetEproxCuts(std::vector<double> cuts);


		//--------------------------------------------------------//
		//! Sets cuts on the cos(zenith) range to use in the calculation of Chi2.
		/*!
			\param cuts a two-element vector with the lower bound first, and the upper bound second.
		*/

		void SetCoszCuts(std::vector<double> cuts);




		//--------------------------------------------------------//
		//! The number of bins to use in expectation along the true energy axis.

		const unsigned int NeutrinoEnergyBins=280;


		//--------------------------------------------------------//
		//! The number of bins to use in expectation and data along the cos(zenith) axis.

    	const unsigned int CosZenithBins=11;

		//--------------------------------------------------------//
		//! The number of bins to use in expectation and data along the energy proxy axis.

    	const unsigned int EnergyProxyBins=50;

		

	protected:


		unsigned int data_years[2];

		int simps_nintervals;

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
		std::vector<double>* energy_edges;
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
		
		ICData* icd[2];
	
	    ConventionalFlux* conv_flux;
	    Tensor* flux;
	    Tensor* detcorr;
	



		std::vector<double>* eprox_cuts;
		std::vector<double>* cosz_cuts;
	
//Nuisance parameters:

	double norm_mean=1;
	double norm_sigma=0.4;
	double gamma_mean=0;
	double gamma_sigma=0.05;
	double r_kpi_mean=1;
	double r_kpi_sigma=0.1;
	double r_nubarnu_mean=1;
	double r_nubarnu_sigma=0.025;


   
// 	   Set_data()
    
//	   Set_simulation()
};

class Chi2_caller
{
	private:
		Verosimilitud* my_verosim;
	public:
		Chi2_caller(Verosimilitud * verosim)
		{
			my_verosim = verosim;
		}


		double operator() (const dlib::matrix<double,0,1>& nuisance) const
		{
			return my_verosim->Chi2(nuisance);
		}

};




class Chi2grad_caller
{
	private:
		Verosimilitud* my_verosim;
	public:

		Chi2grad_caller(Verosimilitud * verosim)
		{
			my_verosim = verosim;
		}


		dlib::matrix<double,0,1> operator() (const dlib::matrix<double,0,1>& nuisance) const
		{
			return my_verosim->Chi2Gradient(nuisance);
		}

};



#endif // __VEROSIMILILTUD_H_INCLUDED__
