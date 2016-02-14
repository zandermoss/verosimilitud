#include "Flux.h"

void Flux::Flux(double arg_normalization, double arg_index)
{
	normalization=arg_normalization;
	index=arg_index;
}

double Flux::EvaluatePowerlawFlux(const powerlawFlux flux, double neutrinoEnergy)
{
    //Note that we include here the correction due to (mis-simulating) the detector
    //optical efficiency. It isn't really part of the flux, but it depends on the
    //spectrum (in this case the spectral index), so this is convenient.
    double DOMEffCorr=pow(1.189944,-2.-flux.index);
    return(flux.normalization*DOMEffCorr
           *pow(neutrinoEnergy/1.e5,flux.index));
}

double Flux::IntegratePowerlawFlux(const powerlawFlux flux, double minEnergy, double maxEnergy)
{
    assert(flux.index!=-1.0 && "Special case of E^{-1} not handled");
    double DOMEffCorr=pow(1.189944,-2.-flux.index);
    double intIndex=1.+flux.index;
    double norm=(flux.normalization*DOMEffCorr)/(intIndex*pow(1.e5,flux.index));
    return(norm*(pow(maxEnergy,intIndex)-pow(minEnergy,intIndex)));
}

void Flux::SetNormalization(double arg_normaliztion)
{
	normalization=arg_normalization;
}

void Flux::SetIndex(double arg_normaliztion)
{
	index=arg_index;
}

double Flux::GetNormalization(void)
{
	return normalization;
}

double Flux::GetIndex(void)
{
	return index;
}
