#ifndef __TOOLS_H_INCLUDED__
#define __TOOLS_H_INCLUDED__

#include <math.h>

//Constants
const double pi=3.141592;
const double log2pi=log(2*pi); //precompute

//Functions
double LogPoissonProbability(unsigned int x, double mu); // simple poisson function
double LogGaussianProbability(double x,double mu,double sigma); // gaussian probability

#endif // __TOOLS_H_INCLUDED__ 







