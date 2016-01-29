#include <math.h>
#include "Tools.h"

double LogPoissonProbability(unsigned int x,double mu) // simple poisson function
{
    double logfact=0;
    for (unsigned int j=1; j<=x; j++)
    {
        logfact+=log((double)j);
    }
    // Use product instead? Numerically stable? Time it?

    return ((double)x)*log(mu)-mu-logfact;
}


double LogGaussianProbability(double x,double mu,double sigma) // gaussian probability
{
    return -log(sigma) -log2pi -0.5*pow((x-mu)/sigma,2);
}


