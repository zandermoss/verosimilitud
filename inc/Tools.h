#ifndef __TOOLS_H_INCLUDED__
#define __TOOLS_H_INCLUDED__

#include <math.h>
#include <vector>

// Constants
const double pi = 3.141592;
const double log2pi = log(2 * pi); // precompute

// Functions
double LogPoissonProbability(unsigned int x,
                             double mu); // simple poisson function
double LogGaussianProbability(double x, double mu,
                              double sigma); // gaussian probability

void PrintArray(std::vector<double> array);
void PrintMatrix(std::vector<std::vector<double>> matrix);

#endif // __TOOLS_H_INCLUDED__
