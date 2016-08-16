#ifndef __CONVENTIONALFLUX_H_INCLUDED__
#define __CONVENTIONALFLUX_H_INCLUDED__



#include "Tensor.h"
#include "H5Cpp.h"
#include <string>
#include <iostream>

class ConventionalFlux
{
	public:
		ConventionalFlux();
		~ConventionalFlux();

		Tensor* GetDetCorr(unsigned int * index);
		Tensor* GetFlux(unsigned int * index);

	private:
		Tensor* detcorr[2];
		Tensor* flux[4];

};


#endif // __CONVENTIONALFLUX_H_INCLUDED__

