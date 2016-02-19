#ifndef __EFFECTIVEAREA_H_INCLUDED__
#define __EFFECTIVEAREA_H_INCLUDED__



#include "Tensor.h"
#include "H5Cpp.h"
#include <string>
#include <iostream>

class EffectiveArea
{
	public:
		EffectiveArea();
		~EffectiveArea();

		Tensor* GetArea(unsigned int * index);
		Tensor* GetEdge(unsigned int * index);
		double GetLivetime(unsigned int index);

	private:
		Tensor* areas[2][2][2];
		Tensor* edges[2][2][2][3];
		double livetime[2];

};


#endif // __EFFECTIVEAREA_H_INCLUDED__

