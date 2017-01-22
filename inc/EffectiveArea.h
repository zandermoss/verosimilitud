#ifndef __EFFECTIVEAREA_H_INCLUDED__
#define __EFFECTIVEAREA_H_INCLUDED__

#include "Tensor.h"
#include "H5Cpp.h"
#include <string>
#include <iostream>
#include <vector>

class EffectiveArea {
public:
  EffectiveArea(std::string effective_area_path);
  ~EffectiveArea();

  Tensor *GetArea(unsigned int *index);
  Tensor *GetEdge(unsigned int index);
  double GetLivetime(void) const;
	unsigned int GetNEff(void);
	double GetEff(unsigned int);

private:

	unsigned int neff = 5;

	//Hard-coding the lifetime from the 2011 IC86 run
	//value extracted from the effective_area.h5 file
	//from Chris Weaver's data release used in the two
	//year runs of this code.
	double livetime = 2.969856e7;

	std::vector<double> eff_vals = {0.90, 0.95, 0.99, 1.089, 1.1979};
	std::vector<std::string> efficiencies = {"0_90","0_95","nominal","1_089","1_1979"};

  Tensor *areas[5][2];
  Tensor *edges[3];

};

#endif // __EFFECTIVEAREA_H_INCLUDED__
