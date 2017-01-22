#ifndef __CONVENTIONALFLUX_H_INCLUDED__
#define __CONVENTIONALFLUX_H_INCLUDED__

#include "Tensor.h"
#include "H5Cpp.h"
#include <string>
#include <iostream>

class ConventionalFlux {
public:
  ConventionalFlux(std::string flux_path);
  ~ConventionalFlux();

  Tensor *GetFlux(unsigned int *index);

private:
  Tensor *flux[4];
};

#endif // __CONVENTIONALFLUX_H_INCLUDED__
