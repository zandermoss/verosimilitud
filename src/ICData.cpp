#include "ICData.h"
#include <iostream>
#include <string>
#include <algorithm>
#include "Tensor.h"
#include <math.h>

ICData::ICData(Tensor *t_eprox_edges, Tensor *t_cosz_edges) {
  cosz_edges = new std::vector<double>(CosZenithBins + 1, 0);
  eprox_edges = new std::vector<double>(EnergyProxyBins + 1, 0);

  for (unsigned int i = 0; i < cosz_edges->size(); i++) {
    (*cosz_edges)[i] = t_cosz_edges->Index(&i);
  }

  for (unsigned int i = 0; i < eprox_edges->size(); i++) {
    (*eprox_edges)[i] = t_eprox_edges->Index(&i);
  }
}

void ICData::OpenCSV(std::string filename) {
  datafile.open(filename.c_str());
  if(!datafile)
    throw std::runtime_error(filename + " not found.");
}

void ICData::ReadCSV(void) {
  std::string value;
  while (datafile.good()) {
    std::getline(datafile, value);
    std::vector<std::string> tokens;
    std::string token;

    while (token != value) {
      token = value.substr(0, value.find_first_of(" "));
      value = value.substr(value.find_first_of(" ") + 1);
      tokens.push_back(token);
    }
    if (tokens.size() < 2) {
      break;
    }
		//Conversion here from zenith to cos(zenith)
    cosz.push_back(cos(atof(tokens[tokens.size() - 1].c_str())));
    eprox.push_back(atof(tokens[tokens.size() - 2].c_str()));
  }
}

void ICData::BinData(std::vector<unsigned int> *cosz_binned,
                     std::vector<unsigned int> *eprox_binned) {
  for (unsigned int i = 0; i < cosz.size(); i++) {
    for (unsigned int bin = 0; bin < cosz_edges->size() - 1; bin++) {
      if ((cosz[i] >= (*cosz_edges)[bin]) &&
          (cosz[i] < (*cosz_edges)[bin + 1])) {
        (*cosz_binned)[bin]++;
      }
    }

    for (unsigned int bin = 0; bin < eprox_edges->size() - 1; bin++) {
      if ((eprox[i] >= (*eprox_edges)[bin]) &&
          (eprox[i] < (*eprox_edges)[bin + 1])) {
        (*eprox_binned)[bin]++;
      }
    }
  }
}

void ICData::BinData(Tensor *binned_data) {
  unsigned int indices[2];
  unsigned int cosbin;
  unsigned int eproxbin;
  double last;

  for (unsigned int i = 0; i < cosz.size(); i++) {
    cosbin = 0;
    eproxbin = 0;

    if (cosz[i] < cosz_edges->front()) {
      cosbin = 0;
    } else if (cosz[i] >= cosz_edges->back()) {
      cosbin = cosz_edges->size() - 2;
    } else {
      while (!((cosz[i] >= (*cosz_edges)[cosbin]) &&
               (cosz[i] < (*cosz_edges)[cosbin + 1]))) {
        cosbin++;
      }
    }

    if (eprox[i] < eprox_edges->front()) {
      eproxbin = 0;
    } else if (eprox[i] >= eprox_edges->back()) {
      eproxbin = eprox_edges->size() - 2;
    } else {
      while (!((eprox[i] >= (*eprox_edges)[eproxbin]) &&
               (eprox[i] < (*eprox_edges)[eproxbin + 1]))) {
        eproxbin++;
      }
    }

    indices[0] = eproxbin;
    indices[1] = cosbin;

    last = binned_data->Index(indices);
    binned_data->SetIndex(indices, last + 1);
  }
}

ICData::~ICData() {
  delete cosz_edges;
  delete eprox_edges;
  datafile.close();
}
