#include "Tensor.h"
#include <iostream>
#include <stdexcept>

// template <class T>
Tensor::Tensor(unsigned int argrank, unsigned int *argdims) {
  rank = argrank;
  dims = new unsigned int[rank];
  strides = new unsigned int[rank];
  datasize = 1;

  for (unsigned int i = 0; i < rank; i++) {
    dims[i] = argdims[i];
    strides[i] = datasize;
    datasize *= argdims[rank - (i + 1)];
  }

  data = new double[datasize];
}

Tensor::Tensor(unsigned int argrank, unsigned int *argdims, double initval) {
  rank = argrank;
  dims = new unsigned int[rank];
  strides = new unsigned int[rank];
  datasize = 1;

  for (unsigned int i = 0; i < rank; i++) {
    dims[i] = argdims[i];
    strides[i] = datasize;
    datasize *= argdims[rank - (i + 1)];
  }

  data = new double[datasize];

  for (unsigned int i = 0; i < datasize; i++) {
    data[i] = initval;
  }
}

// template <class T>
Tensor::~Tensor(void) {
  // delete data[];
  // delete dims[];
  // delete strides[];
  delete data;
  delete dims;
  delete strides;
}

// template <class T>
double Tensor::Index(unsigned int *indices) {
  unsigned int index = 0;
  for (unsigned int i = 0; i < rank; i++) {
    index += strides[i] * indices[rank - (i + 1)];
  }

  return data[index];
}

// template <class T>
void Tensor::SetIndex(unsigned int *indices, double value) {
  unsigned int index = 0;
  for (unsigned int i = 0; i < rank; i++) {
    index += strides[i] * indices[rank - (i + 1)];
  }

  data[index] = value;
}

// template <class T>
double *Tensor::GetDataPointer(void) const { return data; }

unsigned int Tensor::GetDataLength(void) const {
  unsigned int _datalen = 1;
  for (unsigned int i = 0; i < rank; i++) {
    _datalen *= dims[i];
  }

  return _datalen;
}

// template <class T>
unsigned int Tensor::GetRank(void) const { return rank; }

// template <class T>
void Tensor::GetDims(unsigned int *argdims) const {
  for (unsigned int i = 0; i < rank; i++) {
    argdims[i] = dims[i];
  }
}


/*
	Operator overloading for easy tensor summation and
	scalar multiplication. Both operations are commutative,
	and have been so defined below. 
*/

Tensor Tensor::operator+(const Tensor& t) const{	
	//Check for matching rank and dimensions
	if (t.GetRank() != rank){
		throw std::runtime_error("Summand ranks do not match.");
	}
	unsigned int tdims[rank];
	t.GetDims(tdims);
	for (unsigned int i=0; i<rank; i++){
		if (tdims[i] != dims[i]){
			throw std::runtime_error("Summand dimensions do not match.");
		}
	}	
	double* t_data = t.GetDataPointer();
	Tensor sum(rank,dims);
	double* sum_data = sum.GetDataPointer();
	for (unsigned int i=0; i<datasize; i++){
		sum_data[i] = t_data[i] + data[i];
	}	
	return sum;
}

Tensor Tensor::operator*(const double k) const {	
	Tensor product(rank,dims);
	double* prod_data = product.GetDataPointer();
	for (unsigned int i=0; i<datasize; i++){
		prod_data[i] = k*data[i];
	}	
	return product;
}

Tensor operator*(const double k, const Tensor& t) {
	return t*k;
}

