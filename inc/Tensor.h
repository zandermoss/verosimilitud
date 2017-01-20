#ifndef __TENSOR_H_INCLUDED__
#define __TENSOR_H_INCLUDED__

// template <class T>
class Tensor {
public:
  Tensor(unsigned int argrank, unsigned int *argdims);
  Tensor(unsigned int argrank, unsigned int *argdims, double initval);
  ~Tensor(void);
  double Index(unsigned int *indicies);
  void SetIndex(unsigned int *indicies, double value);
  double *GetDataPointer(void) const;
  unsigned int GetDataLength(void) const;
  unsigned int GetRank(void) const;
  void GetDims(unsigned int *dims) const;
	Tensor operator+(const Tensor&) const;
	Tensor operator*(const double) const; //scalar multiplication
	friend Tensor operator*(const double,const Tensor&);
private:
  double *data;
  unsigned int rank;
  unsigned int *dims;
  unsigned int *strides;
	unsigned int datasize;
};

#endif // __TENSOR_H_INCLUDED__
