#ifndef __TENSOR_H_INCLUDED__
#define __TENSOR_H_INCLUDED__


//template <class T>
class Tensor
{
	public:
		Tensor(unsigned int argrank, unsigned int* argdims);
		Tensor(unsigned int argrank, unsigned int* argdims,double initval);
		~Tensor(void);
		double Index(unsigned int* indicies);
		void SetIndex(unsigned int* indicies, double value);
		double* GetDataPointer(void);	
		unsigned int GetDataLength(void);	
		unsigned int GetRank(void);
		void GetDims(unsigned int* dims); 

	private:
		double* data;
		unsigned int rank;
		unsigned int* dims;
		unsigned int* strides;
};


#endif // __TENSOR_H_INCLUDED__
