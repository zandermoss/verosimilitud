#include "Tensor.h"
#include <iostream>

//template <class T>
Tensor::Tensor(unsigned int argrank, unsigned int* argdims)
{
	rank=argrank;
	dims = new unsigned int[rank];
	strides = new unsigned int[rank]; 
	unsigned int datasize=1;
	
	for (unsigned int i=0;i<rank;i++)
	{
		dims[i]=argdims[i];
		strides[i]=datasize;
		datasize*=argdims[rank-(i+1)];
	}

	data = new double[datasize];
}


Tensor::Tensor(unsigned int argrank, unsigned int* argdims,double initval)
{
	rank=argrank;
	dims = new unsigned int[rank];
	strides = new unsigned int[rank]; 
	unsigned int datasize=1;
	
	for (unsigned int i=0;i<rank;i++)
	{
		dims[i]=argdims[i];
		strides[i]=datasize;
		datasize*=argdims[rank-(i+1)];
	}

	data = new double[datasize];

	for (unsigned int i=0; i<datasize; i++)
	{
		data[i]=initval;
	}

}




//template <class T>
Tensor::~Tensor(void)
{
	//delete data[];
	//delete dims[]; 
	//delete strides[]; 
	delete data;
	delete dims; 
	delete strides; 
}


//template <class T>
double Tensor::Index(unsigned int* indices)
{
	unsigned int index=0;
	for (unsigned int i=0;i<rank;i++)
	{
		index+=strides[i]*indices[rank-(i+1)];
		
	}

	return data[index];
}

//template <class T>
void Tensor::SetIndex(unsigned int* indices, double value)
{
	unsigned int index=0;
	for (unsigned int i=0;i<rank;i++)
	{
		index+=strides[i]*indices[rank-(i+1)];
		
	}

	data[index] = value;

}

//template <class T>
double* Tensor::GetDataPointer(void)
{
	return data;
}

unsigned int Tensor::GetDataLength(void)
{
	unsigned int datalen=1;
    for(unsigned int i=0; i<rank; i++)
    {
        datalen*=dims[i];
    }

	return datalen;
}

//template <class T>
unsigned int Tensor::GetRank(void)
{
	return rank;
}


//template <class T>
void Tensor::GetDims(unsigned int* argdims)
{
	for(unsigned int i=0; i<rank; i++)
	{
		argdims[i]=dims[i];
	} 
}
