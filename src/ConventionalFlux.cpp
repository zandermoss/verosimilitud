#include "ConventionalFlux.h"


ConventionalFlux::ConventionalFlux()
{	
	H5::H5File file("/Users/marjon/work/IceCube/neutrino_decay/verosimilitud/data/conventional_flux.h5",H5F_ACC_RDONLY);

	std::string year_names[]={"2010","2011"};
	std::string flavor_names[]={"_mu","_tau"};
	std::string anti_names[]={"","_bar"};
	std::string sl="/";
	std::string nu="nu";

	std::string s_area="area";
	std::string s_edges="bin_edges";
	std::string s_detcorr="detector_correction";
	std::string s_iflux="integrated_flux";

	H5::Group* group;

	hsize_t* dims;
	unsigned int* mydims;


	for (unsigned int year=0; year<2; year++)
	{
		group = new H5::Group (file.openGroup((sl+s_detcorr).c_str()));

		std::string data_name=sl+s_detcorr+sl+year_names[year];
	
		H5::DataSet dataset = file.openDataSet((data_name).c_str());


		H5::DataSpace filespace = dataset.getSpace();

		int rank = filespace.getSimpleExtentNdims();

		dims = new hsize_t[rank];
		rank = filespace.getSimpleExtentDims(dims);
	
		mydims = new unsigned int[rank];	
		for (unsigned int i=0; i<rank; i++)
		{
			//std::cout << (unsigned long)(dims[i]) << std::endl;
			mydims[i] = (unsigned int)(dims[i]);
		}
		
		//std::cout << std::endl;
	
		H5::DataSpace mspace1(rank,dims);
	
		detcorr[year]= new Tensor(rank,mydims);
		
		dataset.read(detcorr[year]->GetDataPointer(), H5::PredType::NATIVE_DOUBLE,mspace1,filespace);

	}

	for (unsigned int anti=0; anti<2; anti++)
	{
		group = new H5::Group (file.openGroup((sl+nu+flavor_names[0]
			+anti_names[anti]).c_str()));

		std::string data_name=sl+nu+flavor_names[0]
            +anti_names[anti]+sl+s_iflux;
	
		H5::DataSet dataset = file.openDataSet((data_name).c_str());


		H5::DataSpace filespace = dataset.getSpace();

		int rank = filespace.getSimpleExtentNdims();

		dims = new hsize_t[rank];
		rank = filespace.getSimpleExtentDims(dims);
	
		mydims = new unsigned int[rank];	
		for (unsigned int i=0; i<rank; i++)
		{
			//std::cout << (unsigned long)(dims[i]) << std::endl;
			mydims[i] = (unsigned int)(dims[i]);
		}
		
		//std::cout << std::endl;
	
		H5::DataSpace mspace1(rank,dims);
	
		flux[anti]= new Tensor(rank,mydims);
		
		dataset.read(flux[anti]->GetDataPointer(), H5::PredType::NATIVE_DOUBLE,mspace1,filespace);

	}

		
	delete dims;
	delete mydims;
	delete group;
	file.close();
}


ConventionalFlux::~ConventionalFlux()
{
	for (unsigned int i=0; i<2; i++)
	{
					delete detcorr[i];
					delete flux[i];
	}
}


Tensor* ConventionalFlux::GetDetCorr(unsigned int * index)
{
	return detcorr[index[0]];
}

Tensor* ConventionalFlux::GetFlux(unsigned int * index)
{
	return flux[index[0]];
}
