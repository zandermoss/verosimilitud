#include "EffectiveArea.h"

EffectiveArea::EffectiveArea(std::string effective_area_path) {


	std::string edge_data_names[] = {"true_energy_bin_edges", "costh_bin_edges", "proxy_energy_bin_edges"};
	std::string anti_data_names[] = {"neutrino_response_array", "antineutrino_response_array"};

	hsize_t *dims;
	unsigned int *mydims;

	for (unsigned int effi=0; effi<neff; effi++){
		std::string fname = effective_area_path+"nufsgen_dom_eff_"
													+efficiencies[effi]+"_response_array.h5";

		H5::H5File file(fname.c_str(), H5F_ACC_RDONLY);
	
			// Now read in the bin edges.
			if (effi==0){
				for (unsigned int edgedim = 0; edgedim < 3; edgedim++) {
					H5::DataSet dataset = file.openDataSet((edge_data_names[edgedim].c_str()));
					H5::DataSpace filespace = dataset.getSpace();
					int rank = filespace.getSimpleExtentNdims();
					dims = new hsize_t[rank];
					rank = filespace.getSimpleExtentDims(dims);
					mydims = new unsigned int[rank];
					for (int i = 0; i < rank; i++) {
						// std::cout << (unsigned long)(dims[i]) << std::endl;
						mydims[i] = (unsigned int)(dims[i]);
					}
					// std::cout << std::endl;
					H5::DataSpace mspace1(rank, dims);
					edges[edgedim] = new Tensor(rank, mydims);
					dataset.read(edges[edgedim]->GetDataPointer(),H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
				}
			}

			for (unsigned int anti = 0; anti < 2; anti++) {

				H5::DataSet dataset =
						file.openDataSet(anti_data_names[anti].c_str());

				H5::DataSpace filespace = dataset.getSpace();

				int rank = filespace.getSimpleExtentNdims();

				dims = new hsize_t[rank];
				rank = filespace.getSimpleExtentDims(dims);

				mydims = new unsigned int[rank];
				for (int i = 0; i < rank; i++) {
					// std::cout << (unsigned long)(dims[i]) << std::endl;
					mydims[i] = (unsigned int)(dims[i]);
				}

				// std::cout << std::endl;

				H5::DataSpace mspace1(rank, dims);

				areas[effi][anti] = new Tensor(rank, mydims);

				dataset.read(areas[effi][anti]->GetDataPointer(),
										 H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
		}
	
		delete dims;
		delete mydims;
		file.close();
	}
}


EffectiveArea::~EffectiveArea() {
	for (unsigned int i = 0; i < neff; i++) {
		for (unsigned int j = 0; j < 2; j++) {
			delete areas[i][j];
		}
	}
	for (unsigned int i=0; i<3; i++){
		delete edges[i];
	}
}

Tensor *EffectiveArea::GetArea(unsigned int *index) {
	return areas[index[0]][index[1]];
}

Tensor *EffectiveArea::GetEdge(unsigned int index) {
	return edges[index];
}

double EffectiveArea::GetLivetime(void) const {
	return livetime;
}

unsigned int EffectiveArea::GetNEff(void){
	return neff;
}

double EffectiveArea::GetEff(unsigned int index){
	return eff_vals[index]; 
}



