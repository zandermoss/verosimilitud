#include "EffectiveArea.h"

EffectiveArea::EffectiveArea(std::string effective_area_path) {
  H5::H5File file(effective_area_path.c_str(),
                  H5F_ACC_RDONLY);

  std::string year_names[] = {"2010", "2011"};
  std::string flavor_names[] = {"_mu", "_tau"};
  std::string anti_names[] = {"", "_bar"};
  std::string sl = "/";
  std::string nu = "nu";

  std::string s_area = "area";
  std::string s_edges = "bin_edges";
  std::string s_dims[3] = {"_0", "_1", "_2"};

  H5::Group *group;
  H5::Attribute *attr;
  H5::DataType *type;

  hsize_t *dims;
  unsigned int *mydims;

  for (unsigned int year = 0; year < 2; year++) {
    group = new H5::Group(file.openGroup((sl + year_names[year]).c_str()));
    attr = new H5::Attribute(
        group->openAttribute("experimental_livetime(seconds)"));
    type = new H5::DataType(attr->getDataType());

    attr->read(*type, &livetime[year]);

    for (unsigned int flavor = 0; flavor < 2; flavor++) {
      for (unsigned int anti = 0; anti < 2; anti++) {

        std::string data_name = sl + year_names[year] + sl + nu +
                                flavor_names[flavor] + anti_names[anti];

        H5::DataSet dataset =
            file.openDataSet((data_name + sl + s_area).c_str());

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

        areas[year][flavor][anti] = new Tensor(rank, mydims);

        dataset.read(areas[year][flavor][anti]->GetDataPointer(),
                     H5::PredType::NATIVE_DOUBLE, mspace1, filespace);

        // Now read in the bin edges.
        for (unsigned int edgedim = 0; edgedim < 3; edgedim++) {

          H5::DataSet dataset = file.openDataSet(
              (data_name + sl + s_edges + s_dims[edgedim]).c_str());

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

          edges[year][flavor][anti][edgedim] = new Tensor(rank, mydims);

          dataset.read(edges[year][flavor][anti][edgedim]->GetDataPointer(),
                       H5::PredType::NATIVE_DOUBLE, mspace1, filespace);
        }
      }
    }
  }

  delete dims;
  delete mydims;
  delete attr;
  delete group;
  delete type;
  file.close();
}


EffectiveArea::~EffectiveArea() {
  for (unsigned int i = 0; i < 2; i++) {
    for (unsigned int j = 0; j < 2; j++) {
      for (unsigned int k = 0; k < 2; k++) {
        delete areas[i][j][k];

        for (unsigned int l = 0; l < 3; l++) {
          delete edges[i][j][k][l];
        }
      }
    }
  }
}

Tensor *EffectiveArea::GetArea(unsigned int *index) {
  return areas[index[0]][index[1]][index[2]];
}

Tensor *EffectiveArea::GetEdge(unsigned int *index) {
  return edges[index[0]][index[1]][index[2]][index[3]];
}

double EffectiveArea::GetLivetime(unsigned int index) const {
  return livetime[index];
}
