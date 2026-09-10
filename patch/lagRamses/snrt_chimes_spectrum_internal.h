#ifndef SNRT_CHIMES_SPECTRUM_INTERNAL_H
#define SNRT_CHIMES_SPECTRUM_INTERNAL_H
#include "snrt_band_spectrum.h"
// Private shared representation; public callers own only an opaque handle.
namespace snrt_chimes_detail {
struct Bank {
  std::vector<int> reaction,shell;
  std::vector<double> binding,sigma,energy;
  std::vector<snrt_band::Grid<128>> grids;
};
struct MoleculeBank {
  std::vector<int> reaction;
  std::vector<double> absorption,dissociation;
};
}
#endif
