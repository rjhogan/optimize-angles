// -*- C++ -*-

#ifndef AtmosProfile_H
#define AtmosProfile_H

#include <string>

#include <adept_arrays.h>

class AtmosProfile {

public:
  AtmosProfile(const std::string& file_name) {
    read(file_name);
  }

  void read(const std::string& file_name);
  
  template <bool IsActive>
  void radiative_transfer_lw(const Array<1,Real,IsActive>& mu,
			     const Array<1,Real,IsActive>& weight,
			     Array<1,Real,IsActive> flux_surf,
			     Array<1,Real,IsActive> flux_toa,
			     Array<2,Real,IsActive> heating_rate) const;

  int nlayer()  const { return optical_depth.size(0); }
  int nspec()   const { return optical_depth.size(1); }
  int ncolumn() const { return optical_depth.size(1); }

  Matrix heating_rate_weight() const;
  
protected:
  // DATA
  
  // Pressure at half levels counting down from top-of-atmosphere (Pa)
  Vector pressure_hl;

  // Planck function at half levels in each spectral interval (W m-2),
  // dimensioned [lev,spec]
  Matrix planck_hl;

  // Planck function at surface in each spectral interval (W m-2)
  Vector planck_surf;

  // Layer optical depth, dimensioned [lay,spec]
  Matrix optical_depth;

};

#endif
