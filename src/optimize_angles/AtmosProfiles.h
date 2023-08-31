// AtmosProfiles.h - Structure to store atmospheric profiles -*- C++ -*-
//
// Copyright (C) 2022- ECMWF.
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
//
// In applying this licence, ECMWF does not waive the privileges and immunities
// granted to it by virtue of its status as an intergovernmental organisation
// nor does it submit to any jurisdiction.
//
// Author:  Robin Hogan
// Email:   r.j.hogan@ecmwf.int

#ifndef AtmosProfiles_H
#define AtmosProfiles_H

#include <string>

#include <adept_arrays.h>

class AtmosProfiles {

public:
  AtmosProfiles(const std::string& file_name) {
    read(file_name);
  }

  void read(const std::string& file_name);
  
  template <bool IsActive>
  void radiative_transfer_lw(const Array<1,Real,IsActive>& mu,
			     const Array<1,Real,IsActive>& weight,
			     Array<1,Real,IsActive> flux_surf,
			     Array<1,Real,IsActive> flux_toa,
			     Array<2,Real,IsActive> heating_rate) const;

  int ncolumn() const { return optical_depth.size(0); }
  int nlayer()  const { return optical_depth.size(1); }
  int nspec()   const { return optical_depth.size(2); }

  Matrix heating_rate_weight() const;
  Matrix pressure() const { return pressure_hl; }
  
protected:
  // DATA
  
  // Pressure at half levels counting down from top-of-atmosphere
  // (Pa), dimensioned [column,level]
  Matrix pressure_hl;

  // Planck function at half levels in each spectral interval (W m-2),
  // dimensioned [column,level,spec]
  Array3D planck_hl;

  // Planck function at surface in each spectral interval (W m-2),
  // dimensioned [column,spec]
  Matrix planck_surf;

  // Layer optical depth, dimensioned [column,layer,spec]
  Array3D optical_depth;

};

#endif
