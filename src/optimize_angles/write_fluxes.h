// write_fluxes.h - Write fluxes to a NetCDF file -*- C++ -*-
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

#include <string>
#include <adept_arrays.h>

void
write_fluxes(const std::string& file_name, const std::string& method_id,
	     const adept::Vector& flux_up_toa, const adept::Vector& flux_dn_surf,
	     const adept::Matrix& heating_rate,
	     const adept::Vector& mu, const adept::Vector& weight,
	     int argc, const char** argv);
