#include "OutputDataFile.h"
#include "write_fluxes.h"
#include "write_standard_attributes.h"

using namespace adept;

void
write_fluxes(const std::string& file_name, const std::string& method_id,
	     const Vector& flux_up_toa, const Vector& flux_dn_surf,
	     const Matrix& heating_rate,
	     const Vector& mu, const Vector& weight,
	     int argc, const char** argv) {
  OutputDataFile file(file_name);
  file.define_dimension("column", flux_up_toa.size());
  file.define_dimension("layer", heating_rate.size(0));
  file.define_dimension("node", mu.size());
  file.define_variable("mu", DOUBLE, "node");
  file.write_long_name("Cosine of zenith angle", "mu");
  file.write_missing_value(-1.0, "mu");
  file.define_variable("weight", DOUBLE, "node");
  file.write_long_name("Weight", "weight");
  file.write_missing_value(0.0, "weight");
  file.define_variable("flux_up_toa", DOUBLE, "column");
  file.write_long_name("Top-of-atmosphere upwelling longwave flux", "flux_up_toa");
  file.write_units("W m-2", "flux_up_toa");
  file.define_variable("flux_dn_surf", DOUBLE, "column");
  file.write_long_name("Surface downwelling longwave flux", "flux_dn_surf");
  file.write_units("W m-2", "flux_dn_surf");
  file.define_variable("heating_rate", DOUBLE, "column", "layer");
  file.write_long_name("Heating rate", "heating_rate");
  file.write_units("K d-1", "heating_rate");
  write_standard_attributes(file,
    "Longwave fluxes and heating rates");
  file.append_history(argc, argv);
  file.write(method_id, "quadrature_method_id");

  file.end_define_mode();

  file.write(mu, "mu");
  file.write(weight, "weight");
  file.write(flux_up_toa, "flux_up_toa");
  file.write(flux_dn_surf, "flux_dn_surf");
  file.write(heating_rate.T(), "heating_rate");
  file.close();
}
