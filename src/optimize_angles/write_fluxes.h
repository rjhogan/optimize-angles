#include <string>
#include <adept_arrays.h>

void
write_fluxes(const std::string& file_name, const std::string& method_id,
	     const adept::Vector& flux_up_toa, const adept::Vector& flux_dn_surf,
	     const adept::Matrix& heating_rate,
	     const adept::Vector& mu, const adept::Vector& weight,
	     int argc, const char** argv);
