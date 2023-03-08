#include "DataFile.h"
#include "AtmosProfiles.h"
#include "Error.h"
#include "write_fluxes.h"

using namespace adept;

int
main(int argc, const char* argv[])
{

  // Read configuration information from command-line and first file
  // on command-line
  DataFile config(argc, argv);
  std::string input_file, quad_file, output_file;
  int iorder;
  config.read(input_file, "1");
  config.read(quad_file,  "2");
  config.read(iorder,     "3");
  config.read(output_file,"4");

  // Load quadrature information
  DataFile quad(quad_file);
  intVector orders;
  quad.read(orders, "order");
  intVector ordind = find(orders == iorder);
  if (ordind.empty()) {
    ERROR << "Quadrature order " << iorder << " not found in " << quad_file;
    THROW(PARAMETER_ERROR);
  }
  Vector mu, wt;
  quad.read(mu, "mu", ordind(0));
  quad.read(wt, "weight", ordind(0));
  std::string method_id;
  quad.read(method_id, DATA_FILE_GLOBAL_SCOPE, "quadrature_method_id");
  
  // Load multiple profiles and perform radiative transfer
  AtmosProfiles profile(input_file);
  int nlay = profile.nlayer();
  int ncol = profile.ncolumn();

  Vector flux_surf(ncol), flux_toa(ncol);
  Matrix heating_rate(nlay,ncol);
  profile.radiative_transfer_lw(mu(range(0,iorder-1)), wt(range(0,iorder-1)),
				flux_surf, flux_toa, heating_rate);
  // Write the results
  write_fluxes(output_file, method_id, flux_toa, flux_surf, 
	       heating_rate, mu(range(0,iorder-1)), wt(range(0,iorder-1)), argc, argv);
}
