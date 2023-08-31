// optimize_angles.cpp - Program for optimizing angles in longwave radiative transfer
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

#include <iostream>
#include <sstream>

#include <omp.h>

#include <adept_arrays.h>
#include <adept_optimize.h>

#include "DataFile.h"
#include "OutputDataFile.h"
#include "OptimizeAngles.h"
#include "AtmosProfile.h"
#include "AtmosProfiles.h"
#include "Error.h"
#include "write_standard_attributes.h"
#include "write_fluxes.h"

using namespace adept;

int
main(int argc, const char** argv)
{
  
  if (!adept::have_linear_algebra()) {
    std::cout << "Adept compiled without linear-algebra support: minimizer not available\n";
    return 0;
  }
  
  adept::set_array_print_style(PRINT_STYLE_MATLAB);
  adept::internal::array_contiguous_separator = ", ";
  adept::internal::vector_separator = ", ";

  // CONFIGURATION

  // Read configuration information from command-line and first file
  // on command-line
  DataFile config(argc, argv);

  std::string log_level;
  if (config.read(log_level, "log_level")) {
    set_log_level(log_level);
  }

  int max_order = 4;
  config.read(max_order, "max_order");

  intVector orders;
  int norder = max_order;
  if (config.read(orders, "orders")) {
    norder = orders.size();
    max_order = orders(end);
  }
  else {
    orders = range(1,max_order);
  }
  
  std::string method_id;
  if (!config.read(method_id, "method_id")) {
    method_id = "optimization";
  }

  /*
  // Optional forced ratio between nodes in order-2 quadrature
  int integer_ratio = -1;
  if (config.read(integer_ratio, "integer_ratio")) {
    if (integer_ratio < 1) {
      ERROR << "\"integer_ratio\" must be 2 or greater";
      THROW(PARAMETER_ERROR);
    }
  }
  */

  Matrix tmp_ratios;
  intMatrix integer_ratios;
  if (config.read(tmp_ratios, "integer_ratios")) {
    if (tmp_ratios.size(0) != norder
	|| tmp_ratios.size(1) != max_order-1) {
      ERROR << "integer_ratios must be dimensioned [" << norder
	    << "," << max_order-1 << "]";
      THROW(PARAMETER_ERROR);
    }
    integer_ratios = tmp_ratios;
    if (any(integer_ratios != tmp_ratios)) {
      ERROR << "At least one of integer_ratios is not an integer";
      THROW(PARAMETER_ERROR);
    }
  }
  
  // Modify search path
  /*
  std::string mod_path;
  if (config.read(mod_path, "prepend_path")) {
    prepend_search_directory(mod_path);
  }
  if (config.read(mod_path, "append_path")) {
    append_search_directory(mod_path);
  }
  */
  
  // Allow debugging
  //  set_trace_exceptions(true);

  std::string input_file;
  std::string output_file;
  if (!config.read(output_file, "2") || !config.read(input_file, "1"))  {
    ERROR << "Usage: " << argv[0] << " [arguments] <input_profiles>.nc <output_coefficients>.nc";
    THROW(PARAMETER_ERROR);
  }

  //  if (!config.read(input_file, "1")) {
  //    ERROR << "\"input=<file>.nc\" not provided on command line";
  //    THROW(PARAMETER_ERROR);
  //  }
  
  // Load one profile and optimize on spectral fluxes
  //AtmosProfile  profile(input_file);
  // Load multiple profiles and optimize on broadband fluxes
  AtmosProfiles profile(input_file);

  Real flux_weight = 0.02;
  config.read(flux_weight, "flux_weight");
  
  Real broadband_weight = 1.0;
  config.read(broadband_weight, "broadband_weight");

  int nlayer = profile.nlayer();
  int ncol  = profile.ncolumn();

  // Angles and weights of reference calculation
  int nmu_ref;
  Vector mu_ref, weight_ref;

  bool save_fluxes = false;
  config.read(save_fluxes, "save_fluxes");
  
  std::string ref_quadrature, ref_method_id;
  if (config.read(ref_quadrature, "ref_quadrature")) {
    // Read the final entry of the named data file, presumably the one
    // with the most number of nodes and therefore the most accurate
    DataFile ref_quad(ref_quadrature);
    intVector ref_quad_size = ref_quad.size("mu");
    int nmu_ref_entries = ref_quad_size(0);
    int ref_index = nmu_ref_entries-1;
    config.read(ref_index, "ref_index");
    if (ref_index >= nmu_ref_entries) {
      ERROR << "ref_index exceeds the number of entries in the quadrature file";
      THROW(PARAMETER_ERROR);
    }
    ref_quad.read(mu_ref, "mu", ref_index);
    ref_quad.read(weight_ref, "weight", ref_index);
    ref_quad.read(ref_method_id, DATA_FILE_GLOBAL_SCOPE, "quadrature_method_id");
    ref_quad.close();
    nmu_ref = mu_ref.size();
    if (nmu_ref < 8) {
      WARNING << "Reference calculation using only " << nmu_ref
	      << " angles: not likely to be accurate";
      ENDWARNING;
    }
  }
  else {
    // Or use this many evenly-spaced nodes
    nmu_ref = 512;
    mu_ref = linspace(1.0/(2.0*nmu_ref), 1.0-1.0/(2.0*nmu_ref), nmu_ref);
    weight_ref.resize(nmu_ref);
    weight_ref = 1.0/nmu_ref;
    ref_quadrature = "evenly-spaced nodes";
  }

  std::string init_quadrature;
  intVector init_orders;
  Matrix init_mu;
  Matrix init_wt;
  bool use_init_quad = false;
  if (config.read(init_quadrature, "init_quadrature")) {
    DataFile init_quad(init_quadrature);
    init_quad.read(init_orders, "order");
    init_quad.read(init_mu, "mu");
    init_quad.read(init_wt, "weight");
    use_init_quad = true;
  }

  std::string prior_quadrature;
  intVector prior_orders;
  Matrix prior_mu;
  Matrix prior_wt;
  bool use_prior_quad = false;
  if (config.read(prior_quadrature, "prior_quadrature")) {
    DataFile prior_quad(prior_quadrature);
    prior_quad.read(prior_orders, "order");
    prior_quad.read(prior_mu, "mu");
    prior_quad.read(prior_wt, "weight");
    use_prior_quad = true;
  }

  Vector prior_weight;
  if (use_prior_quad) {
    if (!config.read(prior_weight, "prior_weight")) {
      ERROR << "prior_quadrature must be accompanied by prior_weight";
      THROW(PROCESSING_ERROR);
    }
    if (prior_weight.size() == 1 && orders.size() > 1) {
      // If user provides only a single prior weight, expand it up to
      // the number of orders to be used
      Real pw = prior_weight(0);
      prior_weight.resize(orders.size());
      prior_weight = pw;
    }
    else if (prior_weight.size() != orders.size()) {
      ERROR << "prior_weight must be same size as orders";
      THROW(PROCESSING_ERROR);
    }
  }
  else {
    prior_weight.resize(orders.size());
    prior_weight = 0.0;
  }

  
  
  Vector y_ref((nlayer+2)*ncol+2);
  // Soft links to y_ref
  Matrix heating_rate_ref = y_ref(range(0,ncol*nlayer-1)).reshape(nlayer,ncol);
  Vector flux_surf_ref    = y_ref(range(ncol*nlayer,ncol*(nlayer+1)-1));
  Vector flux_toa_ref     = y_ref(range(ncol*(nlayer+1),ncol*(nlayer+2)-1));
  //Real& flux_surf_bb_ref  = y_ref(end-1);
  //Real& flux_toa_bb_ref   = y_ref(end);
  
  LOG << "Performing radiative transfer on reference profiles with "
      << nmu_ref << " angles (from " << ref_quadrature << ")\n";
  profile.radiative_transfer_lw(mu_ref, weight_ref, flux_surf_ref,
				flux_toa_ref, heating_rate_ref);
  if (save_fluxes) {
    write_fluxes("reference_fluxes.nc", ref_method_id, 
		 flux_toa_ref, flux_surf_ref,
		 heating_rate_ref, mu_ref, weight_ref,
		 argc, argv);
  }
  y_ref(end-1) = sum(flux_surf_ref);
  y_ref(end)   = sum(flux_toa_ref);  
  
  //LOG << "  Surface downwelling flux: " << sum(flux_surf_ref) << " W m-2\n";
  //LOG << "  TOA upwelling flux:       " << sum(flux_toa_ref)  << " W m-2\n";
  /*
  LOG << "  Layer, central pressure (hPa), heating rate (K d-1):\n";
  for (int ilay = 0; ilay < nlayer; ++ilay) {
    LOG << "    " << ilay << ", "
	<< 0.01*0.5*(profile.pressure_hl(ilay)+profile.pressure_hl(ilay+1))
	<< ", " << sum(heating_rate_ref[ilay]) << "\n";
  }
  */

  Vector y_wt((nlayer+2)*ncol+2);
  Matrix heating_rate_wt = y_wt(range(0,ncol*nlayer-1)).reshape(nlayer,ncol);
  Vector flux_surf_wt    = y_wt(range(ncol*nlayer,ncol*(nlayer+1)-1));
  Vector flux_toa_wt     = y_wt(range(ncol*(nlayer+1),ncol*(nlayer+2)-1));

  //  Vector hr_wt = (sqrt(profile.pressure_hl(range(1,end)))-sqrt(profile.pressure_hl(range(0,end-1))))
  //    / sqrt(profile.pressure_hl(end));
  //  spread<1>(hr_wt, ncol);
  heating_rate_wt = profile.heating_rate_weight();
  flux_surf_wt = flux_weight;
  flux_toa_wt  = flux_weight;
  y_wt(end-1) = flux_weight*broadband_weight;
  y_wt(end)   = flux_weight*broadband_weight;

  Matrix mu_save(norder, max_order);
  Matrix wt_save(norder, max_order);
  mu_save = -1.0;
  wt_save = 0.0;

  // More efficient use of parallel threads to start with the slowest
  // calculations
  std::cout << "Processing quadrature orders from highest to lowest for most efficient use of parallel threads\n";
#pragma omp parallel for schedule(dynamic,1)
  for (int iang = norder-1; iang >= 0; --iang) {
    int iorder = orders(iang);

#pragma omp critical
    {
      if (iorder > 4) {
	std::cout << "[Thread " << omp_get_thread_num() << " optimizing "
		  << iorder << " angles: this could take some time...]\n";
      }
      else {
	std::cout << "[Thread " << omp_get_thread_num() << " optimizing "
		  << iorder << " angles]\n";
      }
    }

    OptimizeAngles<decltype(profile)> opt_problem(profile, iorder, y_ref, y_wt);

    // Levenberg-Marquardt is considerably slower (including needing
    // more iterations)
    //Minimizer minimizer(MINIMIZER_ALGORITHM_LEVENBERG_MARQUARDT);
    Minimizer minimizer(MINIMIZER_ALGORITHM_LIMITED_MEMORY_BFGS);
    minimizer.set_converged_gradient_norm(0.0000000000001);
    minimizer.set_levenberg_damping_start(0.25);
    minimizer.ensure_updated_state(2);
    minimizer.set_max_iterations(20000);

    int nx;
    /*
    if (integer_ratio > 1) {
      if (iorder != 2) {
	ERROR << "\"integer_ratio\" can only be used when optimizing two angles";
	THROW(PROCESSING_ERROR);
      }
      // The ratio of the two nodes has been specified, so we only
      // need to optimize one node and one weight
      nx = 2;
      opt_problem.set_integer_ratio(integer_ratio);
    }
    */
    Vector integer_ratio;
    if (iorder > 1 && !integer_ratios.empty()) {
      integer_ratio = integer_ratios(iang,range(0,iorder-2));	
      opt_problem.set_integer_ratios(integer_ratio);
      nx = iorder;
    }
    else {
      nx = iorder*2-1;
    }

    if (prior_weight(iang) > 0.0) {
      intVector ordindex = find(prior_orders == iorder);
      if (ordindex.empty()) {
	ERROR << "Quadrature order " << iorder << " not found in " << prior_quadrature;
	THROW(PARAMETER_ERROR);
      }
      int iord = ordindex(0);
      opt_problem.set_prior(prior_weight(iang), prior_mu(iord,range(0,iorder-1)),
			    prior_wt(iord,range(0,iorder-1)));
    }

    Vector x(nx);
    // soft links
    Vector mu;
    Vector wt;

    // Index to last mu in state vector
    int imuend = iorder-1;
    /*
    if (integer_ratio > 1) {
      mu >>= x(range(0,0));
      imuend = 0;
      mu(0) = 1.0/(integer_ratio+1.0);
      wt >>= x(range(1,1));
      wt(0) = 0.5;
    }
    */
    if (iorder > 1 && !integer_ratios.empty()) {
      mu >>= x(range(0,0));
      imuend = 0;
      mu(0) = 1.0/(maxval(integer_ratio)+1.0);
      wt >>= x(range(imuend+1,end));
      wt = 0.5 / (mu(0) + mu(0)*sum(integer_ratio));
    }
    else if (use_init_quad) {
      mu >>= x(range(0,iorder-1));
      // Read initial quadrature from file
      intVector ordindex = find(init_orders == iorder);
      if (ordindex.empty()) {
	ERROR << "Quadrature order " << iorder << " not found in " << init_quadrature;
	THROW(PARAMETER_ERROR);
      }
      int iord = ordindex(0);
      mu = init_mu(iord,range(0,imuend));
      if (iorder > 1) {
	wt >>= x(range(imuend+1,end));
	wt = init_wt(iord,range(0,iorder-2));
      }
    }
    else {
      mu >>= x(range(0,iorder-1));
      // Use evenly spaced quadrature points to start
      if (nx > 1) {
	wt >>= x(range(imuend+1,end));
	wt = 1.0 / iorder;
      }
      if (iorder == 1) {
	mu(0) = 0.5;
      }
      else {
	mu = linspace(0.5/iorder,1-0.5/iorder,iorder);
      }
    }
   
    Vector x_lower, x_upper;
    adept::minimizer_initialize_bounds(nx, x_lower, x_upper);
    Vector mu_lower = x_lower(range(0,imuend));
    Vector mu_upper = x_upper(range(0,imuend));
    mu_lower = 1.0e-4;
    if (!integer_ratio.empty()) {
      mu_upper = 1.0/maxval(integer_ratio);
    }
    else {
      mu_upper = 1.0 - 1.0e-4;
    }
    if (nx > 1) {
      Vector wt_lower = x_lower(range(imuend+1,end));
      Vector wt_upper = x_upper(range(imuend+1,end));
      wt_lower = 1.0e-4;
      wt_upper = 10.0;
    }
    Vector mu2(iorder), wt2(iorder);
    opt_problem.split_state_vector(x, mu2, wt2);
    MinimizerStatus status = minimizer.minimize(opt_problem,
						x, x_lower, x_upper);
    mu.resize(iorder);
    wt.resize(iorder);
    opt_problem.split_state_vector(x, mu, wt);
#pragma omp critical
    {
      std::cout << "*** " << iorder << " ANGLES ***\n";
      std::cout << "  Initial:  mu=" << mu2 << ", wt=" << wt2 << "\n";
      std::cout << "  Solution: mu=" << mu << ", wt=" << wt << "\n";
      std::cout << "  Iterations: " << minimizer.n_iterations()
		<< ", samples: " << minimizer.n_samples()
		<< " (" << minimizer_status_string(status) << ")\n";
      std::cout << "  Cost function reduced from " << sqrt(minimizer.start_cost_function())
		<< " to " << sqrt(minimizer.cost_function()) << "\n";
      std::cout << "\n";
    }

    mu_save(iang,range(0,iorder-1)) = mu;
    wt_save(iang,range(0,iorder-1)) = wt;

    if (save_fluxes) {
      std::stringstream ss;
      ss << "fluxes_" << iorder << ".nc";

      Vector flux_up_toa(ncol), flux_dn_surf(ncol);
      Matrix heating_rate(nlayer,ncol);
      profile.radiative_transfer_lw(mu, wt, flux_dn_surf,
				    flux_up_toa, heating_rate);
      write_fluxes(ss.str(), method_id, flux_up_toa,flux_dn_surf,
		   heating_rate, mu, wt, argc, argv);
    }

  }

  OutputDataFile output(output_file);
  output.define_dimension("order", norder);
  output.define_dimension("node", max_order);

  output.define_variable("order", SHORT, "order");
  output.write_long_name("Order of quadrature", "order");
  
  output.define_variable("mu", DOUBLE, "order", "node");
  output.write_long_name("Cosine of zenith angle", "mu");
  output.write_missing_value(-1.0, "mu");
  
  output.define_variable("weight", DOUBLE, "order", "node");
  output.write_long_name("Weight", "weight");
  output.write_missing_value(0.0, "weight");

  write_standard_attributes(output,
    "Optimized quadrature angles for longwave radiative transfer");
  output.append_history(argc, argv);
  output.write("Hogan, R. J., 2023: What are the optimum discrete angles to use in thermal-infrared radiative transfer calculations? Submitted to Q. J. R. Meteorol. Soc. Available from http://www.met.rdg.ac.uk/~swrhgnrj/publications/discrete_ordinate_angles.pdf", "reference");

  if (!method_id.empty()) {
    output.write(method_id, "quadrature_method_id");
  }
  output.write("Optimization using clear-sky profiles", "quadrature_method");
  output.write(input_file, "training_file");

  std::string config_str;
  config.read(config_str);
  output.write(config_str, "config");
  output.write("The weights in this file are defined such that to compute an irradiance F from a set of radiances I at discrete angles mu, use F=2*pi*sum(weight*mu*I). The weights w presented in Table 1 of Hogan (2023) may be obtained from w=2*mu*weight.","comment");
  
  output.end_define_mode();

  output.write(orders, "order");
  output.write(mu_save, "mu");
  output.write(wt_save, "weight");

  output.close();
}
