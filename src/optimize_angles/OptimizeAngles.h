// OptimizeAngles.h - Optimizing class inheriting from Adept's Optimizable -*- C++ -*-
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

#include <adept_arrays.h>
#include <adept_optimize.h>

#include "AtmosProfile.h"

using namespace adept;

template <class Profile>
class OptimizeAngles : public Optimizable {
public:
  OptimizeAngles(const Profile& prof, int nangs, Vector yref, Vector wt)
    : profile(prof), y_ref(yref), weight(wt), nangles(nangs) { }
   
  virtual bool provides_derivative(int order) {
    if (order >= 0 && order <= 2) {
      return true;
    }
    else {
      return false;
    }
  }

  void set_verbose(int iverb) { iverbose = iverb; }

  void set_integer_ratios(intVector ratios) {
    integer_ratios = ratios;
    if (any(integer_ratios <= 1)) {
      ERROR << "User-specified integer ratios must all be more than 1";
      THROW(PARAMETER_ERROR);
    }
    LOG << "The second and subsequent angles are the following multiples of the first: "
	<< integer_ratios << "\n";
  }

  void set_prior(Real prior_weight_, Vector prior_mu_, Vector prior_wt_) {
    prior_weight = prior_weight_;
    prior_mu     = prior_mu_;
    prior_wt     = prior_wt_;
  }
  
  // Split state vector into cosine of zenith angles, mu, and weights,
  // wt. The final value of wt must be computed for energy
  // conservation given that sum(wt*mu)=0.5.
  template <bool IsActive>
  void split_state_vector(const Array<1,Real,IsActive>& x,
			  Array<1,Real,IsActive>& mu,
			  Array<1,Real,IsActive>& wt) {
    mu = x(range(0,nangles-1));
    /*
    if (nangles == 2 && integer_ratio > 1) {
      mu(1) = mu(0)*integer_ratio;
      wt(0) = x(1);
      wt(1) = (0.5-wt(0)*mu(0))/mu(1);
    }
    */
    if (nangles > 1 && !integer_ratios.empty()) {
      // Only the first angle is in the state vector - the others are
      // scaled from the first
      //mu(range(1,nangles-1)) = min(mu(0) * integer_ratios, 1.0);
      mu(range(1,nangles-1)) = integer_ratios;
      //mu(1) = 3.7529;
      mu(range(1,nangles-1)) *= mu(0);
      wt(range(0,nangles-2)) = x(range(1,end));
      wt(end) = (0.5-sum(wt(range(0,nangles-2))*mu(range(0,nangles-2)))) / mu(end);
    }
    else if (nangles > 1) {
      wt(range(0,nangles-2)) = x(range(nangles,end));
      wt(end) = (0.5-sum(wt(range(0,nangles-2))*mu(range(0,nangles-2)))) / mu(end);
    }
    else {
      wt(0) = 0.5 / mu(0);
    }
  }
  
  // Calculate spectral fluxes and heating rates for all optical depths
  template <bool IsActive>
  Array<1,Real,IsActive> calc_y(const Array<1,Real,IsActive>& mu,
				const Array<1,Real,IsActive>& wt) {
    int nlayer = profile.nlayer();
    int ncol   = profile.ncolumn();
    int ny     = (nlayer+2)*ncol+2;
    Array<1,Real,IsActive> y(ny);
    Array<2,Real,IsActive> heating_rate = y(range(0,ncol*nlayer-1)).reshape(nlayer,ncol);
    Array<1,Real,IsActive> flux_surf    = y(range(ncol*nlayer,ncol*(nlayer+1)-1));
    Array<1,Real,IsActive> flux_toa     = y(range(ncol*(nlayer+1),ncol*(nlayer+2)-1));

    profile.radiative_transfer_lw(mu, wt, flux_surf, flux_toa, heating_rate);
    y(end-1) = sum(flux_surf);
    y(end)   = sum(flux_toa);
    return y;
  }

  template <bool IsActive>
  typename internal::active_scalar<Real,IsActive>::type calc_cost_function_active(const Array<1,Real,IsActive>& x) {
    Array<1,Real,IsActive> mu(nangles), wt(nangles);
    split_state_vector(x, mu, wt);
    typename internal::active_scalar<Real,IsActive>::type cost = 0.0;
    // If prior weight is very large then we skip the (slow) profile
    // radiative transfer calculations entirely and only use the prior
    // - this will obviously return the prior unless we are using
    // integer ratios
    if (prior_weight < 1.0e10) {
      Array<1,Real,IsActive> y = calc_y(mu, wt);
      Array<1,Real,IsActive> dy = y - y_ref;
      //    LOG << "mu=" << mu << "\n";
      //    LOG << "yref=" << y_ref << "\n";
      //    LOG << "y=" << y << "\n";
      cost += 0.5*sum(dy*dy*weight);
    }
    if (prior_weight > 0.0) {
      cost += prior_weight * sum((mu-prior_mu)*(mu-prior_mu)
				      +(wt-prior_wt)*(wt-prior_wt));
    }
    return cost;
  }

  virtual Real calc_cost_function(const Vector& x) {
    return calc_cost_function_active(x);
  }

  virtual Real calc_cost_function_gradient(const Vector& x,
					   Vector gradient) {
    Stack stack;
    aVector xactive = x;
    stack.new_recording();
    aReal cost = calc_cost_function_active(xactive);
    //stack.print_status();
    cost.set_gradient(1.0);
    stack.reverse();
    gradient = xactive.get_gradient();
    //    LOG << cost << " " << x << " " << gradient << "\n";
    return value(cost);
  }

  virtual Real calc_cost_function_gradient_hessian(const Vector& x,
						   Vector gradient,
						   SymmMatrix& hessian) {
    Stack stack;
    aVector xactive = x;
    stack.new_recording();
    aVector mu(nangles), wt(nangles);
    split_state_vector(xactive, mu, wt);
    aVector y = calc_y(mu,wt);
    aVector dy = y-y_ref;
    aReal cost = 0.5*sum(dy*dy*weight);
    stack.independent(xactive);
    stack.dependent(y);
    Matrix jac = stack.jacobian();
    hessian  = jac.T() ** diag_matrix(weight) ** jac;
    gradient = jac.T() ** (weight*value(dy));
    return value(cost);
  }

private:
  // Data
  const Profile& profile;
  Vector y_ref;
  Vector weight;
  int nangles;
  int iverbose = 1;
  intVector integer_ratios;

  // Prior values of nodes and corresponding weights
  Vector prior_mu, prior_wt;
  // Weight of prior term in cost function
  Real prior_weight = 0.0;
};
