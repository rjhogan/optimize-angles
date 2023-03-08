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

  void set_integer_ratio(int ratio) { integer_ratio = ratio; }

  // Split state vector into cosine of zenith angles, mu, and weights,
  // wt. The final value of wt must be computed for energy
  // conservation given that sum(wt*mu)=0.5.
  template <bool IsActive>
  void split_state_vector(const Array<1,Real,IsActive>& x,
			  Array<1,Real,IsActive>& mu,
			  Array<1,Real,IsActive>& wt) {
    mu = x(range(0,nangles-1));
    if (nangles == 2 && integer_ratio > 1) {
      mu(1) = mu(0)*integer_ratio;
      wt(0) = x(1);
      wt(1) = (0.5-wt(0)*mu(0))/mu(1);
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
    Array<1,Real,IsActive> y = calc_y(mu, wt);
    Array<1,Real,IsActive> dy = y - y_ref;
    //    LOG << "mu=" << mu << "\n";
    //    LOG << "yref=" << y_ref << "\n";
    //    LOG << "y=" << y << "\n";
    typename internal::active_scalar<Real,IsActive>::type cost = 0.5*sum(dy*dy*weight);
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
  int integer_ratio = -1;
};
