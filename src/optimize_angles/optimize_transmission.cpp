#include <iostream>

#include <adept_arrays.h>
#include <adept_optimize.h>

using namespace adept;

class Transmission : public Optimizable {
public:
  Transmission(int nangs) : nangles(nangs) {
    initialize();
  }

  void initialize() {
    int nod = 351;
    //int nod = 13;
    int nmu = 1000;
    optical_depth = pow(10.0, linspace(-5.0,2.0,nod));
    transmittance.resize(nod);
    Vector mu = linspace(0.5/nmu,1.0-0.5/nmu,nmu);
    Real recip_sum_mu = 1.0 / sum(mu);
    for (int iod = 0; iod < optical_depth.size(); ++iod) {
      transmittance(iod) = sum(mu*exp(-optical_depth(iod)/mu))*recip_sum_mu;
    }
    weight.resize(nod);
    weight(range(1,end-1))
      = 0.5 * (transmittance(range(0,end-2))-transmittance(range(2,end)));
    weight(0)   = 1.0-0.5*(transmittance(0)+transmittance(1));
    weight(end) = 0.5*(transmittance(end-1)+transmittance(end));
    if (iverbose > 1) {
      std::cout << "O: " << optical_depth << "\n";
      std::cout << "T: " << transmittance << "\n";
      std::cout << "W: " << weight << "\n";
    }
	
  }
    
  virtual bool provides_derivative(int order) {
    if (order >= 0 && order <= 2) {
      return true;
    }
    else {
      return false;
    }
  }

  void set_verbose(int iverb) { iverbose = iverb; }

  // Split state vector into cosine of zenith angles, mu, and weights,
  // wt. The final value of wt must be computed for energy
  // conservation given that sum(wt*mu)=0.5.
  template <bool IsActive>
  void split_state_vector(const Array<1,Real,IsActive>& x,
			  Array<1,Real,IsActive>& mu,
			  Array<1,Real,IsActive>& wt) {
    mu = x(range(0,nangles-1));
    if (nangles > 1) {
      wt(range(0,nangles-2)) = x(range(nangles,end));
      wt(end) = (0.5-sum(wt(range(0,nangles-2))*mu(range(0,nangles-2)))) / mu(end);
    }
    else {
      wt(0) = 0.5 / mu(0);
    }
  }
  
  // Calculate transmittance for all optical depths
  template <bool IsActive>
  Array<1,Real,IsActive> calc_y(const Array<1,Real,IsActive>& mu,
				const Array<1,Real,IsActive>& wt) {
    Array<1,Real,IsActive> y(optical_depth.size());
    y = 0.0;
    for (int imu = 0; imu < nangles; ++imu) {
      y = y + 2.0*mu(imu)*wt(imu)*exp(-optical_depth/mu(imu));
    }
    return y;
  }

  template <bool IsActive>
  typename internal::active_scalar<Real,IsActive>::type calc_cost_function_active(const Array<1,Real,IsActive>& x) {
    Array<1,Real,IsActive> mu(nangles), wt(nangles);
    split_state_vector(x, mu, wt);
    Array<1,Real,IsActive> y = calc_y(mu, wt);
    Array<1,Real,IsActive> dy = y - transmittance;
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
    cost.set_gradient(1.0);
    stack.reverse();
    gradient = xactive.get_gradient();
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
    aVector dy = y-transmittance;
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
  Vector transmittance;
  Vector optical_depth;
  Vector weight;
  int nangles;
  int iverbose = 1;
};

int main()
{
  
  if (!adept::have_linear_algebra()) {
    std::cout << "Adept compiled without linear-algebra support: minimizer not available\n";
    return 0;
  }
  
  adept::set_array_print_style(PRINT_STYLE_MATLAB);
  
  for (int iangs = 1; iangs <= 8; ++iangs) {
    std::cout << "*** ANGLES: " << iangs << " ***\n";
    Transmission transmission_problem(iangs);
    //Minimizer minimizer(MINIMIZER_ALGORITHM_LEVENBERG_MARQUARDT);
    Minimizer minimizer(MINIMIZER_ALGORITHM_LIMITED_MEMORY_BFGS);
    minimizer.set_converged_gradient_norm(0.00000000000001);
    minimizer.set_levenberg_damping_start(0.25);
    minimizer.ensure_updated_state(2);
    minimizer.set_max_iterations(1000000);

    int nx = iangs*2-1;
    Vector x(nx);
    // soft links
    Vector mu = x(range(0,iangs-1));
    Vector wt;
    if (nx > 1) {
      wt >>= x(range(iangs,end));
      wt = 1.0 / iangs;
    }
    if (iangs == 1) {
      mu(0) = 0.5;
    }
    else {
      mu = linspace(0.5/iangs,1-0.5/iangs,iangs);
    }

    Vector x_lower, x_upper;
    adept::minimizer_initialize_bounds(nx, x_lower, x_upper);
    Vector mu_lower = x_lower(range(0,iangs-1));
    Vector mu_upper = x_upper(range(0,iangs-1));
    mu_lower = 1.0e-4;
    mu_upper = 1.0 - 1.0e-4;
    if (nx > 1) {
      Vector wt_lower = x_lower(range(iangs,end));
      Vector wt_upper = x_upper(range(iangs,end));
      wt_lower = 1.0e-4;
      wt_upper = 10.0;
    }
    MinimizerStatus status = minimizer.minimize(transmission_problem,
						x, x_lower, x_upper);
    std::cout << "  Status: " << minimizer_status_string(status) << "\n";
    wt.resize(iangs);
    transmission_problem.split_state_vector(x, mu, wt);
    std::cout << "  Solution: mu=" << mu << "\n";
    std::cout << "  Solution: wt=" << wt << "\n";
    std::cout << "  Number of iterations: " << minimizer.n_iterations() << "\n";
    std::cout << "  Number of samples: " << minimizer.n_samples() << "\n";
    std::cout << "  RMSE: " << sqrt(minimizer.cost_function()) << "\n";
    //std::cout << "  T: " << transmission_problem.calc_y(mu,wt) << "\n";
    std::cout << "\n";
    std::cerr << "      nodes(1:" << iangs << ") = " << mu << "\n";
    std::cerr << "      weights(1:" << iangs << ") = " << wt << "\n";
  }
}
