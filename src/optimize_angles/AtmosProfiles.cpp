#include "DataFile.h"
#include "AtmosProfiles.h"

static const adept::Real ACCEL_GRAVITY     = 9.80665; // m s-2
static const adept::Real SPECIFIC_HEAT_AIR = 1004.0;  // J kg-1 K-1

template <bool IsActive>
static void
calc_heating_rate(const Matrix& pressure_hl, //< Half-level pressure (Pa)
		  const Array<2,Real,IsActive>& flux_dn,//< Flux down (W m-2)
		  const Array<2,Real,IsActive>& flux_up,//< Flux up (W m-2)
		  Array<2,Real,IsActive>& heating_rate) {//< Heating rate (K d-1)
  // Factor to convert from difference in net flux across a layer to
  // heating rate in K/day
  Matrix conversion = -(86400.0*ACCEL_GRAVITY/SPECIFIC_HEAT_AIR)
    / transpose(pressure_hl(__,range(1,end))-pressure_hl(__,range(0,end-1)));

  heating_rate = conversion
    * (  flux_dn(range(1,end),__) - flux_dn(range(0,end-1),__)
       - flux_up(range(1,end),__) + flux_up(range(0,end-1),__)  );

}

Matrix
AtmosProfiles::heating_rate_weight() const {
  return transpose((sqrt(pressure_hl(__,range(1,end))) - sqrt(pressure_hl(__,range(0,end-1))))
		   / spread<1>(sqrt(pressure_hl(__,end)),nlayer()));
}
  
void
AtmosProfiles::read(const std::string& file_name) {
  DataFile file(file_name);
  file.read(pressure_hl, "pressure_hl");
  file.read(planck_hl, "planck_hl");
  if (file.exist("optical_depth")) {
    file.read(optical_depth, "optical_depth");
  }
  else {
    file.read(optical_depth, "od_lw");
  }
  if (file.exist("planck_surf")) {
    file.read(planck_surf,"planck_surf");
  }
  else {
    Matrix lw_emissivity, lw_emission;
    file.read(lw_emissivity, "lw_emissivity");
    file.read(lw_emission,   "lw_emission");
    planck_surf = lw_emission / lw_emissivity;
  }
  file.close();
}

template <bool IsActive>
void
AtmosProfiles::radiative_transfer_lw(const Array<1,Real,IsActive>& mu,
				     const Array<1,Real,IsActive>& weight,
				     Array<1,Real,IsActive> flux_surf,
				     Array<1,Real,IsActive> flux_toa,
				     Array<2,Real,IsActive> heating_rate) const {
  typedef typename internal::active_scalar<Real,IsActive>::type areal;
  int nmu = mu.size();
  int ncol  = optical_depth.size(0);
  int nlay  = optical_depth.size(1);
  int nspec = optical_depth.size(2);
  
  flux_surf = 0.0;
  flux_toa = 0.0;
  heating_rate = 0.0;
  Array<2,Real,IsActive> flux_dn(nlay+1,nspec), flux_up(nlay+1,nspec);
  Array<2,Real,IsActive> emissivity(nlay,nspec), factor(nlay,nspec);
  Array<2,Real,IsActive> flux_dn_col(nlay+1,ncol), flux_up_col(nlay+1,ncol);
  flux_dn_col = 0.0;
  flux_up_col = 0.0;
  for (int icol = 0; icol < ncol; ++icol) {
    for (int imu = 0; imu < nmu; ++imu) {
      areal secant = 1.0/mu(imu);
      emissivity = 1.0 - exp(-secant*optical_depth(icol,__,__));
      factor.where(emissivity > 1.0e-5)
	= either_or(1.0 - mu(imu)*(emissivity/optical_depth(icol,__,__)),
		    0.5 * emissivity);
      // Work down from top
      flux_dn(0,__) = 0.0;
      for (int ilay = 0; ilay < nlay; ++ilay) {
	flux_dn(ilay+1,__) = flux_dn(ilay,__) * (1.0 - emissivity(ilay,__))
	  + planck_hl(icol,ilay,__)   * (emissivity(ilay,__)-factor(ilay,__))
	  + planck_hl(icol,ilay+1,__) * factor(ilay,__);
      }
      // Work up from surface
      flux_up(nlay,__) = planck_surf(icol,__);
      for (int ilay = nlay-1; ilay >= 0; --ilay) {
	flux_up(ilay,__) = flux_up(ilay+1,__) * (1.0 - emissivity(ilay,__))
	  + planck_hl(icol,ilay+1,__) * (emissivity(ilay,__)-factor(ilay,__))
	  + planck_hl(icol,ilay,__)   * factor(ilay,__);
      }
      flux_up *= 2.0*weight(imu)*mu(imu);
      flux_dn *= 2.0*weight(imu)*mu(imu);
      flux_up_col(__,icol) += sum(flux_up,1);
      flux_dn_col(__,icol) += sum(flux_dn,1);
      flux_surf(icol) += sum(flux_dn(end,__));
      flux_toa (icol) += sum(flux_up(0,__));
    }
  }
  calc_heating_rate(pressure_hl, flux_dn_col, flux_up_col, heating_rate);
}

// Explicit instantiations
template
void
AtmosProfiles::radiative_transfer_lw<false>(const Array<1,Real,false>& mu,
					    const Array<1,Real,false>& weight,
					    Array<1,Real,false> flux_surf,
					    Array<1,Real,false> flux_toa,
					    Array<2,Real,false> heating_rate) const;
template
void
AtmosProfiles::radiative_transfer_lw<true>(const Array<1,Real,true>& mu,
					   const Array<1,Real,true>& weight,
					   Array<1,Real,true> flux_surf,
					   Array<1,Real,true> flux_toa,
					   Array<2,Real,true> heating_rate) const;
