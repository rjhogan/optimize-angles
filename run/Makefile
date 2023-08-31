# This Makefile enables the user to recreate the quadratures presented
# in the Hogan (2024) paper. You should have already compiled the code
# as described in the README in the top-level directory.

# Executables
OPT_TRANSMISSION = ../src/optimize_angles/optimize_transmission
OPT_ANGLES = ../src/optimize_angles/optimize_angles
RT = ../src/optimize_angles/radiative_transfer_lw

# 50 CKDMIP clear-sky Evaluation-1 profiles, used for training
INFILE1 = optical-depth/ecckd-1.0_evaluation1_lw_climate_fsck-32b_optical-depth_present.nc

# 50 CKDMIP clear-sky Evaluation-2 profiles, used for independent evaluation
INFILE2 = optical-depth/ecckd-1.0_evaluation2_lw_climate_fsck-32b_optical-depth_present.nc

# Weighting between fluxes and heating rates in the optimization
FLUXWEIGHT=0.02

# Type "make TAG=xxx ..." to use a different tag for output files
TAG = default
OPTTAG = optimized-$(TAG)

OUTBASE = quadrature/quadrature_$(OPTTAG)
OUTFILE = $(OUTBASE).nc

# The optimization is trained against some "truth", calculated using
# an accurate quadrature scheme with the largest number of angles
# available - in this case we use Gauss-Jacobi-5 quadrature with 64
# angles (equivalent to 128 streams).
REF=quadrature/quadrature_gauss-jacobi-5.nc

# Numbers of angles to compute quadratures for
ANGLES_OPT="1 2 3 4 5 6 7 8 10 12 14 16"
ANGLES_SUBSET="1 2 3 4"
ANGLES_RATIO="2 3 4"

# Numbers of angles to use for radiative transfer with Gaussian
# quadrature schemes
ANGLES_GAUSS="1 2 3 4 5 6 7 8 10 12 14 16 32 64"

# "make" and "make all" just print the following message, in order
# that "make" in the top-level directory only compiles the code.
all:
	@echo '***********************************************************************'
	@echo '  Type "make check" to produce quadratures with up to four angles      '
	@echo '  Type "make full" in the "run" directory for many more angles (slow)  '
	@echo '***********************************************************************'

# Do all optimization used in the paper (slow)
full: angles_full rt_full

# Do only a limited number of angles that can be compared to the
# REFERENCE data in the quadrature directory
check: angles_subset rt_opt_subset
	@echo '***********************************************************************'
	@echo '  Now compare run/quadrature/quadrature_optimized-default-subset.nc to '
	@echo '    run/quadrature/REFERENCE-quadrature_optimized-default-subset.nc    '
	@echo '***********************************************************************'

# The Optimized, Optimized-IR and Optimized-IRJP quadratures in the
# paper
angles_full: angles_opt angles_ratio angles_ratio_jacobi_prior

# Reproduce the numbers in Table 1 of the paper (only up to four
# angles)
table: angles_subset angles_ratio angles_ratio_jacobi_prior

# Compute the Optimized quadrature for all angles used in the paper
angles_opt: $(OPT_ANGLES)
	$(OPT_ANGLES) flux_weight=$(FLUXWEIGHT) broadband_weight=0 \
		orders=$(ANGLES_OPT) \
		ref_quadrature=$(REF) init_quadrature=$(REF) \
		$(INFILE1) $(OUTFILE)

# Compute the Optimized quadrature for the angles shown in Table 1
angles_subset: $(OPT_ANGLES)
	$(OPT_ANGLES) flux_weight=$(FLUXWEIGHT) broadband_weight=0 \
		orders=$(ANGLES_SUBSET) \
		ref_quadrature=$(REF) init_quadrature=$(REF) \
		$(INFILE1) $(OUTBASE)-subset.nc

# Compute the "Optimized-IR" quadrature: the angles are in specified
# integer ratios of each other, in order for a downstream application
# to reduce the number of exponentials
angles_ratio: $(OPT_ANGLES)
	$(OPT_ANGLES) flux_weight=$(FLUXWEIGHT) broadband_weight=0 \
		orders="2 3 4" integer_ratios[3][3]="4 0 0 \
	5 12 0 \
	5 16 32" \
		$(INFILE1) $(OUTBASE)-ratio.nc

# Compute the "Optimized-IRJP" quadrature: the angles are in specified
# integer ratios of each other, and are constrained by the
# Gauss-Jacob-5 quadrature, producing quadratures appropriate for
# all-sky (scattering) atmospheres that can also be used in downstream
# applications with a reduced number of exponentials
angles_ratio_jacobi_prior: $(OPT_ANGLES)
	$(OPT_ANGLES) flux_weight=$(FLUXWEIGHT) broadband_weight=0 \
		orders="2 3 4" integer_ratios[3][3]="3 0 0 \
	4 8 0 \
	5 13 20" \
		ref_quadrature=$(REF) \
		prior_quadrature=$(REF) prior_weight="0.001" \
		$(INFILE1) $(OUTBASE)-jacobi-5-ratio-prior.nc

# Compute quadratures using only the transmission curve in Fig. 1a of
# the paper (not well maintained)
transmission: $(OPT_TRANSMISSION)
	$(OPT_TRANSMISSION) > opt.txt

# Perform all radiative transfer
rt_full: rt_gauss rt_opt rt_opt_ratio rt_opt_ratio_jacobi_prior rt_lacis rt_elsasser

# Run clear-sky radiative transfer for the Evaluation-1 and
# Evaluation-2 datasets using the three Gaussian quadrature schemes
rt_gauss:
	mkdir -p rt_out
	./run_rt.sh $(RT) evaluation1 $(INFILE1) gauss-legendre $(ANGLES_GAUSS)
	./run_rt.sh $(RT) evaluation1 $(INFILE1) gauss-laguerre $(ANGLES_GAUSS)
	./run_rt.sh $(RT) evaluation1 $(INFILE1) gauss-jacobi-5 $(ANGLES_GAUSS)
	./run_rt.sh $(RT) evaluation2 $(INFILE2) gauss-legendre $(ANGLES_GAUSS)
	./run_rt.sh $(RT) evaluation2 $(INFILE2) gauss-laguerre $(ANGLES_GAUSS)
	./run_rt.sh $(RT) evaluation2 $(INFILE2) gauss-jacobi-5 $(ANGLES_GAUSS)

# Run clear-sky radiative transfer using Optimized quadrature with a
# subset of angles shown in Table 1 of the paper
rt_opt_subset:
	mkdir -p rt_out
	./run_rt.sh $(RT) evaluation1 $(INFILE1) $(OPTTAG)-subset $(ANGLES_SUBSET)
	./run_rt.sh $(RT) evaluation2 $(INFILE2) $(OPTTAG)-subset $(ANGLES_SUBSET)

# Run clear-sky radiative transfer using Optimized quadrature for all
# angles
rt_opt:
	mkdir -p rt_out
	./run_rt.sh $(RT) evaluation1 $(INFILE1) $(OPTTAG) $(ANGLES_OPT)
	./run_rt.sh $(RT) evaluation2 $(INFILE2) $(OPTTAG) $(ANGLES_OPT)

# Run clear-sky radiative transfer using the Optimized-IR or
# Optimized-IRJP quadratures
rt_opt_ratio:
	mkdir -p rt_out
	./run_rt.sh $(RT) evaluation1 $(INFILE1) $(OPTTAG)-ratio $(ANGLES_RATIO)
	./run_rt.sh $(RT) evaluation2 $(INFILE2) $(OPTTAG)-ratio $(ANGLES_RATIO)
rt_opt_ratio_jacobi_prior:
	mkdir -p rt_out
	./run_rt.sh $(RT) evaluation1 $(INFILE1) $(OPTTAG)-jacobi-5-ratio-prior $(ANGLES_RATIO)
	./run_rt.sh $(RT) evaluation2 $(INFILE2) $(OPTTAG)-jacobi-5-ratio-prior $(ANGLES_RATIO)

# Run clear-sky radiative transfer using the Elsasser or Lacis
# quadratures (only for 1 or 3 angles, respectively)
rt_elsasser:
	mkdir -p rt_out
	./run_rt.sh $(RT) evaluation1 $(INFILE1) elsasser 1
	./run_rt.sh $(RT) evaluation2 $(INFILE2) elsasser 1
rt_lacis:
	mkdir -p rt_out
	./run_rt.sh $(RT) evaluation1 $(INFILE1) lacis 3
	./run_rt.sh $(RT) evaluation2 $(INFILE2) lacis 3

# "make clean" does not delete any data produced by the user - need to
# run "make dist-clean"
clean:

# Clear fluxes and quadratures produced by this Makefile
dist-clean: clean-fluxes clean-quadratures

clean-fluxes:
	rm -f fluxes/evaluation?_*_fluxes.nc

clean-quadratures:
	rm -f quadrature/quadrature_optimized*.nc
