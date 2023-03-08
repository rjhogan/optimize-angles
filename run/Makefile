OPT_TRANSMISSION = ../src/optimize_angles/optimize_transmission
OPT_ANGLES = ../src/optimize_angles/optimize_angles
RT = ../src/optimize_angles/radiative_transfer_lw

INFILE1 = optical-depth/ecckd-1.0_evaluation1_lw_climate_fsck-32b_optical-depth_present.nc
INFILE2 = optical-depth/ecckd-1.0_evaluation2_lw_climate_fsck-32b_optical-depth_present.nc

FLUXWEIGHT=0.02

# Type "make TAG=xxx check" to use a different tag for output files
TAG = default

OUTBASE = quadrature/quadrature_optimized-$(TAG)
OUTFILE = $(OUTBASE).nc
LOGFILE = $(OUTBASE).log

OUTFILEx4 = $(OUTBASE)-ratio4.nc
OUTFILEx3 = $(OUTBASE)-ratio3.nc

REF=quadrature/quadrature_gauss-jacobi-5.nc
all:
		"orders=1 2 3 4 5 6 7 8 10 12 14 16" \

check: run
run: run_angles
#ref_quadrature=$(REF) 
run_angles: $(OPT_ANGLES)
	$(OPT_ANGLES) flux_weight=$(FLUXWEIGHT) broadband_weight=0 \
		"orders=1 2 3 4 5 6 7 8 10 12 14 16" \
		ref_quadrature=$(REF) init_quadrature=$(REF) \
		$(INFILE1) $(OUTFILE)

run_angles_x4: $(OPT_ANGLES)
	debug $(OPT_ANGLES) flux_weight=$(FLUXWEIGHT) broadband_weight=0 \
		"orders=2" integer_ratio=4 \
		ref_quadrature=$(REF) init_quadrature=$(REF) \
		$(INFILE1) $(OUTFILEx4)
run_angles_x3: $(OPT_ANGLES)
	debug $(OPT_ANGLES) flux_weight=$(FLUXWEIGHT) broadband_weight=0 \
		"orders=2" integer_ratio=3 \
		ref_quadrature=$(REF) init_quadrature=$(REF) \
		$(INFILE1) $(OUTFILEx3)

run_transmission: $(OPT_TRANSMISSION)
	$(OPT_TRANSMISSION) > opt.txt

run_rt: run_rt_gauss run_rt_opt

ANGLES="1 2 3 4 5 6 7 8 10 12 14 16 32 64"
run_rt_gauss:
	mkdir -p rt_out
	./run_rt.sh $(RT) evaluation1 $(INFILE1) gauss-legendre $(ANGLES)
	./run_rt.sh $(RT) evaluation1 $(INFILE1) gauss-laguerre $(ANGLES)
	./run_rt.sh $(RT) evaluation1 $(INFILE1) gauss-jacobi-5 $(ANGLES)
	./run_rt.sh $(RT) evaluation2 $(INFILE2) gauss-legendre $(ANGLES)
	./run_rt.sh $(RT) evaluation2 $(INFILE2) gauss-laguerre $(ANGLES)
	./run_rt.sh $(RT) evaluation2 $(INFILE2) gauss-jacobi-5 $(ANGLES)

ANGLES_OPT="1 2 3 4 5 6 7 8 10 12 14 16"
run_rt_opt:
	./run_rt.sh $(RT) evaluation1 $(INFILE1) $(TAG) $(ANGLES_OPT)
	./run_rt.sh $(RT) evaluation2 $(INFILE2) $(TAG) $(ANGLES_OPT)
	./run_rt.sh $(RT) evaluation1 $(INFILE1) $(TAG)-ratio4 2
	./run_rt.sh $(RT) evaluation2 $(INFILE2) $(TAG)-ratio4 2
	./run_rt.sh $(RT) evaluation1 $(INFILE1) $(TAG)-ratio3 2
	./run_rt.sh $(RT) evaluation2 $(INFILE2) $(TAG)-ratio3 2

run_rt_elsasser:
	./run_rt.sh $(RT) evaluation1 $(INFILE1) elsasser 1
	./run_rt.sh $(RT) evaluation2 $(INFILE2) elsasser 1

run_rt_lacis:
	./run_rt.sh $(RT) evaluation1 $(INFILE1) lacis 3
	./run_rt.sh $(RT) evaluation2 $(INFILE2) lacis 3

clean: clean_fluxes clean_quadratures

clean_fluxes:
	rm -f fluxes/evaluation?_*_fluxes.nc

clean_quadratures:
	rm -f quadrature/quadrature_optimized*.nc
