bin_PROGRAMS = optimize_transmission optimize_angles radiative_transfer_lw

optimize_transmission_SOURCES = optimize_transmission.cpp

optimize_angles_SOURCES = optimize_angles.cpp \
	AtmosProfile.cpp AtmosProfile.h \
	AtmosProfiles.cpp AtmosProfiles.h \
	write_fluxes.cpp write_fluxes.h

radiative_transfer_lw_SOURCES = radiative_transfer_lw.cpp \
	AtmosProfiles.cpp AtmosProfiles.h \
	write_fluxes.cpp write_profiles.h

LDADD = ../tools/libtools.a
