# Optimizing discrete angles for longwave radiative transfer

Author: Robin Hogan <r.j.hogan@ecmwf.int>

This document was last updated 31 August 2023

This package computes the optimal angles for longwave atmospheric
radiative transfer by minimizing errors in a set of clear-sky
atmospheric profiles. It is an implementation of the algorithm
described in this paper:

Hogan, R. J., 2023: What are the optimum discrete angles to use in
thermal-infrared radiative transfer calculations? Submitted to
Q. J. R. Meteorol. Soc.

The latest version of the paper is available [here](http://www.met.rdg.ac.uk/~swrhgnrj/publications/discrete_ordinate_angles.pdf).

## COMPILING AND RUNNNING

The optimize_angles package uses the autotools build system but two
pre-requisites: the NetCDF C library and the
[Adept](http://www.met.reading.ac.uk/clouds/adept) combined automatic
differentiation, array and optimization library.

If you obtained the software from GitHub, you will need autotools
installed, in which case you can generate the `configure` script with

`autoreconf -i`

Then run

`./configure`
`make`

When running at ECMWF we use the `./configure_ecmwf.sh` script which
calls `./configure` with extra arguments - this may be useful when
compiling the code on your own system.

You can run the optimization with

`make check`

See the `Makefile` in the `run` directory for more options.

## LICENCE

Files in the `src/optimize_angles` directory (the core algorithm):
Copyright (C) 2022- ECMWF.

Files in the `src/tools` and `src/include` directories: Copyright (C)
2015- ECMWF, Copyright (C) 2006-2015 University of Reading

See also the copyright statements at the top of each source file.

This software is licensed under the terms of the Apache Licence Version 2.0
which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

In applying this licence, ECMWF does not waive the privileges and immunities
granted to it by virtue of its status as an intergovernmental organisation
nor does it submit to any jurisdiction.
Copyright statements are given in the file NOTICE.

## CONTACT

Please email Robin Hogan <r.j.hogan@ecmwf.int> with any queries or bug
fixes, but note that ECMWF does not commit to providing support to
users of this software.
