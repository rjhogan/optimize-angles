# Optimizing discrete angles for longwave radiative transfer

Author: Robin Hogan <r.j.hogan@ecmwf.int>

This document was last updated 8 March 2023

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

You can run the optimization with

`make check`

See the `Makefile` in the `run` directory for more options.

## LICENCE

(C) Copyright 2019- ECMWF.

This software is licensed under the terms of the Apache Licence Version 2.0
which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

In applying this licence, ECMWF does not waive the privileges and immunities
granted to it by virtue of its status as an intergovernmental organisation
nor does it submit to any jurisdiction.
Copyright statements are given in the file NOTICE.

Note that while the core code in `src/optimize_angles` is owned
entirely by ECMWF, the copyright of several support files in the
`src/include` and `src/tools` directories is either solely or jointly
held with the University of Reading, as stated in the copyright
statements at the top of each file.

## CONTACT

Please email Robin Hogan <r.j.hogan@ecmwf.int> with any queries or bug
fixes, but note that ECMWF does not commit to providing support to
users of this software.
