# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])

# Remember to update the version number for a new release
AC_INIT([optimize_angles], [1.0], [r.j.hogan@ecmwf.int], [optimize_angles])

# Most code is C++
AC_LANG([C++])

AC_CONFIG_SRCDIR([src/include/readconfig.h])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([-Wall -Wno-error foreign])
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs
AC_PROG_CXX
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
AC_PROG_LIBTOOL

# Check for OpenMP
AC_OPENMP
AC_SUBST(AM_CXXFLAGS,"$OPENMP_CXXFLAGS")
AC_SUBST(AM_CFLAGS,"$OPENMP_CXXFLAGS")

# Header directories needed during the compilation; also need to run
# multiple Adept stacks simultaneously requiring the thread-safe
# option
AC_SUBST(AM_CPPFLAGS,"-I../include -DADEPT_STORAGE_THREAD_SAFE")

# Check for Adept library, using m4/adept.m4
AX_CHECK_ADEPT([have_adept=yes], [have_adept=no])

# Check for NetCDF library, using m4/netcdf.m4
have_netcdf=yes
AX_CHECK_NETCDF([have_netcdf=yes], [have_netcdf=no])
AM_CONDITIONAL([HAVE_NETCDF], [test "$have_netcdf" = yes])
AM_COND_IF([HAVE_NETCDF], AC_SUBST([NETCDF_DIR], [$netcdf_prefix]))

# Makefiles to create
AC_CONFIG_FILES([Makefile
                 src/tools/Makefile
                 src/optimize_angles/Makefile])

# Provide a summary to the user
AC_MSG_NOTICE([********************* Summary **************************************])
AC_MSG_NOTICE([  CXX      = $CXX ])
AC_MSG_NOTICE([  CPPFLAGS = $CPPFLAGS])
AC_MSG_NOTICE([  CXXFLAGS = $CXXFLAGS $OPENMP_CXXFLAGS])
AC_MSG_NOTICE([  LDFLAGS  =  $LDFLAGS])
AC_MSG_NOTICE([  LIBS     = $LIBS])
AC_MSG_NOTICE([Typing "make; make install" will install the code and auxiliary files])
AC_MSG_NOTICE([********************************************************************])

AC_OUTPUT
