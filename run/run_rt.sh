#!/bin/bash

set -e

EXE=$1
shift
TAG=$1
shift
INFILE=$1
shift
QUADRATURE=$1
shift

mkdir -p fluxes

for ORDER in $@
do
  OUTFILE=fluxes/${TAG}_${QUADRATURE}_${ORDER}_fluxes.nc
  $EXE "$INFILE" quadrature_${QUADRATURE}.nc "$ORDER" "$OUTFILE"
done
