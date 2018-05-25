#!/bin/bash
CPPFLAGS="-traditional-cpp -DLINUX -DNO_SHR_VMATH -DMPI_P"
NETCDFINC="  -g -I/opt/netCDF_old/include "
NETCDFLIB="  -L/opt/netCDF_old/lib -lnetcdff -lnetcdf "
MPIINC="  -I/opt/intel_old/impi/3.2.2.006/include64 "
MPILIB="  -L/opt/intel_old/impi/3.2.2.006/lib64 "

export CC=mpiicc
export CXX=mpiicpc
export FC=mpiifort
#export CPP=/usr/bin/cpp
export CFLAGS="-O2 -DFORTRANUNDERSCORE -g"
export CXXFLAGS="-O2 -c -DFORTRANUNDERSCORE -g"
export FFLAGS="-g -free -O2 -c -i4  -r8 -convert big_endian -assume byterecl -fp-model precise"
export INCLDIR=" ${NETCDFINC} ${MPIINC} "
export SLIBS=" ${NETCDFLIB} ${MPILIB} "
export CPPFLAGS="${CPPFLAGS} ${INCLDIR} "

make -j 8
