#!/bin/bash
#
set -x
  export machine=${1:-wcoss_dell_p3}

. $MODULESHOME/init/sh
module purge

#
#  Make sure the directory "CMakeModules/Modules" exists!
#  if not do "git submodule update --init --recursive"
#
 if [ $machine = wcoss_dell_p3 ] ; then
#  module load ips/18.0.1.163
#  module load impi/18.0.1
#  module load lsf/10.1
#  module load python/3.6.3
#  module load cmake/3.20.2

   module use /usrx/local/nceplibs/dev/hpc-stack/libs/hpc-stack/modulefiles/stack

   module load hpc/1.1.0
   module load hpc-ips/18.0.1.163
   module load hpc-impi/18.0.1

   module load cmake/3.16.2
   module load jasper/2.0.22
   module load zlib/1.2.11
   module load png/1.6.35
   module load hdf5/1.10.6
   module load netcdf/4.7.4
   module load bacio/2.4.1
   module load crtm/2.3.0
   module load g2/3.4.1
   module load g2tmpl/1.10.0
   module load ip/3.3.3
   module load nemsio/2.5.2
   module load sp/2.3.3
   module load w3emc/2.7.3
   module load w3nco/2.4.1
   module load sigio/2.3.2 sfcio/1.4.1 gfsio/1.4.1
   module load wrf_io/1.1.1

 elif [ $machine = hera.intel ] ; then

   module load cmake/3.20.1
   module use /scratch2/NCEPDEV/nwprod/hpc-stack/libs/hpc-stack/modulefiles/stack
   module load hpc/1.1.0
   module load hpc-intel/18.0.5.274
   module load hpc-impi/2018.0.4
   module load jasper/2.0.22
   module load zlib/1.2.11
   module load png/1.6.35
   module load hdf5/1.10.6
   module load netcdf/4.7.4
   module load bacio/2.4.1
   module load crtm/2.3.0
   module load g2/3.4.1
   module load g2tmpl/1.10.0
   module load ip/3.3.3
   module load nemsio/2.5.2
   module load sp/2.3.3
   module load w3emc/2.7.3
   module load w3nco/2.4.1
   module load sigio/2.3.2 sfcio/1.4.1 gfsio/1.4.1
   module load wrf_io/1.1.1
 else
   echo 'unsupported machine at this time'
   exit
 fi
   module list

  export CMAKE_C_COMPILER=mpiicc
  export CMAKE_CXX_COMPILER=mpiicpc
  export CMAKE_Fortran_COMPILER=mpiifort
  export CMAKE_Platform=$machine

# module show netcdf/4.7.4
# setenv CMAKE_C_COMPILER mpiicc
# setenv CMAKE_CXX_COMPILER mpiicpc
# setenv CMAKE_Fortran_COMPILER mpiifort
# setenv CMAKE_Platform $machine
  
#here is how to build post with cmake
  mkdir build
  cd build
  cmake .. -DBUILD_POSTEXEC=ON -DCMAKE_INSTALL_PREFIX=../install
  make -j8
  make install
  echo ' Installation of nceppost complete'
  exit

