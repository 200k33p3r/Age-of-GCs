#!/bin/bash

GC_name=$1
file_dir="${GC_name}_inputfiles"
file_temp="${GC_name}_inputfiles_temp"
mkdir ${file_dir}
mkdir ${file_temp}

cp define_mags.f90 ${file_temp}
cp PhotoPropI.f90 ${file_temp}
cp RadialDensity.f90 ${file_temp}
cp qsort.f90 ${file_temp}
cp sort.sh ${file_temp}
cd ${file_temp}

gfortran -c define_mags.f90
gfortran PhotoPropI.f90 define_mags.o -o PhotoPropI
gfortran -c qsort.f90
gfortran RadialDensity.f90 qsort.o -o RadialDensity

./PhotoPropI 99999 ../../M55_data/M55artstars.dat
./RadialDensity 13361 ../../M55_data/M55_fitstars.dat

./sort.sh

cp -r *s.dat ../${file_dir}
cp -r *Boundary.dat ../${file_dir}
cp fit_id.dat ../${file_dir}
cp Distance.dat ../${file_dir}
cd ..
rm -rf ${file_temp}



