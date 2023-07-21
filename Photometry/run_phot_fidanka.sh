#!/bin/bash
#this code use fiducial isochrones instead of AS test to determine the photometric error

GC_name=$1
file_dir="${GC_name}_fidanka_inputfiles"
file_temp="${GC_name}_fidanka_inputfiles_temp"
mkdir ${file_dir}
mkdir ${file_temp}

#write the result of fiducial isochrone fitting into AS test file format
python phot_fidanka.py $GC_name

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

AS_path=../../${GC_name}_data/${GC_name}_fiducial_artstars.dat
obs_path=../../${GC_name}_data/${GC_name}_fitstars.dat

AS_stars=$((`wc --lines < $AS_path` - 1))
obs_stars=$((`wc --lines < $obs_path` - 1))

./PhotoPropI $AS_stars $AS_path False
./RadialDensity $obs_stars $obs_path

./sort.sh

cp -r *s.dat ../${file_dir}
cp -r *Boundary.dat ../${file_dir}
cp fit_id.dat ../${file_dir}
cp Distance.dat ../${file_dir}
cd ..
rm -rf ${file_temp}
