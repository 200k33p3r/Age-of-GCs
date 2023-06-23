#!/bin/bash

GC_name=$1
feh=$2

progdir=/home/mying/Desktop/Age-of-GCs/CalibrationStarTest
data_dir=/data/
isodir=${data_dir}${GC_name}_data/${GC_name}_iso


gfortran $progdir/test_iso.f90 -o $progdir/test_iso

echo "#chi2 FeH Yprim alphafe CMIXLA FGRY FGRZ KTTAU ALPHAE ALPHAC PP He3+He3 He3+He4 P+C12 P+C13 P+N14 P+O16 alexcoef opalcoef2 talphacoef plascoef cocoef MCnumber" >chi2.dat
for filename in $isodir/*l.* ; do
    $mc_num=${filename#*.}
    $progdir/test_iso  $GC_name $feh $mc_num $data_dir >>chi2.dat
done
