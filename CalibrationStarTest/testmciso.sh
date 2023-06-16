#!/bin/bash

progdir=/Users/chaboyer/gaia/dr3/hst
isodir=$progdir/outiso


gfortran $progdir/test_iso.f90 -o $progdir/test_iso

echo "#chi2 FeH Yprim alphafe CMIXLA FGRY FGRZ KTTAU ALPHAE PP He3+He3 He3+He4 P+C12 P+C13 P+N14 P+O16 alexcoef opalcoef2 talphacoef plascoef cocoef Colour MCnumber" >chi2.dat
for filename in $isodir/feh230isol.* ; do
    $progdir/test_iso  $filename >>chi2.dat
done
