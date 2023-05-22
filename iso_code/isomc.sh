#! /bin/csh
#calculate isochrones using *.iso output files 
#from a red giant branch monte carlo simulation

 setenv local  /home/mying/Desktop/iso_code
 setenv track  /data/M55_tracks
 setenv var   /data/M55_var
 setenv run   $local
 setenv out   $local/outiso

#define other variables
 @ feh = 190

#read from the command line which runs to do
 if ($#argv != 2) then
    echo "usage: isomc2.com 'first run' 'last run'"
    exit 1
 endif

 @ frun = $1
 @ lrun = $2

 rm fort.*
 @ r = $frun

 while ($r <= $lrun)
 @ good_track = 0

#input namelist files
    ln -s  $run/isomc.nml                     fort.95

    ln -s $var/varfeh-${feh}.$r                      fort.75

    tar -xf $track/iso_mcfeh${feh}.$r.tar 

#input evolutionary track files which are evolved up the RGB
    @ fn = 10
    @ tn = 1
    foreach m (0200 0240 0280 0320 0360 0400 0440 0480 0520 0560 0600 0640 0680 0700 0750 0800 0850 0900 0950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1500 1600 1700 1800 1900 2000 2200 2400 2600 2800 3000)
      @ fn = 10 + $tn
      gzip -d iso_m$m.feh${feh}.$r.gz
      @ size=`ls -l iso_m$m.feh${feh}.$r | awk '{print $5}'`
      #exclude failed models
      if ($size > 1000) then
         ln -s  iso_m$m.feh${feh}.$r               fort.$fn
         @ tn = $tn + 1
         @ good_track = $good_track + 1
      endif
    end
    #update the number of track in input.nml
    sed -i "8s/.*/  NTRK = ${good_track}/" isomc.nml

#output file
    ln -s  $out/feh${feh}.$r                     fort.92

    time $run/isomc2

    rm fort.*

# now, replace the low mass points with the FreeEOS runs
    ln -s $var/varfeh-${feh}.$r                      fort.75
    ln -s  $out/feh${feh}.$r                     fort.11

    @ tn = 1
    @ good_track = 0
    foreach m (0200 0240 0280 0320 0360 0400 0440 0480 0520 0560 0600 0640 0680)
      @ fn = 20 + $tn
      @ size=`ls -l iso_m$m.feh${feh}.$r | awk '{print $5}'`
      #exclude failed models
      if ($size > 1000) then
         ln -s  iso_m$m.feh${feh}.$r               fort.$fn
         @ tn = $tn + 1
         @ good_track = $good_track + 1
      endif
    end
#update the number of track in input.nml
    sed -i "9s/.*/  NTRKLOW = ${good_track}/" isomc.nml
#output file
    ln -s  $out/feh${feh}l.$r       fort.12

    set ntrklow=`grep NTRKLOW isomc.nml  | awk '{print $3}' `
    set numage=`grep Numage isomc.nml   | awk '{print $3}' `


    time $run/lowmass2 $ntrklow $numage 0.691

    rm fort.*
    rm iso_m*.feh${feh}.$r


# finally, create a file suitable for fake CMD generation
    ln -s $out/feh${feh}l.$r              fort.11
# output file
     ln -s  $out/feh${feh}cmd.$r       fort.12

     time $run/prepiso $numage 0.691

     rm fort.*

  

 @ r++
end


