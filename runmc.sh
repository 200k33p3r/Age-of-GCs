#! /bin/bash
#
# script which automatically runs the Dartmouth stellar evolution code  
# for a Monte Carlo simulation
# requires program newpoly 
# requires subsrunmc.sh 
# requires gs98getz to change [Fe/H] & [alpha/Fe]
#

LD_LIBRARY_PATH=/dartfs-hpc/rc/home/f/f004qnf/dsep3/libs/free_eos-2.2.1/build/:$LD_LIBRARY_PATH

local=/dartfs-hpc/rc/lab/C/ChaboyerB/Catherine/dsep3

poly=$local/poly
zams=$local/zams
prems=$local/prems
opac=$local/opac
atm=$local/surfBC
opal=$opac/OPAL
phx=$opac/phoenix
bcphx=$atm/phoenix/GS98
bckur=$atm/atmk

#mc directory
mcdir=/dartfs-hpc/rc/home/f/f004qnf/mcdir
nml=$mcdir/nml
prog=/dartfs-hpc/rc/lab/C/ChaboyerB/Catherine/dsep3/dsepX
#use variable files Rowan generated
vardir=/dartfs-hpc/rc/home/v/f0056bv/MC_var
rundir=/dartfs-hpc/rc/lab/C/ChaboyerB/Catherine/run

#read from command line which run to do
if [ $# -ne 2 ]; then
    echo "You need to input a var file number"
    exit 1
fi

# create a  random working directory

number=$RANDOM
num2=$(date +%s)
num2=$((num2 - 1265900000))

out=/scratch/$number$num2
echo "Scratch dir: $out"
mkdir $out

# directory which output gets copied to from the temp directory
storeoutput=/dartfs-hpc/rc/lab/C/ChaboyerB/Catherine/M55
echo "Output dir: $storeoutput"



ID=GS98
#generate mass from 0.2 to 3 solar mass
Lmass=($(seq 200 40 690))
MLmass=($(seq 700 50 1400))
MHmass=($(seq 1500 100 2000))
Hmass=($(seq 2000 200 3000))
masses=( "${Lmass[@]}" "${MLmass[@]}" "${MHmass[@]}" "${Hmass[@]}" )
ovshoot=0.0
nmod=9999
endage=1.8D10

runum=$1
tmpfeh=${2#-}
absfeh=`echo "$tmpfeh * 100" | bc`
intfeh=${absfeh%.*}
if (($tmpfeh > 0)); then
    varfile=$vardir/varfeh${intfeh}.$runum
else
    varfile=$vardir/varfeh-${intfeh}.$runum
fi
echo "$varfile"

feh=$(cat $varfile | awk 'NR==1 {print $1}')
echo "[Fe/H] = $feh"
afe=$(cat $varfile | awk 'NR==3 {print $1}')
echo "[a/Fe] = $afe"
yy=$(cat $varfile | awk 'NR==2 {print $1}')
echo "Yp = $yy"
mixlen=$(cat $varfile | awk 'NR==4 {print $1}')
echo "mixL = $mixlen"

zz=$($rundir/gs98getz $feh $afe)
echo "Z= $zz"

xx=$(awk "BEGIN {print 1.0-$zz-$yy}")
echo "X= $xx"

#input various subroutines

. $rundir/subsrunmc.sh


cp $rundir/newpoly $out/.
cd $out

#setup a bunch of directories, filenames and parameters for newpoly
setups

#setup mixture specific variables
setupmix

#calculate adjustment to starting Teff & Lum based upon Z of the model
# note that Z is calculated in the subroutine setupmix
#bash can't do real number arithmetic, so have to use awk

# These are the originals
#Ladjust=$(awk "BEGIN {print -0.115*log($zz/$zsolar)/log(10.0)}")
#Tadjust=$(awk "BEGIN {print -0.022*log($zz/$zsolar)/log(10.0)}")

# These are Karabo's and they don't work for high metallicity
#Ladjust=$(awk "BEGIN {print -0.015*log($zz/$zsolar)/log(10.0)}")
#Tadjust=$(awk "BEGIN {print -0.010*log($zz/$zsolar)/log(10.0)}")
#Tadjust=$(awk "BEGIN {print -0.010*log($zz/$zsolar)/log(10.0) + 0.1*($mixlen-1.9)}")

# These are mine and I'm not sure what they are for...
Ladjust=$(awk "BEGIN {print 0.13141*$mixlen+1.5197}")
Tadjust=$(awk "BEGIN {print 0.04724*$mixlen+3.5936}")

ln -s $varfile    fort.85


for mass in "${masses[@]}"; do

echo "mass = $mass"

#calculate Teff, and luminosity  of the starting pre-MS model
    sumass=$(awk "BEGIN {print $mass/1000.0}")
    get_teff_lum
    echo "Mass= $sumass  lumonisty=$lum  Teff=$teff"

#write out the namelist and  run the program
    wrtpoly
    echo "Finished wrtpoly"
    ./newpoly
    rm poly.nml

#run stellar evolution code

#assign physics namelist based upon mass
    getphysnml

    if (($mass < 1000 )); then
	mprefix=0$mass
    else
	mprefix=$mass
    fi


    fname=m$mprefix.feh$intfeh.$runum
    ln -s $out/$fname fort.35

    echo "about to start DSEP"
    time $prog/dsepX
    gzip $out/$fname
#    mv $out/$fname.gz $storeoutput/.
 
    rm fort.{11,12,13,35,19,20,22}
 
done
 
cd $out
 
tar -cf mcfeh$intfeh.$runum.tar m*.*.$runum.gz
 
mv mcfeh$intfeh.$runum.tar  $storeoutput/.
 
rm -rf $out
 
exit 0
