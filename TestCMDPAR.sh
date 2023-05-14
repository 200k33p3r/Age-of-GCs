#! /bin/bash 
# generate simulated CMDs from MC isochornes
# assumes input isochrones are in subdirectory ./inputiso
# assumes input photmetric error files are in subdirectory ./inputfiles
# MY 01/22 Code is modified to be able to run in parallel on Discovery

frun=$1
lrun=$2
pdmf=$3
binaryfrac=$4
numstars=$5
cmdage=$6
GC_name=$7

local=/work2/08819/mying/$GC_name/simulateCMD
#inputiso=$local/inputiso
inputiso=/work2/08819/mying/$GC_name/outiso
out=$local/outcmd

#create a folder in scratch
number=$RANDOM
num3=$RANDOM
num2=$(date +%s)
num2=$((num2 - 1265900000))

tmp=/scratch/08819/mying/$number$num3$num2
#check if folder exists
while [ -f "$tmp" ]
do
    $num2=$((num2 + 1))
    $tmp=/scratch/08819/mying/$number$num2
done

mkdir $tmp

cp define_errors.f90 $tmp
cp random.f90 $tmp
cp MakeTestCMD.f90 $tmp
#cp random.mod $tmp
#cp real_precision.mod $tmp
cp real_precision.f $tmp
cp -r inputfiles $tmp

cd $tmp
pwd

gfortran -c real_precision.f
gfortran -c define_errors.f90
gfortran -c random.f90
gfortran MakeTestCMD.f90 define_errors.o random.o real_precision.o -o MakeTestCMD



#read from the command line which runs to do and simulation parameters
if [ $# -ne 6 ]; then
    echo "usage: SimulateCMD.sh 'first run' 'last run' 'PDMF slope' 'binary fraction' 'number simulated stars' 'cmdage'"
    exit 1
fi 

#if [[ -e fort.* ]]; then
#    rm fort.*
#fi


for ((r=frun; r <= lrun; r++ )); do
   echo "$inputiso/feh230cmd.$r    $pdmf"
  time ./MakeTestCMD $inputiso/feh230cmd.$r  $pdmf $binaryfrac  $numstars $cmdage
  mv $tmp/mc$r.* $out/.
done

cd $local
rm -r $tmp