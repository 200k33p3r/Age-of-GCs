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
feh=$8
deltaMag=$9
resample_id=${10}

local=/work2/08819/mying/$GC_name/simulateCMD
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

cp $local/define_errors.f90 $tmp
cp $local/random.f90 $tmp
cp $local/MakeTestCMD.f90 $tmp
#cp random.mod $tmp
#cp real_precision.mod $tmp
cp $local/real_precision.f $tmp
cp -r $local/inputfiles $tmp

cd $tmp
pwd

gfortran -c real_precision.f
gfortran -c define_errors.f90
gfortran -c random.f90
gfortran MakeTestCMD.f90 define_errors.o random.o real_precision.o -o MakeTestCMD

for ((r=frun; r <= lrun; r++ )); do
   echo "$inputiso/feh${feh}cmd.$r    $pdmf"
  time ./MakeTestCMD $inputiso/feh${feh}cmd.$r  $pdmf $binaryfrac  $numstars $cmdage $deltaMag
  mv $tmp/mc$r.a$cmdage $out/mc$r.a$cmdage$resample_id
done

cd $local
rm -r $tmp