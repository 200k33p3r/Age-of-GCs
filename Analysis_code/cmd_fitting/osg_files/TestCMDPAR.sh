#! /bin/bash 
# generate simulated CMDs from MC isochornes
# assumes input isochrones are in subdirectory ./inputiso
# assumes input photmetric error files are in subdirectory ./inputfiles
# MY 01/22 Code is modified to be able to run in parallel on Discovery
#Modify the code to run on OSG

frun=$1
lrun=$2
pdmf=$3
binaryfrac=$4
numstars=$5
cmdage=$6
GC_name=$7
feh=$8

local=.
inputiso=.
out=.

#read from the command line which runs to do and simulation parameters
if [ $# -ne 8 ]; then
    echo "usage: SimulateCMD.sh 'first run' 'last run' 'PDMF slope' 'binary fraction' 'number simulated stars' 'cmdage' 'GC_names' 'feh'"
    exit 1
fi 

for ((r=frun; r <= lrun; r++ )); do
   echo "$inputiso/feh${feh}cmd.$r    $pdmf"
  time ./MakeTestCMD $inputiso/feh${feh}cmd.$r  $pdmf $binaryfrac  $numstars $cmdage
done