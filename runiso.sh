#! /bin/bash
#
# scipt to run Isochrone code (Dotter 2016) using stellar evolution 
# models generated by Dartmouth Stellar Evolution Program
# Used particularly for HSTACS data and require the corresponding 
# color tables
# makecmd is bugged currently because the bc_tables only have [alpha/Fe] = 0
#

#define working directory
local='/home/mying/Desktop/iso_D'
DSEP_output='/data/M55_tracks'
tracks="$local/data/tracks"
eeps="$local/data/eeps"
isochrones="$local/data/isochrones"
bc_tables="$local/HST_ACSWF"

#read mc_num
start=$1
end=$2

#define variables
num_tracks=40
feh=190
Lmass=($(seq 200 40 690))
MLmass=($(seq 700 50 1400))
MHmass=($(seq 1500 100 2000))
Hmass=($(seq 2000 200 3000))
masses=( "${Lmass[@]}" "${MLmass[@]}" "${MHmass[@]}" "${Hmass[@]}" )

if [ $# -ne 2 ]; then
    echo "You need to input the start id and end id"
    exit 1
fi

#write bc_table.list
if [ ! -f bc_table.list ]; then
    echo "18" >bc_table.list
    echo "${bc_tables}/fehm400" >>bc_table.list
    echo "${bc_tables}/fehm350" >>bc_table.list
    echo "${bc_tables}/fehm300" >>bc_table.list
    echo "${bc_tables}/fehm275" >>bc_table.list
    echo "${bc_tables}/fehm250" >>bc_table.list
    echo "${bc_tables}/fehm225" >>bc_table.list
    echo "${bc_tables}/fehm200" >>bc_table.list
    echo "${bc_tables}/fehm175" >>bc_table.list
    echo "${bc_tables}/fehm150" >>bc_table.list
    echo "${bc_tables}/fehm125" >>bc_table.list
    echo "${bc_tables}/fehm100" >>bc_table.list
    echo "${bc_tables}/fehm075" >>bc_table.list
    echo "${bc_tables}/fehm050" >>bc_table.list
    echo "${bc_tables}/fehm025" >>bc_table.list
    echo "${bc_tables}/fehp000" >>bc_table.list
    echo "${bc_tables}/fehp025" >>bc_table.list
    echo "${bc_tables}/fehp050" >>bc_table.list
    echo "${bc_tables}/fehp075" >>bc_table.list
fi

ids=($(seq $start 1 $end))
for id in "${ids[@]}"; do
    echo "MC_num = $id"
    #extract files to 
    tar -xf ${DSEP_output}/mcfeh${feh}.${id}.tar
    gzip -d *.feh${feh}.${id}.gz
    mv *.feh${feh}.${id} $tracks

    #write input files
    if [[ -e input.isoD ]]; then
	    rm input.isoD
    fi
    echo "#version string, max 8 characters" >input.isoD
    echo "isoD" >>input.isoD
    echo "#initial Y, initial Z, [Fe/H], [alpha/Fe], v/vcrit (space separated)" >>input.isoD
    echo "        0           0   -1.90        0.80       0" >>input.isoD
    echo "#data directories: 1) history files, 2) eeps, 3) isochrones" >>input.isoD
    echo "$tracks" >>input.isoD
    echo "$eeps" >>input.isoD
    echo "$isochrones" >>input.isoD
    echo "#read history_columns" >>input.isoD
    echo "my_history_columns_basic.list" >>input.isoD
    echo "#specify tracks" >>input.isoD
    echo "$num_tracks"  >>input.isoD
    for mass in "${masses[@]}"; do
        if (($mass < 1000 )); then
	        mprefix=0$mass
        else
	        mprefix=$mass
        fi
        echo "m${mprefix}.feh${feh}.${id}"  >>input.isoD
    done
    echo "#specify isochrones"  >>input.isoD
    echo "isochrones_${id}.txt"  >>input.isoD
    echo "min_max"  >>input.isoD
    echo "linear"  >>input.isoD
    echo "41"  >>input.isoD
    echo "8e9"  >>input.isoD
    echo "16e9"  >>input.isoD
    #run isochrone code
    ./make_eep input.isoD
    ./make_iso input.isoD
    ./make_cmd HST_ACSWF "${isochrones}/isochrones_${id}.txt"
done

#remove tracks
rm -rf $tracks/*.feh${feh}.${id}
