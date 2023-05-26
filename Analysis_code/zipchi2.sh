#! /bin/bash
#this code is used to create zipped files from the chi2 files
#which can make files to be moved easier

GC_name='M55'

chi2_dir=/work2/08819/mying/$GC_name/outchi2
zipped_dir=/work2/08819/mying/$GC_name/chi2_zips

ages=($(seq 8000 200 16001))
prefix='chi2_a'
mc_prefix='_mc'

frun=$1
lrun=$2

if [ $# -ne 2 ]; then
    echo "usage: 'first mc_num' 'last mc_num'"
    exit 1
fi 
cd ${chi2_dir}
for ((r=frun; r <= lrun; r++ )); do
    echo "Working on mc_num=$r"
    for age in "${ages[@]}"; do
        if (($age < 10000)); then
        age_str=0$age
        else
        age_str=$age
        fi
        fname=${prefix}${age_str}${mc_prefix}${r}
        gzip $chi2_dir/$fname
    done
    
    tar -cf chi2_mc${r}.tar *${mc_prefix}${r}.gz
    mv chi2_mc${r}.tar  $zipped_dir/.
done
