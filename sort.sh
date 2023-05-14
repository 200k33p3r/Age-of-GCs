#!/bin/bash
for filename in V*.dat; do
 outfile=${filename%.dat}
 outfile="${outfile}s.dat"
 echo "File: $filename $outfile"
 (head -n 3 $filename && tail -n +4 $filename | sort -g ) > $outfile
done

for filename in I*.dat; do
 outfile=${filename%.dat}
 outfile="${outfile}s.dat"
 echo "File: $filename $outfile"
 (head -n 3 $filename && tail -n +4 $filename | sort -g ) > $outfile
done