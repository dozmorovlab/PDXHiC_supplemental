#!/bin/bash
#PBS -S /bin/bash
#PBS -l nodes=1
#PBS -M mdozmorov@vcu.edu
#PBS -m abe
#PBS -N HiCrep_A
#PBS -j oe
#PBS -q workq
# PBS -V

cd $PBS_O_WORKDIR

source hicrepfix/bin/activate

export HDF5_USE_FILE_LOCKING='FALSE'

DIRIN=/home/sequencing/juicer/Mikhail/PDXproject_analysis/01.Arima_replicates_hic2cool
DIROUT=/home/sequencing/juicer/Mikhail/PDXproject_analysis/01.Arima_replicates_juicer

array=( 2437_01 2437_02 2437_03 2437_04 )

total=$((${#array[@]}-1)) # Minus 1 to match array indices
# total=17
for i in $(seq 0 $total); do
    for j in $(seq 0 $total); do
        if [ $j -gt $i ] 
        then
            FILEIN1=$DIRIN"/"${array[$i]}".mcool"
            FILEIN2=$DIRIN"/"${array[$j]}".mcool"
            FILEOUT=${array[$i]}"_"${array[$j]}".txt"
            hicrep  --binSize 10000 --h 10 --dBPMax 5000000 $FILEIN1 $FILEIN2 $FILEOUT
        fi
    done
done

