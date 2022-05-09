# Script to use LoopWidth https://github.com/magmarshh/LoopWidth to compute loop widths for 
# Mustache loops from preprocessing_any folder for each condition and resolution

# Working directory is /Users/MarshallMa4

# Conditions, correspond to file names 
FILE=( PR_unique CR_unique PR_common CR_common )

# Mustache Loops Data Directory 
LOOPDIR=Google\ Drive/My\ Drive/HiC_files/results/mustache_results/preprocessing_any
# Results directory
DIROUT=Google\ Drive/My\ Drive/HiC_files/results/Maggie/Mustache_results/preprocessing_any/Loops/Sizes_and_CTCF_Orientation
# Loop resolution
LOOPRES=( 10kb 25kb 50kb 100kb )
LOOPNUM=( 10000 25000 50000 100000 )

# Loop through every condition and resolution and output loop piecharts
for res in {1..4}; do
	for f in ${FILE[@]}; do
		python Desktop/LoopWidth/LoopWidth_piechart.py --loop ${LOOPDIR}/loops_${f}_${LOOPRES[$res]}.bedpe --res ${LOOPNUM[$res]} \
		--output ${DIROUT}/${LOOPRES[$res]}/Mustache_${f}_${LOOPRES[$res]}_loop_sizes_piechart.svg;
	done;
done




