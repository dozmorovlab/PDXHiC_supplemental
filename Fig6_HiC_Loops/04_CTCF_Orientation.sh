# Script to compute CTCF orientation using Maggie and Pratiik's CTCF orientation tool at https://github.com/magmarshh/CTCF_orientation 
# for Mustache loops from preprocessing_any folder 

# hg38 CTCF fimo.bed file will be used 

### Working directory is ~/Users/MarshallMa4/

# Files, same for each data type 
FILE=( PR_unique CR_unique PR_common CR_common )
# hg38 fimo file 
FIMODIR=Google\ Drive/My\ Drive/HiC_files/results/Maggie/CTCF_orientation_misc/fimo.bed

# Mustache Loops Data Directory where the loop files are stored 
LOOPDIR=Google\ Drive/My\ Drive/HiC_files/results/mustache_results/preprocessing_any
# Mustache results directory 
DIROUT=Google\ Drive/My\ Drive/HiC_files/results/Maggie/Mustache_results/preprocessing_any/Loops/Sizes_and_CTCF_Orientation
# Loop resolution
LOOPRES=( 10kb 25kb 50kb 100kb )

# Loop through every condition and resolution and produce the CTCF piechart
for res in {1..4}; do
	for f in ${FILE[@]}; do
		python Google\ Drive/My\ Drive//CTCF_orientation/ctcf_orientation.py -l ${LOOPDIR}/loops_${f}_${LOOPRES[$res]}.bedpe \
		-m ${FIMODIR} -o ${DIROUT}/${LOOPRES[$res]}/Mustache_${f}_${LOOPRES[$res]}_CTCF_orientation.svg;
	done;
done
