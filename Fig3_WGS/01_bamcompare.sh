#!/bin/bash
#PBS -S /bin/bash
# PBS -V
#PBS -l nodes=1:ppn=24
#PBS -M mdozmorov@vcu.edu
#PBS -m ae
#PBS -N bamcompare
#PBS -j oe

cd $PBS_O_WORKDIR

source activate deeptools
# Deduplicated and sorted bam files were converted to bigwig files for visualization in IGV

DIRIN=/home/sequencing/juicer/Mikhail/WGS/03_dedup_sorted
FILEIN1=${DIRIN}/UCD52PR/UCD52PR.dedup.sorted.bam
FILEIN2=${DIRIN}/UCD52CR/UCD52CR.dedup.sorted.bam
DIROUT=/home/sequencing/juicer/Mikhail/WGS/03_dedup_sorted_bw
mkdir -p $DIROUT

EXCLUDE=/home/sequencing/juicer/Mikhail/WGS/exclude.cnvnator_100bp.GRCh38.20170403.bed
THREADS=24

# Multiple resolutions, in bases
RESOLUTION=( 10000 100 1000 ) # 100000

for RES in ${RESOLUTION[@]}; do
  # log2ratio calculation
	FILEOUT=${DIROUT}/UCD52CR_vs_PR_wgs_${RES}_log2ratio.bw
	bamCompare -b1 ${FILEIN2} -b2 ${FILEIN1} --binSize ${RES} --blackListFileName ${EXCLUDE} --numberOfProcessors ${THREADS} --effectiveGenomeSize 2913022398 -o ${FILEOUT}
  # Coverage for individual files
	bamCoverage --bam $FILEIN1 -o ${DIROUT}/$(basename $FILEIN1 .dedup.sorted.bam)_wgs_${RES}.bw --binSize ${RES};
	bamCoverage --bam $FILEIN2 -o ${DIROUT}/$(basename $FILEIN2 .dedup.sorted.bam)_wgs_${RES}.bw --binSize ${RES};
done

conda deactivate

# Convert all files to bedGraph format
source activate samtools

for RES in ${RESOLUTION[@]}; do
	FILEOUT=${DIROUT}/UCD52CR_vs_PR_wgs_${RES}_log2ratio.bw
	FILEOUT2=${DIROUT}/`basename ${FILEOUT} .bw`.bedGraph
	bigWigToBedGraph ${FILEOUT} ${FILEOUT2}
	echo ${FILEOUT} ${FILEOUT2}
	gzip ${FILEOUT2}

	FILEOUT=${DIROUT}/$(basename $FILEIN1 .dedup.sorted.bam)_wgs_${RES}.bw
	FILEOUT2=${DIROUT}/`basename ${FILEOUT} .bw`.bedGraph
	bigWigToBedGraph ${FILEOUT} ${FILEOUT2}
	echo ${FILEOUT} ${FILEOUT2}
	gzip ${FILEOUT2}

	FILEOUT=${DIROUT}/$(basename $FILEIN2 .dedup.sorted.bam)_wgs_${RES}.bw
	FILEOUT2=${DIROUT}/`basename ${FILEOUT} .bw`.bedGraph
	bigWigToBedGraph ${FILEOUT} ${FILEOUT2}
	echo ${FILEOUT} ${FILEOUT2}
	gzip ${FILEOUT2}
done

