#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=1
#PBS -M mdozmorov@vcu.edu
#PBS -m ae
#PBS -N htseq
#PBS -q workq
#PBS -j oe
# PBS -o /home/glasser/Dissertation/Dozmorov/RNA-seq/cuffdiffout3

cd $PBS_O_WORKDIR

source activate htseq

# FILEGTF=/home/sequencing/data/ExtData/UCSC/rn6/Rattus_norvegicus.Rnor_6.0.85_chr.gtf
FILEGTF=/home/sequencing/data/ExtData/UCSC/hg38gdc/gencode.v39.annotation.gtf

# SAMPLE=UCD52PR
# # SAMPLE=UCD52CR

# DIRIN=/home/sequencing/juicer/Mikhail/WGS/03_dedup_sorted
# FILEIN=${DIRIN}/${SAMPLE}/${SAMPLE}.dedup.sorted.bam

# DIROUT=/home/sequencing/juicer/Mikhail/WGS/05_htseq
# mkdir -p ${DIROUT}
# FILEOUT=${DIROUT}/${SAMPLE}.txt

# htseq-count -f bam -r pos -s no ${FILEIN} ${FILEGTF} > ${FILEOUT}

# SAMPLE=UCD52PR
SAMPLE=UCD52CR

DIRIN=/home/sequencing/juicer/Mikhail/WGS/03_dedup_sorted
FILEIN=${DIRIN}/${SAMPLE}/${SAMPLE}.dedup.sorted.bam

DIROUT=/home/sequencing/juicer/Mikhail/WGS/05_htseq_CR
mkdir -p ${DIROUT}
FILEOUT=${DIROUT}/${SAMPLE}.txt

htseq-count -f bam -r pos -s no ${FILEIN} ${FILEGTF} > ${FILEOUT}
