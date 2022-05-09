# Script to run pyGenomeTracks to show coverage log2ratio signal with Del/Dup regions and add additional gene track to show ABC transporters 

# install pyGenomeTracks via conda 
conda create -n pygenometracks -c bioconda -c conda-forge pygenometracks python=3.7


# Working directory is Google\ Drive/My\ Drive/HiC_files

# activate environment 
source activate pygenometracks

# set data directory for bigwig files 
DIRIN=results/WGS/03_dedup_sorted_bw
# set output directory, also directory holding bed files 
DIROUT=results/WGS/results

# resolution parameter 
RES=10000
# data type parameter
DATA=wgs

# create configuration file describing the tracks using make_tracks_file for chromosome-specific plots 
# add the coverage tracks (PR then CR), then the log2ratio BW, then bed file 
### We will use the bed file with the prefix new_color_ because we changed the colors to red and blue 
make_tracks_file --trackFiles ${DIRIN}/UCD52PR_wgs_10000.bw ${DIRIN}/UCD52CR_wgs_10000.bw ${DIRIN}/UCD52CR_vs_PR_${DATA}_${RES}_log2ratio.bw ${DIROUT}/new_color_UCD52CR_vs_PR_${DATA}_${RES}_log2ratio.bed -o ${DIROUT}/wgs_tracks.ini
# comment out the line that allows for the bed file for the duplications and deletions to be colored based on the color column 
# comment out min_value = 0 for bigwig 
# 



## Deletions Visualizations 
### CHR3 0-198,295,559

pyGenomeTracks --tracks ${DIROUT}/wgs_tracks.ini --region chr3:0-198,295,559  --height 7 --width 20 --outFileName ${DIROUT}/wgs_chr3_DEL_DUP_visualization.svg


### CHR5 0-181,000,000

pyGenomeTracks --tracks ${DIROUT}/wgs_tracks.ini --region chr5:0-181,000,000 --height 7 --width 20 --outFileName ${DIROUT}/wgs_chr5_DEL_visualization.svg


### CHR20 chr20	0-64,000,000
pyGenomeTracks --tracks ${DIROUT}/wgs_tracks.ini --region chr20:0-64,000,000 --height 7 --width 20 --outFileName ${DIROUT}/wgs_chr20_DEL_visualization.svg



### CHR 22 chr22	0-50,810,000
pyGenomeTracks --tracks ${DIROUT}/wgs_tracks.ini --region chr22:0-50,810,000 --height 7 --width 20 --outFileName ${DIROUT}/wgs_chr22_DEL_visualization.svg


## Duplications Visualizations 
### CHR4 0-190,000,000
pyGenomeTracks --tracks ${DIROUT}/wgs_tracks.ini --region chr4:0-190,000,000 --height 7 --width 20 --outFileName ${DIROUT}/wgs_chr4_DUP_DEL_visualization.svg



### CHR6 0-170,000,000
pyGenomeTracks --tracks ${DIROUT}/wgs_tracks.ini --region chr6:0-170,000,000 --height 7 --width 20 --outFileName ${DIROUT}/wgs_chr6_DUP_visualization.svg


### CHR17 0-83,000,000
# chr17	64,930,001	77,320,001 --> dups
# chr17	56,520,001	64,920,001 --> dups

pyGenomeTracks --tracks ${DIROUT}/wgs_tracks.ini --region chr17:0-83,000,000 --height 7 --width 20 --outFileName ${DIROUT}/wgs_chr17_DEL_DUP_visualization.svg



##### *** Zoom up chr17 duplication, add genes track and highlight ABC transporters 
# Create new tracks configuration file adding in the gene tracks 
make_tracks_file --trackFiles ${DIRIN}/UCD52PR_wgs_10000.bw ${DIRIN}/UCD52CR_wgs_10000.bw ${DIRIN}/UCD52CR_vs_PR_${DATA}_${RES}_log2ratio.bw ${DIROUT}/new_color_UCD52CR_vs_PR_${DATA}_${RES}_log2ratio.bed ${DIROUT}/hg38_genes.bed -o ${DIROUT}/wgs_genes_tracks.ini

# comment out the line that allows for the bed file for the duplications and deletions to be colored based on the color column 
# comment out in_value = 0 for bigwig 

# test one really zoomed in to look at genes
pyGenomeTracks --tracks ${DIROUT}/wgs_genes_tracks.ini --region chr17:68,000,000-70,000,000 --outFileName ${DIROUT}/wgs_chr17_DUP_Genes_zoomed_visualization.svg
