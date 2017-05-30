#!/bin/bash

###############################################################################################################################
### Preliminary analysis for the Data-Integration (ENCODE/Roadmap Epigenomics) project: #######################################
###############################################################################################################################
### statistical survey on the overlap between different ChIP-seq data and genes from the GenCode database   @Francesco Gandolfi
###############################################################################################################################

#NB (Task); Convert processed BED files into Normalized Signal tracks.


# 0.1 define root directorties
PROJECT_DIR=/home/fgandolfi/projects/ENCODE_project
DATA_DIR=${PROJECT_DIR}/data;
SCRIPT_DIR=/home/fgandolfi/scripts/ENCODE_project

# 0.2 Analysis directory
ANALYSIS_DIR=${PROJECT_DIR}/genomic_survey;
mkdir ${ANALYSIS_DIR};

# 0.3 Sample-datasheet
DATASHEET_SAMPLES_FILE=${DATA_DIR}/ENCODE_project_IMR90.tsv;
#DATASHEET_SAMPLES_FILE=${DATA_DIR}/ENCODE_project_datasheet_input_9-2015.tsv


# 0.4 Uniqueness mappability track
MAPPABILITY_TRACK=/home/fgandolfi/templates/wgEncodeDukeMapabilityUniqueness35bp.uniqueMapRegions.bedGraph

# 0.5 Chrominfo file
CHROMINFO_FILEPATH=${DATA_DIR}/chromInfo.txt

# 0.5 Be sure that you already have the k-binned genome file (k=bin size)
#/usr/bin/bedtools makewindows -g ${CHROMINFO_FILEPATH} -w 200 > ${ANALYSIS_DIR}/hg19binned.200bp.bed

# 0.6 Genomic bin size
GENOME_BINSIZE=200;

# 0.7 Define Binned hg19 genome path
BINNED_GENOME=${ANALYSIS_DIR}/hg19binned.200bp.bed









##############################
####### LOOP-1 ###############
##############################

# 1. Calculate the total number of informative reads mapped across all the datasets

# Take note of the number of informative reads in that sample to estimate total number of tags over all experiments in a cell line
# Define cell types
CELL_LINES=IMR90;
# Define signal-types
SIGNAL_TRACKS=(H3K27ac H3K9me3 H3K79me1 Pol2);

TAG_TOTALS=();  # should be the vector of totals across all ChIPseq experiments for each cell line.
#sep=' ';

for CELL_LINE in ${CELL_LINES[@]};do
	

mkdir ${ANALYSIS_DIR}/${CELL_LINE}
	
	for SIGNAL_TRACK in ${SIGNAL_TRACKS[@]};do

# 2.1 Define the output directory
OUT_DIR=${ANALYSIS_DIR}/${CELL_LINE}/${SIGNAL_TRACK};
mkdir ${OUT_DIR};
	
	awk -F "\t" '$1=="'${CELL_LINE}'" && $3=="'${SIGNAL_TRACK}'" {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' ${DATASHEET_SAMPLES_FILE} > ${ANALYSIS_DIR}/ENCODE_project_datasheet_sub.tsv;
	
	TAG_TOTALS=();  # should be the total of informative reads across all ChIP-seq samples in a given sample group (cell line/signal-track)
	sep=' ';
		
	
			while read LINE;do
				
SAMPLE_ID=$(echo ${LINE} | awk '{split($0,a," ");print a[4]}');
# 1.3.3 Identify and retrieve processed.bed file and its path
PROCESSED_FILEDIR=${DATA_DIR}/${CELL_LINE}/${SIGNAL_TRACK}/${SAMPLE_ID};	
cd ${PROCESSED_FILEDIR};
PROCESSED_FILE=*.processed.bed
filelines=$(wc -l ${PROCESSED_FILE})
N_INFORMATIVE_TAGS=$(echo ${filelines} | awk '{split($0,a," ");print a[1]}');
echo "Increase the total count of informative reads across all ChIP-seq assay samples in ${CELL_LINE}: $SIGNAL_TRACK";
TAG_TOTALS+=$N_INFORMATIVE_TAGS;
TAG_TOTALS+=$sep;	

			done < ${ANALYSIS_DIR}/ENCODE_project_datasheet_sub.tsv;
	
	echo ${TAG_TOTALS} > ${OUT_DIR}/${CELL_LINE}_${SIGNAL_TRACK}_allExp_totals.txt;  # txt file contains the total read counts for all samples Bi in the group P.
	done;


done;









##############################
####### LOOP-2 ###############
##############################


# 2. Once the cell-type and the signal-type are defined, extract all samples from the cell-type X and signal-type Y from the ENCODE_project_datasheet_samples.tsv
#     then loop row by row on the subset 

for CELL_LINE in ${CELL_LINES[@]};do


	for SIGNAL_TRACK in ${SIGNAL_TRACKS[@]};do

# 2.1 Regenerate OUTPUT directory:
OUT_DIR=${ANALYSIS_DIR}/${CELL_LINE}/${SIGNAL_TRACK};
	
	
# 2.2 Create an empty the list of files
> ${ANALYSIS_DIR}/input_data_files.txt
	
# 2.3 Select from the full-datasheet file only raw corresponding to CELL_LINE and SIGNAL_TRACK
awk -F "\t" '$1=="'${CELL_LINE}'" && $3=="'${SIGNAL_TRACK}'" {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9}' ${DATASHEET_SAMPLES_FILE} > ${ANALYSIS_DIR}/ENCODE_project_datasheet_sub.tsv	
	
	
	
# Read through the datasheet subset (foreach sample Bi belonging to the group X...)
while read LINE;do

# 2.3.1 get values in single fields;
SAMPLE_ID=$(echo ${LINE} | awk '{split($0,a," ");print a[4]}')
ALIGN_FILENAME=$(echo ${LINE} | awk '{split($0,a," ");print a[5]}')
	
# 2.3.2 check if the alignment file is available or not
if [ -z "$ALIGN_FILENAME" ] 
then
continue
fi 	
	
# 2.3.3 Identify and retrieve processed.bed file and its path
PROCESSED_FILEDIR=${DATA_DIR}/${CELL_LINE}/${SIGNAL_TRACK}/${SAMPLE_ID};	
cd ${PROCESSED_FILEDIR};
PROCESSED_FILE=*processed.bed
PROCESSED_FILEPATH=${PROCESSED_FILEDIR}/${PROCESSED_FILE};

# 2.3.4 Estimate raw counts in each genomic bin;
echo "intersect genomic bins with processed bed ${PROCESSED_FILEPATH} to estimate rawCounts"
/usr/bin/bedtools intersect -a ${ANALYSIS_DIR}/hg19binned.200bp.bed -b ${PROCESSED_FILEPATH} -c > ${OUT_DIR}/${SAMPLE_ID}.rawCounts.bed
# 2.3.5 Extract only informative bins in 'temp.bed' , and replace 'rawCounts' by this one.
awk -F "\t" '$4 > 0 {print $1 "\t" $2 "\t" $3 "\t" $4}' ${OUT_DIR}/${SAMPLE_ID}.rawCounts.bed > ${OUT_DIR}/temp.bed
mv ${OUT_DIR}/temp.bed ${OUT_DIR}/${SAMPLE_ID}.rawCounts.bed
# 2.3.6 Put the rawCount file name in a list for normalized signal estimation step
echo "${OUT_DIR}/${SAMPLE_ID}.rawCounts.bed" >> ${ANALYSIS_DIR}/input_data_files.txt

done < ${ANALYSIS_DIR}/ENCODE_project_datasheet_sub.tsv;
	
	
	
	
# 2.4 Process rawCounts; the R script takes the list of rawCount files just generated.
#INPUT: rawCounts files for a given group (cell-lineX/signaltrackY) /// OUTPUT: normalizedSignal files for the group (cell-lineX/signalrackY)

# arguments are: 1) datafile with raw counts 
				# 2) datasheet_samples_file.tsv 
				# 3) Mappability track file 
				# 4) genomic bin sizes 
				# 5) total read counts per cell type.txt 
/usr/bin/Rscript ${SCRIPT_DIR}/binCount2normalizedSignal.R ${ANALYSIS_DIR}/input_data_files.txt ${DATASHEET_SAMPLES_FILE} ${MAPPABILITY_TRACK} ${BINNED_GENOME} ${GENOME_BINSIZE} ${OUT_DIR}/${CELL_LINE}_${SIGNAL_TRACK}_allExp_totals.txt;
# 2.5. Estimate normalizedSignal values from all samples (replicates/experiments from different labs) pooled together
cd ${OUT_DIR}
BEDG_EXPECTED_FILES=*expectedScore.bg
/usr/bin/bedtools unionbedg -i ${BEDG_EXPECTED_FILES[@]} > ${OUT_DIR}/expectedScoreSum.bg;
BEDG_OBSERVED_FILES=*observedScore.bg
/usr/bin/bedtools unionbedg -i ${BEDG_OBSERVED_FILES[@]} > ${OUT_DIR}/observedScoreSum.bg;

/usr/bin/Rscript ${SCRIPT_DIR}/getPooledNormalizedSignal.R ${OUT_DIR}/expectedScoreSum.bg ${OUT_DIR}/observedScoreSum.bg;
# Remove intermediate files:

done
done
	


### End script ####
