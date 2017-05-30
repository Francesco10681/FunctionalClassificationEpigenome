#!/bin/bash

###############################################################################################################################
### Preliminary analysis for the Data-Integration (ENCODE/Roadmap Epigenomics) project: #######################################
###############################################################################################################################
### statistical survey on the overlap between different ChIP-seq data and genes from the GenCode database   @Francesco Gandolfi
###############################################################################################################################


#NB; (Task); Process raw BAM files into filtered BED files (duplicates removal, quality alignment filtering, 200bp-tag extension).




#################################################
## phase 1: prepare bed files from .bam files####
#################################################

# keep only uniquely mapped reads in .bam file;;
# remove duplicate reads from .bam file;;
# sort and index .bam file;;
# convert each .bam into .bed file;;


# 0. define root directorties
PROJECT_DIR=/home/fgandolfi/projects/ENCODE_project
DATA_DIR=${PROJECT_DIR}/data;
DATASHEET_SAMPLES_FILE=${DATA_DIR}/ENCODE_project_IMR90.tsv;
cd ${DATA_DIR};


# 1.1 read the datasheet-tab-delimited file by line
while read LINE;do

# get values in single fields;
	CELL_LINE=$(echo ${LINE} | awk '{split($0,a," ");print a[1]}');
	ASSAY=$(echo ${LINE} | awk '{split($0,a," ");print a[3]}');
	SAMPLE_ID=$(echo ${LINE} | awk '{split($0,a," ");print a[4]}');
	ALIGN_FILENAME=$(echo ${LINE} | awk '{split($0,a," ");print a[5]}');
	ALIGN_FILETYPE=$(echo ${LINE} | awk '{split($0,a," ");print a[6]}');

# 1.2 check if the alignment file is available or not
if [ -z "$ALIGN_FILENAME" ] 
then
continue
fi 


#progress
echo "working on sample ${SAMPLE_ID} assay ${ASSAY} cell-line ${CELL_LINE}";



# 1.3 check if processed.bed file is already provided and in this case re-format this.
if [ $ALIGN_FILETYPE = "processed" ]
then
	echo "processed bed file $ALIGN_FILENAME";

# 1.4 Replace empty-space SF by tabs
	echo "replace empty spaces with tabs..";
	sed 's/ /\t/g' ${DATA_DIR}/${CELL_LINE}/${ASSAY}/${SAMPLE_ID}/${ALIGN_FILENAME} > ${DATA_DIR}/${CELL_LINE}/${ASSAY}/${SAMPLE_ID}/temp0.bed

	echo "verify columns in the bed format..";
	NCOL=$(awk -F '\t' '{print NF;exit}' ${DATA_DIR}/${CELL_LINE}/${ASSAY}/${SAMPLE_ID}/temp0.bed)
	if [ "$NCOL" -eq 5 ]
	then
	echo "adjust number of columns";
		sed -i "s/$/\t255/" ${DATA_DIR}/${CELL_LINE}/${ASSAY}/${SAMPLE_ID}/temp0.bed;
		awk -F "\t" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $6 "\t" $5}' ${DATA_DIR}/${CELL_LINE}/${ASSAY}/${SAMPLE_ID}/temp0.bed > ${DATA_DIR}/${CELL_LINE}/${ASSAY}/${SAMPLE_ID}/temp1.bed;
	else
	echo "processed bed file $ALIGN_FILENAME is correctly formatted";
	cp ${DATA_DIR}/${CELL_LINE}/${ASSAY}/${SAMPLE_ID}/temp0.bed ${DATA_DIR}/${CELL_LINE}/${ASSAY}/${SAMPLE_ID}/temp1.bed
	fi



# 1.5 Add 'processed' suffix to the processed .bed file
PROCESSED_FILENAME=${ALIGN_FILENAME/".bed"/".processed.bed"};
# Remove anomalous position
echo "remove anomalous position along the chromosomes..";
awk -F "\t" '$2 >= 0 {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' ${DATA_DIR}/${CELL_LINE}/${ASSAY}/${SAMPLE_ID}/temp1.bed > ${DATA_DIR}/${CELL_LINE}/${ASSAY}/${SAMPLE_ID}/${PROCESSED_FILENAME};
# Cancel temporary files
rm ${DATA_DIR}/${CELL_LINE}/${ASSAY}/${SAMPLE_ID}/temp0.bed;rm ${DATA_DIR}/${CELL_LINE}/${ASSAY}/${SAMPLE_ID}/temp1.bed;
echo "alignment file is ${ALIGN_FILETYPE} and skip";
continue
fi  # close block for 'processed bed files'



# @bam files or .bed not-processed@:
echo "raw bed OR bam file $ALIGN_FILENAME";
CELL_LINE_DIR=${DATA_DIR}/${CELL_LINE};
SAMPLE_DIR=${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID};
# get the path to the alignment file
ALIGN_FILEPATH=${SAMPLE_DIR}/${ALIGN_FILENAME};


			# take note of every bam file contained;
			cd ${SAMPLE_DIR};
			#echo "DIRECTORY is ${DATA_DIR}/${CELL_LINE}/${ASSAY}/${SAMPLE}";
			
			# If there is a .bed file (this is in raw format), convert it into .bam
			bedcount=`ls -1 *.bed 2>/dev/null | wc -l`
			if [ "$bedcount" -gt 0 ] 
			then
			echo "convert bed to bam..";
            # define the output file name just generated
			BAM_FILE=${ALIGN_FILENAME/".bed"/".bam"};
			/usr/bin/bedtools bedtobam -i ${ALIGN_FILEPATH} -g ${DATA_DIR}/chromInfo.txt > ${BAM_FILE};
			else
			BAM_FILE=${ALIGN_FILENAME};
			fi
			
			
#NB: Starting from here we always have a .bam file!
			
			# Retrieve only the sample name:
			BAM_NAME=${BAM_FILE/".bam"/""};

			# sort and index bam file;
			echo "sort and index ${BAM_NAME}";
		    /usr/bin/samtools sort ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_FILE} ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.sorted  
		    /usr/bin/samtools index ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.sorted.bam 

			# remove duplicate reads;
			echo "mark and remove duplicates in ${BAM_NAME}";
		    /usr/bin/picard-tools MarkDuplicates INPUT=${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.sorted.bam OUTPUT=${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.nodup.bam ASSUME_SORTED=True REMOVE_DUPLICATES=True METRICS_FILE=${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.metrics.log VALIDATION_STRINGENCY=SILENT
		
		    # take only reliable alignments ( only matches with 99% probability of being real matches are kept )
			echo "make a new bam-index and keep only high-quality alignments in ${BAM_NAME}";
			/usr/bin/samtools index ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.nodup.bam
			/usr/bin/samtools view -q 20 -b ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.nodup.bam > ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.filtered.bam
			# Alternatively, if Bowtie was used, you select uniquematch reads:
			# /usr/bin/samtools view -f 0x2 -q 255 -b tumor1.bam -o unique.bam
			
			# convert .bam to .bed (bedtools bamtobed);
			echo "convert ${BAM_NAME} into .bed format";
			/usr/bin/bedtools bamtobed -i ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.filtered.bam  > ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.processed.bed
		
			# For TFs and HistMod signals apply tag extension of 200bp 
			## split by strand +/- ; 
			### apply tag extension of 200bp; 
			#### merge extended bed files
			
			cd ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID};
			
			if [ "$ASSAY" != "DNase-seq" ]
			then
			echo "apply tag extension on ${BAM_NAME}";
		    awk -F "\t" '$6=="+" {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.processed.bed > ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.pos.bed
			awk -F "\t" '$6=="-" {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.processed.bed > ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.neg.bed	
			awk -F "\t" '{print $1,$2,$2+199,$4,$5,$6}' ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.pos.bed > ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.extpos.bed
			awk -F "\t" '{print $1,$3-199,$3,$4,$5,$6}' ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.neg.bed > ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.extneg.bed
			cat ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.extneg.bed ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.extpos.bed > ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.processed.bed		
			
			# remove intermediate files
			rm ${BAM_NAME}.ext*; rm ${BAM_NAME}.pos.bed;rm ${BAM_NAME}.neg.bed;
			# otherwise sort the .bed file only
			else
			echo "sort by strand,chromosome,pos $BAM_NAME";
			sort -k6,6 -k1,1 -k2,2n ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.processed.bed -o ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.processed.bed
			
			fi
		
			# introduce tab-delimited separator
			awk '$1=$1' FS=" " OFS="\t" ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.processed.bed > ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/out.bed
			
			# Remove anomalous position
			awk -F "\t" '$2 >= 0 {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/out.bed > ${CELL_LINE_DIR}/${ASSAY}/${SAMPLE_ID}/${BAM_NAME}.processed.bed
			
			
			# remove intermediate files
			rm ${BAM_NAME}.nodup*; rm ${BAM_NAME}.sorted*; rm ${BAM_NAME}.filtered.bam; rm out.bed;
					
					
done<${DATASHEET_SAMPLES_FILE};


###################
### End script ####
##################



