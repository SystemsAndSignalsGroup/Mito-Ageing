#!/bin/sh
#PBS -N Alignment
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=1:mem=1gb


: '
UPDATES to 2.1

-Asks and checks for the name of the mitochondrial chromosome for the species you align to
	This is needed and checked for by mitogetter
-There was an error in the flow control when the pipeline checked if the STAR expression matrix existed
	This has been fixed
-The run_remap function erroneosly had the genome it used hardcoded. Now behaves as a more flexible functiont
-MAJOR CHANGE: sequencing_type=custom
	This should generalise the pipeline beyond the rigid framework of 10x
	Now makes use of umi_tools to Align with regular STAR, not STARsolo. Assuming 2 read files, one with cell barcode + umi and the other file with the genomic read, it works as follows:
	1) Takes the barcode fastq, strips the barcode and umi and places them, "_" delimited, into the corresponding genome read names. This is done using umi_tools extract
	2) Take the newly renomd genomic fastq and align it with STAR

	umi_tools allows the user a range of options for how to deduplicate the UMI reads. You need to think hard and carefully about the sequencing type you are using, look through the UMI tools manual, and select the appropriate option. You will need to go inside the runSTAR() function and enter these parameters yourself
	You will need to make these changes for both umi_tools count, and umi_tools dedup
	Steps 3) and 4) both make use of the barcode extraction in step 1

	3) Use featureCounts & umi_tools counts to obtain the expression matrix, using the deduplication procedure that you have just thought really carefully about
	4) Use umi_tools dedup to obtain the deduplicated BAM file also using the deduplication procedure you have just thought really carefully about

	In order to interact with the remapping procedure using minimap2, we need to make further modifications to the deduplicated bam we produced in 4)

	5) Use the readnames, which now contain the barcode and umi because of step 1), to put the cell barcode and umi as tagged fields in the bam file using custom_retag.py
	6) Run the remapper
	-Under the hood, the remapping.py will now use custom_renamer.py instead of renamer.py when doing this sequencing_type=custom


'



#To run the script use qsub -v chemistry="insert_10x_chemistry",remap="TRUE(TRUE is default/set FALSE if needed) 10x_full_pipeline.sh
#If tired of entering parameters at the command line, fix their values below and just qsub 10x_full_pipeline.sh


#Easily modifiable user parameters

#sequencing type. Choose either 10x or custom
seq_type="10x"
chemistry="2"

#if seq_type=custom custom, try some of the following defaults to modify
: '

seq_type="custom"
#In the barcode read file, we need to give barcode_umi structure. Use Letter C for every position with cell barcode and Letter N for every position with UMI" See the example below
barcode_structure="CCCCCCCCNNNN"
#example of whitelist you might use for a sequencing type"
whitelist=$PBS_O_WORKDIR/CEL_SEQ2_whitelist.txt

'


#remap="(*default TRUE)/FALSE"
#skip_freebayes="TRUE/(default FALSE)"
#Give the correct name of the mitochondrial chromosome in species you're aligning to
#Needed for mitogetter
MitoChromosome="chrM"



#set up the main variables before the main pipeline runs


#The output directory of the STAR aligned BAM files

HOME_DIR="$WORK/sc_mito_genet_pheno/Ma/wat"
LISTS_DIR="$HOME_DIR/file_finder"
#SAMPLE_GSM=$(head -$PBS_ARRAY_INDEX $LISTS_DIR/unique_gsm_list.txt | tail -1 )
SAMPLE_GSM="GSM4331851"
#Navigate to the GSM file containing all of the SRRs
DATA_DIR="$RDS_PROJECT/sc-ageing/ephemeral/Datasets/AidanMa/wat/${SAMPLE_GSM}"


: '
Ensure that you point to the correct refernce genome. It is easy to accidentally align to the wrong species. (Sigh)
'

# Align the reads in the Fastq files using STARsolo
GENOME_DIR="${RDS_PROJECT}/sc-ageing/live/refs/referenceGenomeRat"
GENOME_NAME="/rds/general/project/sc-ageing/live/refs/referenceGenomeRat/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa"
# .gtf file used by STARsolo
GTF="/rds/general/project/sc-ageing/live/refs/referenceGenomeRat/Rattus_norvegicus.Rnor_6.0.99.gtf"
OUT_DIR_ALIGN="${SAMPLE_GSM}/Align"
TEMP_DIR="${RDS_PROJECT}sc-ageing/ephemeral/STARalignerTemp/${SAMPLE_GSM}"
rm ${TEMP_DIR}* -rf

THEBAM="$DATA_DIR/${SAMPLE_GSM}"
#minimap2 bam file name
miniBAM="$DATA_DIR/$SAMPLE_GSM/minimap2/minimap_tagged_sorted.bam"
#minimap2 vcf files
miniVCF="$DATA_DIR/$SAMPLE_GSM/minimap2/merged_sorted_vcf.vcf.gz"


: '
STEP 0 - Preliminaries

The initial steps here check that the relevant paths have been set up and that the command line inputs make sense.
Namely, it checks:
1) Samtools check
2) The user has specified the sequencing type. If 10x chemistry needs to be specified aswell
3) If remapping checks the paths for the remapping done by minimap2
4) Check the relevant python modules are installed for the remap
5) Check the relevant paths are available for the MitoGetter.sh script to run
'

set -e

module load anaconda3/personal
source activate pipeline2.0

#Define the main functions that will be used to run the pipeline

preliminary_check () {
	#samtools check
	# Only the latest version of Samtools worked and conda update wouldn't work either. Hence the export
	echo "Prelim check 1: Samtools"
	export PATH=/rds/general/user/asm119/home/sc_mito_genet_pheno/samtools_install/bin:$PATH
	samcheck=$(samtools --version | grep "samtools")
	echo "	Using $samcheck"
	echo "	If samtools older than 1.10, please update, then proceed"



	echo "Prelim check 2: seq check"

	if [ $seq_type == 10x ]; then
		echo "	10x Sequencing protocol chosen"
		echo " 	Checking chemistry settings"
		#STAR bam file name
		STARBAM="$THEBAM/AlignAligned.sortedByCoord.out.bam"
		Expression_matrix="$THEBAM/AlignSolo.out/Gene/filtered/matrix.mtx"
		if [ $chemistry == 2 ]; then
			whitelist="$PBS_O_WORKDIR/10xv2_whitelist.txt"
			UMIlen="10"
			echo "	Pipeline assuming 10x chemistry v2"
			echo "	Using Whitelist: $whitelist"
			echo " 	prelim check 2 complete. Proceeding"
		elif [ $chemistry == 3 ]; then
			whitelist="$PBS_O_WORKDIR/10xv3_whitelist.txt"
			UMIlen="12"
			echo "	Pipeline assuming 10x chemistry v3"
			echo "	Using Whitelist: $whitelist"
			echo "	prelim check2 complete. Proceeding"
		else
			echo "	No 10x chemistry specified. ABORTING!"
			echo "	Specify a chemistry"
			exit 111
		fi
	elif [ $seq_type == custom ]; then
		echo "	Custom sequencing protocol chosen"
		echo "	Checking minimal settings selected"
		#STAR bam file name
		STARBAM="$THEBAM/sorted_retagged_deduped.bam"
		Expression_matrix="$THEBAM/counts.tsv.gz"
		if [ -n $whitelist ]; then
			echo "	Using whitelist: $whitelist"
		else
			echo "	No Whitelist Provided. Assuming error. ABORTING!"
			echo "	Please provide a whitelist or modify pipeline to ignore"
			exit 111
		fi
		if [ -n $barcode_structure ]; then
			echo "	Barcode Structure: $barcode_structure"
			echo "	prelim check2 complete. Proceeding"
		else
			echo "	No barcode structure provided. ABORTING!"
			echo "	Please give barcode+umi structure"
			exit 111
		fi
	else
		echo "	Sequencing type not specified. ABORTING!"
		echo "	Specify a sequencing type. Either 10x or custom"
		exit 111
	fi




	echo "Prelim check 3 & 4: Remapping"
	#remapping path environment checking
	DO_REMAP=${remap:-TRUE}
	if [ $DO_REMAP == TRUE ]; then
		#export the path to the remapping tools
		export PATH=$PBS_O_WORKDIR/remap_package:$PATH
		#Check the export worked
		hash remapping.py 2>/dev/null || { echo >&2 "remapping.py not found in the path. ABORTING!"; exit 1; }
		hash renamer.py 2>/dev/null || { echo >&2 "renamer.py not found in the path. ABORTING!"; exit 1; }
		hash retag.py 2>/dev/null || { echo >&2 "retag.py not found in the path. ABORTING!"; exit 1; }
		echo "	Remap sequence: Activated"
		echo "	prelim check 3 complte. proceeding"
		module_check.py
		echo "	all python modules installed"
		#The default should lead to FREEBAYES running
		SKIP_FREEBAYES=${skip_freebayes:-FALSE}
			echo "	skip_freebayes = $SKIP_FREEBAYES"
		echo "	prelim check 4 complete. proceeding"
	elif [ $DO_REMAP == FALSE ]; then
		echo "	Remap sequence: Deativated"
		echo "	prelim steps 3 & 4 complete. proceeding"
	else
		echo "	Invalid remap parameter value. ABORTING!"
		echo "	Remap parameter must either be TRUE or FALSE in all caps"
	fi


	echo "Prelim check 5: MitoGetter.sh"
	export PATH=$PBS_O_WORKDIR:$PATH
	hash split_script.py 2>/dev/null || { echo >&2 "split_script.py not found in the path. Aborting!"; exit 1; }
	hash MitoGetter.sh 2>/dev/null || { echo >&2 "MitoGetter.sh not found in the path. Aborting!"; exit 1; }
	if [ -n "$MitoChromosome" ]; then
		echo "	Mitochromosome name: $MitoChromosome"
	else
		echo "	Mitochromosome not specified: ABORTING!"
		echo "	Specify Mitchondrial Chromosome name!"
		exit 111
	fi
	echo "	prelim check 5 complete. proceeding"
}







: '
STEP 1 - ALIGNMENT

Now perform an alignment to the Genome using STARsolo.

One of the read files contains reads Cell Barcode and UMIs. Another contains the Genome reads.

For 10x Genomics, the fastq files should be in the order --readFilesIn genome.fastq.gz Barcode+UMI.fastq.gz. See the STARmanual, section 13 STARsolo

The whitelist of Barcodes that you provide to STARsolo depends on the "chemistry" used in sequencing. Google around until you get the right whitelist for the sequencing protocol that was used to produce the FASTQ files. We should have a copy of this whitelist available locally.

STARsolo should output a BAM file of the genome reads, each read tagged with Barcode and UMI separately

The parameter --soloCBstart and --soloCBlen require the start barcode start position and length respectively. This also depends on the sequencing technique used. Google around till you find the appropriate values. Dito fot the UMI equivalents
'


runSTAR () {
	echo "STEP 1:STARsolo Alignment - commencing"

	# I'm time stamping the different steps to see where the bottlenecks in this process are
	now="$(date +'%m/%d/%Y')"
	now2="$(date +'%r')"
	now="$(date +'%m/%d/%Y')"
	now2="$(date +'%r')"
	echo "$now"
	echo "$now2"


	: '
	In this version of the pipeline STARsolo aligns all of the SRRs within a single GSM file simultaneously. These are assumed to be from multiple runs on the same librar. Therefore sequence them all together and make sure you deduplicate all the UMIs (For the love of god make sure you deduplicate)
	The long, rather forbidding looking pipeing of commands below gets the addresses of all the SRRs within the same GSM and puts it in a format that STARsolo can work with
	'



	#Get the addresses of all the read1 files in the GSM. This will be fed into STARsolo
	read1s=$(ls ${DATA_DIR}/SRR*/*_1.fastq.gz | sed '$!s/$/,/' | tr -d '\n')
	#ditto for the read2 files
	read2s=$(ls ${DATA_DIR}/SRR*/*_2.fastq.gz | sed '$!s/$/,/' | tr -d '\n')
	echo "	Read file destinations:"
	echo "$read1s"
	echo "$read2s"

	: '
	In version 2.1 of the pipeline STAR alignment forks depending on the whether 10x or custom sequencing option is chosen. If 10x is chosen, then the conventional STARsolo we know and love is used. If custom is used, then a new process is used, which makes use of some other tools, mainly umi_tools
	'
	if [ $seq_type == 10x ]; then
		#Run STARsolo
		mkdir -p $OUT_DIR_ALIGN
		STAR --runThreadN 16 \
		--genomeDir $GENOME_DIR \
		--sjdbGTFfile $GTF \
		--readFilesCommand gunzip -c \
		--soloBarcodeReadLength 0 \
		--readFilesIn $read2s $read1s \
		--outTmpDir $TEMP_DIR \
		--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
		--outSAMtype BAM SortedByCoordinate \
		--soloUMIdedup 1MM_All \
		--outFileNamePrefix $OUT_DIR_ALIGN \
		--soloType Droplet \
		--soloCBwhitelist $whitelist \
		--soloCBstart 1 \
		--soloCBlen 16 \
		--soloUMIstart 17 \
		--soloUMIlen $UMIlen
		echo "Done."

		echo "	STAR Alignment - complete"


		: '
		STARsolo automatically outputs 3 files required to make an expression matrix, barcodes.tsv  features.tsv  matrix.mtx.
		In order for scanpy to work with these files, features.tsv needs to be renamed genes.tsv (I have no idea why, it just does)
		I therefore rename the file here
		'


		FILTERED_MATRIX="$SAMPLE_GSM/AlignSolo.out/Gene/filtered/"

		mv ${FILTERED_MATRIX}features.tsv ${FILTERED_MATRIX}genes.tsv
		echo "	Copying STAR outputs to $DATA_DIR"
		#rm -r $DATA_DIR/$SAMPLE_SRR
		#Copy the outputs of STAR back to the DATA_DIR
		cp -r ${SAMPLE_GSM} $DATA_DIR/
		echo "	Copy to STAR outputs to ephemeral space: success"
		echo "STEP 1 - complete"
	fi

	if [ $seq_type == custom ]; then
		#step 1, get umi tools to extract the umis from the read to the read name
		#cel-seq2 as used by Muraro
		#make the output into a new folder
		mkdir -p $DATA_DIR/${SAMPLE_SRR}/umitools
		#extract the barcode and UMI into the read names into the genomic read file
		umi_tools extract -I $DATA_DIR/${SAMPLE_SRR}/${SAMPLE_SRR}_1.fastq.gz \
		--bc-pattern=CCCCCCCCNNNN --log=processed.log \
		-S $DATA_DIR/${SAMPLE_SRR}/umitools/readname_${SAMPLE_SRR}_1.fastq.gz \
		--read2-in=$DATA_DIR/${SAMPLE_SRR}/${SAMPLE_SRR}_2.fastq.gz \
		--read2-out=$DATA_DIR/${SAMPLE_SRR}/umitools/readname_${SAMPLE_SRR}_2.fastq.gz \
		--filter-cell-barcode \
		--whitelist=$PBS_O_WORKDIR/Muraro_whitelist.txt

		#Get the addresses of the read files. These are the files with the extracted names we just got from umi_tools
		read1s=$(ls $DATA_DIR/${SAMPLE_SRR}/umitools/readname_${SAMPLE_SRR}_1.fastq.gz | sed '$!s/$/,/' | tr -d '\n')
		#ditto for the read2 files
		read2s=$(ls $DATA_DIR/${SAMPLE_SRR}/umitools/readname_${SAMPLE_SRR}_2.fastq.gz | sed '$!s/$/,/' | tr -d '\n')
		echo "	Read file destinations:"
		echo "$read1s"
		echo "$read2s"

		#Run STARsolo
		mkdir -p $OUT_DIR_ALIGN
		STAR --runThreadN 16 \
		--genomeDir $GENOME_DIR \
		--sjdbGTFfile $GTF \
		--readFilesCommand gunzip -c \
		--readFilesIn $read2s \
		--outTmpDir $TEMP_DIR \
		--outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix $OUT_DIR_ALIGN
		echo "Done."

		echo "	STAR Alignment - complete"

		#We've made it this far. If all is well, should have a BAM file with where read names contain the UMI
		#Can now use a couple of steps to get the deduped bam and expression matrix
		#First use feature counts to output a bam with a tag of where a gene maps to
		featureCounts -a $GTF -o gene_assigned -R BAM $THEBAM/AlignAligned.sortedByCoord.out.bam -T 16
		#sort and index again
		samtools sort AlignAligned.sortedByCoord.out.bam.featureCounts.bam -o assigned_sorted.bam
		samtools index assigned_sorted.bam
		#obtain the Expression matrix using count, to get counts per gene per cell for CEL-SEQ2 protocol
		umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I assigned_sorted.bam -S $THEBAM/counts.tsv.gz

		#We also need a deduplicated BAM file for variant calling and heteroplasmy estimation
		#We also place tags into this so that the file can be used for minimap2 remapping

		#Deduplicate the bam for variant calling using dedup
		umi_tools dedup --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I assigned_sorted.bam --output-stats=deduplicated -S $THEBAM/deduplicated.bam
		#Now there's a multistep sequence to put cell barcode and UMI into tags in the bam file
		#1) convert the bam to sam
		samtools view -h -o out.sam assigned_sorted.bam
		#2) Use the retag script to put the tags into a new bam
		celseq2_retag.py -s out.sam -o retagged_deduped.bam
		#3) The usual sort and index again
		samtools sort retagged_deduped.bam -o sorted_retagged_deduped.bam
		samtools index sorted_retagged_deduped.bam
		#4) Save this to the output space
		cp sorted_retagged_deduped.bam $THEBAM/
	fi


}




run_remap () {

	: '
	STEP 2 - Remapping with minimap2 and variant calling
	'

	echo "STEP 2: Minimap2 Remapping & Variant calling - commencing"
	# second time stamp
	now="$(date +'%m/%d/%Y')"
	now2="$(date +'%r')"
	echo "$now"
	echo "$now2"

	: '
	Making use of some scripts inspired by Souporcell pipeline, we realign the BAM files output from STARsolo with minimap2. Minimap2 is argued to be produce fewer false positive variants than STAR for variant calling, which will be done later in the pipeline

	This will only run if the DO_REMAP, paramter is TRUE, as determined by the user.
	'

	#Give the arguments for the remapper
	#Get the STAR aligned .bam to be remapped with minimap
	remap_BAMFILE=$1
	#Get the Barcodes to be remapped
	remap_barcodes=$2
	#Get the location of the .fa file which is used by minimap2
	remap_GENOMEfile=$3

	if [ $DO_REMAP == TRUE ]; then
		#run the remapping script. This should do a bunch of stuff to get the minimap2 realigned bam files
		#make make minimap2 output and a temporary file for bcftools to interact with
		mkdir -p $DATA_DIR/$SAMPLE_GSM/minimap2/temporary
		echo "	temporary file created"
		echo "	--skip_freebayes = $SKIP_FREEBAYES"
		echo "	commencing script"
		remapping.py -i $remap_BAMFILE \
		-b $remap_barcodes \
		-f $remap_GENOMEfile \
		-t 16 \
		-o $DATA_DIR/$SAMPLE_GSM/minimap2 \
		--skip_freebayes $SKIP_FREEBAYES \
		-s $seq_type
		echo "	script complete, and files should be in remap folder"
		echo "STEP 2 - complete"
	elif [ $DO_REMAP == FALSE ]; then
		echo "	Remapping deactivated. Skipping to STEP 3"
		echo "STEP 2 - complete"
	else
		echo "Invalid remap parameter value. ABORTING!"
		echo "If the pipeline got this far, you shouldn't be able to read this!"
		echo "Something terrible must have happened. Good Luck!"
		exit 111
	fi

}




run_MitoGetter () {
	# third time stamp
	now="$(date +'%m/%d/%Y')"
	now2="$(date +'%r')"
	echo "$now"
	echo "$now2"

	echo "STEP 3 - commencing"

	: ' STEP 3 - Getting only the mitochondrial reads
	'

	: '
	Here we call the MitoGetter.sh script.
	Should go through the full bam files and create a new bamfile for each cell containing only the mitochondrial reads

	This will be run for the STAR aligned and Minimap2 aligned BAM files. Modify this according to your neads

	Ensure the mitochondrial chromosome is named correctly for your species

	'

	#Get the .bam file to be worked on
	BAMFILE=$1
	#Get the name of the mitochondrial chromosome
	ChrMt=$2
	#Output Directory
	Output=$3
	#barcode whitelist used
	whitelist=$4

	echo "	Getting STAR aligned mitochondrial reads per cell"
	#split the star BAM, and put it into the star folder
	MitoGetter.sh $BAMFILE $ChrMt $Output $whitelist
	echo "	STAR mitoReads done"




	echo "run_MitoGetter - Complete"
}




#MAIN PIPELINE RUN
#0 Run the preliminary sequence
preliminary_check

#set the counter
c="0"
#1 Run STARsolo
if [ ! -d $THEBAM ]; then
	runSTAR
else
	echo "Pipeline restarting"
	echo "	BAM folder detected"
	c="1"
	echo "	Checking for STAR .bam file"
	if [ -f $STARBAM ]; then
		echo "	STAR Aligned .bam file detected"
	else
		echo "	STAR Aligned .bam file not detected. Aborting!"
		exit 111
	fi

	echo "	Checking for STAR expression matrix"

	if [ -f $Expression_matrix ]; then
		echo "	STAR expression matrix found"
	else
		echo " 	STAR expression matrix not detected. Aborting!"
		exit 111
	fi

	echo "STEP 1 already complete. Moving to STEP 2"
fi

#2 Run the remap
#Set the Freebayes parameter
if [ $DO_REMAP == "TRUE" ]; then
	if [ ! -f $miniBAM ]; then
	STARbarcodes="$DATA_DIR/$SAMPLE_GSM/AlignSolo.out/Gene/filtered/barcodes.tsv"
		run_remap
	else
		echo "	minimap2 Aligned .bam file detected"
		if [ $SKIP_FREEBAYES == "TRUE" ]; then
			echo "	skip_freebayes selected"
			echo "	Moving to STEP 3"
		else
			if [ -f $miniVCF ]; then
				echo "	minimap2 vcfs detected"
				echo " 	Moving to STEP 3"
			else
				echo " 	No vcf detected"
				echo "	Running Freebayes"
				run_remap $STARBAM $STARbarcodes $GENOME_NAME
			fi
		fi
	fi
fi

#3 Run MitoGetter for the STAR aligned bam
if [ ! -d $THEBAM/MitoReads ]; then
	echo "	Getting STAR aligned mitochondrial reads per cell"
	run_MitoGetter $STARBAM $MitoChromosome $THEBAM $whitelist
	echo "	STAR mitoReads done"
else
	echo "	STAR MitoReads directory detected"
	count=$(ls $THEBAM/MitoReads | wc -l)
	if [ $counts == "0" ]; then
		echo "	MitoReads file found, but empty"
		echo "	Assuming timeout error occurred in previous run"
		echo "	Re-running MitoGetter"
		run_MitoGetter $STARBAM $MitoChromosome $THEBAM $whitelist
	else
		echo "	STAR MitoReads already split"

	fi

fi

#4 Run MitoGetter for the minimap2 Aligned bam file

if [ $DO_REMAP == "TRUE" ]; then
	if [ ! -d $DATA_DIR/$SAMPLE_GSM/minimap2/MitoReads ];then
		echo "	Getting minimap aligned mitochondrial reads per cell"
		run_MitoGetter $miniBAM $MitoChromosome $THEBAM/minimap2 $whitelist
		echo "	minimap2 mitoReads done"
	else
		echo "	minimap2 MitoReads directory detected"
		count=$(ls $DATA_DIR/$SAMPLE_GSM/minimap2/MitoReads | wc -l)
		if [ $counts == "0" ]; then
			echo "	MitoReads file found, but empty"
			echo "	Assuming timeout error occurred in previous run"
			echo "	Re-running MitoGetter"
			run_MitoGetter $miniBAM $MitoChromosome $THEBAM/minimap2 $whitelist
		else
			echo "	minimap2 MitoReads already split"
		fi
	fi
elif [ $DO_REMAP == "FALSE" ]; then
	echo "	Remap sequence: Deativated"
	echo "	Ending"
else
	echo "	Invalid remap parameter value. ABORTING!"
	echo "	If the pipeline got this far, you shouldn't be able to read this!"
	echo "	Something terrible must have happened. Good Luck!"
fi

#all complete
echo "PIPELINE COMPLETE!"
# successful exit code
exit 0

