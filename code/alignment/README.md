# Alignment steps

The scripts in this directory take the user from the raw fastq files and outputs a directory in the containing a bam file for each and every valid barcode that contains that cell's mitochondrial reads. The main skeleton of the pipeline is implemented within 10x_PlusCustomBarcode_full_pipeline_perGSM.sh. This long file name reflects that this script is capable of taking 10x chromium data, as well as other library preparation approaches --- such as cellseq2 --- with distinct barcode structures, hence *PlusCustomBarcode*. The file is called per GSM to note that for many datasets downloaded from GEO, which are grouped with a unique GSM number, typically need to have all fastq files from that GSM aligned and demultiplexed by STARSolo as they are often different runs from the same experiment. 

The relevant conda environment should work with the packages versions as specified in the pipeline2.yml file. Additionally, a version of samtools later than 1.10. If not, then the functionality for sorting cells by their barcode will not work. If a user runs 10x_PlusCustomBarcode_full_pipeline_perGSM.sh, the pipeline should actually check if the environment has the necessary scripts to run. If not, an error will be thrown saying what is missing. Users will also need to get the barcode whitelists for the relevant protocols they work on. Eg, for 10x genomics, details about whitelists can be found here https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist-.

Once the relevant packages have been installed a user should be able to modify the script to suit their needs. At the beginning of the 10x_PlusCustomBarcode_full_pipeline_perGSM.sh script there are numerous parameters. Each of these should be commented sufficiently for the user to identify their meaning. The most important starting point is the DATA_DIR variable which should point at a directory containing all of the fastq files from within a GEO GSM. Another things to be wary of is the variable **MitoChromosome**. This pesky variable takes the name of the chromosome for your organism of interest and has different names for different species.

Once these preliminaries have begun the pipeline operates on a 5-step procedure. 

0. Preliminaries - Checks all the relevant installs and helpfully points out errors if there are any and offers solutions
1. Alignment - with STARSolo and modifying using the parameters from the beginning of the pipeline as inputs.
2. Remapping -There is a second optional remapping step to realign with bam file with FREEBAYES. This is not used for the purposes of the manuscript. This will make use of the other scripts found within ./remap_package directory
3. MitoGetter - Executes the MitoGetter script which takes the aligned bam file output from STARSolo, strips out the mitochondrial reads of this bam file, sorts the mitochondrial bam by the cell barcode,  and uses split_script.py to split the mitochondrial bam file into separate bam files for each cell. These are the files used for variant calling and heteroplasmy estimates in downstream analysis for the rest of the manuscript.
4. Remapped Mitogetter - Exectues Mitogetter on the realigned bam file from FREEBAYES if there users wanted this option.


These steps are packaged into functions within the script and are then executed according to the flow control sequence at the end of the script. This flow control factors in parameter choices at the beginning of the script (eg. whether to remap or not). Additionally, if certain output files are found up to certain steps, this flow control will take skip to the latest part of the sequence without the relevant scripts.



We ran the scripts on a PBS cluster scheduler on the Imperial College London hpc. Operability may differ on other cluster types.

