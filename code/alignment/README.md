# Alignment steps

The scripts in this directory take the user from the raw fastq files and outputs a directory in the containing a bam file for each and every valid barcode that contains that cell's mitochondrial reads. The main skeleton of the pipeline is implemented within 10x_PlusCustomBarcode_full_pipeline_perGSM.sh. This long file name reflects that this script is capable of taking 10x chromium data, as well as other library preparation approaches --- such as cellseq2 --- with distinct barcode structures, hence *PlusCustomBarcode*. The file is called per GSM to note that for many datasets downloaded from GEO, which are grouped with a unique GSM number, typically need to have all fastq files from that GSM aligned and demultiplexed by STARSolo as they are often different runs from the same experiment. 

The relevant conda environment should work with the packages versions as specified in the pipeline2.yml file. Additionally, a version of samtools later than 1.10. If not, then the functionality for sorting cells by their barcode will not work. If a user runs 10x_PlusCustomBarcode_full_pipeline_perGSM.sh, the pipeline should actually check if the environment has the necessary scripts to run. If not, an error will be thrown saying what is and isn't missing.

Once the relevant packages have been installed a user should be able to modify the script to suit their needs. At the beginning of the 10x_PlusCustomBarcode_full_pipeline_perGSM.sh script there are numerous parameters. Each of these should be commented sufficiently for the user to identify their meaning. The most important starting point is the DATA_DIR variable which should point at a directory containing all of the fastq files from within a GEO GSM. Other things to be wary of include the variable **MitoChromosome**. 

Once these preliminaries have begun the pipeline operates on a 5-step procedure. 

0. Preliminaries - Checks all the relevant installs and helpfully points out errors if there are any and offers solutions
1. Runs alignment with STARSolo and modifying using the parameters from the beginning of the pipeline as inputs.
2. 



We ran the scripts on a PBS cluster scheduler on the Imperial College London. Operability may differ on other cluster types.

