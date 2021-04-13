#!/usr/bin/env python

##### Code has not been tested on unsorted bam files, sort on barcode (CB):
##### samtools sort -t CB unsorted.bam > sorted_tags.bam
###
##### INPUT: .bam file to be sorted and output directory to place split BC
##### OUTPUT: .bam file for each unique barcode, best to make a new directory

### INPUT and OUTPUT first and second arguments resepctively to be given when you run the script


### Python 3.6.8
import pysam
import sys
### Input varibles to set
# file to split on
unsplit_file = sys.argv[1]
# where to place output files
out_dir = sys.argv[2] + "/"

# variable to hold barcode index
CB_hold = 'unset'
itr = 0
# read in upsplit file and loop reads by line
samfile = pysam.AlignmentFile( unsplit_file, "rb")
for read in samfile.fetch( until_eof=True):
    # barcode itr for current read
    CB_itr = read.get_tag( 'CB')
    # if change in barcode or first line; open new file
    if( CB_itr!=CB_hold or itr==0):
        # close previous split file, only if not first read in file
        if( itr!=0):
            split_file.close()
        CB_hold = CB_itr
        itr+=1
        split_file = pysam.AlignmentFile( out_dir + "CB_{}.bam".format( CB_hold), "wb", template=samfile)

    # write read with same barcode to file
    split_file.write( read)
split_file.close()
samfile.close()
