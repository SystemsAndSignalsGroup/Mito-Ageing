import numpy as np
import pysam
import os
import re
import pandas as pd
import pickle

def check_read(read):
    """Check if a read passes mapping quality and has none of the following flags set:
        BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP"""
    return (((read.flag & (1024+512+256+4)) == 0))		 

def filter_indels(pileupcolumn, min_base_qual, flank = 5):
	"""
	Takes a pileup column and filters for any indels that are at least flank away from the read ends and that have a median quality in flank of min_base_qual

	Parameters
	-------------
	pileupcolumn - a pysam pileup column
	min_base_qual - the minimum base quality required for flanking reads
	flank - the required length of the flanking region for indels
	"""    
	sequences = pileupcolumn.get_query_sequences(mark_matches=False,add_indels=True)
	
	index_ins = []
	index_del = []
	indel_len = []

	for index, s in enumerate(sequences):
		indel_len.append(0)
		if "-" in s:
			index_del.append(index)
			indel_len.append(0)
		elif "+" in s:
			index_ins.append(index)
			indel_len.append(int(''.join([c for c in s.split('+')[1] if c.isdigit()])))

	if (len(index_del)==0) & (len(index_ins)==0):
		return [], []
	
	
	Filtered_Inserts = []
	Filtered_Deletes = []

	for index, read in enumerate(pileupcolumn.pileups):
		if sequences[index] in '<>*':
			pass
		elif (read.query_position > flank) & ((read.query_position + indel_len[index]) < (read.alignment.query_alignment_length - flank)) & (index in index_del):
			upstream_qual = read.alignment.query_qualities[read.query_position-flank: read.query_position]
			downstream_qual = read.alignment.query_qualities[read.query_position+indel_len[index]+1: read.query_position+indel_len[index]+1+flank]
			if (np.median(upstream_qual)>min_base_qual) & (np.median(downstream_qual)>min_base_qual):
				Filtered_Deletes.append(sequences[index])
				
		elif (read.query_position > flank) & ((read.query_position + indel_len[index]) < (read.alignment.query_alignment_length - flank)) & (index in index_ins):
			upstream_qual = read.alignment.query_qualities[read.query_position-flank: read.query_position]
			downstream_qual = read.alignment.query_qualities[read.query_position+indel_len[index]+1: read.query_position+indel_len[index]+1+flank]
			if (np.median(upstream_qual)>min_base_qual) & (np.median(downstream_qual)>min_base_qual):
				Filtered_Inserts.append(sequences[index])
				
	return Filtered_Deletes, Filtered_Inserts

def call_mutations(bam_file, ref, mito_name, min_base_qual, min_map_qual, min_hf, min_depth, flank):
	"""
	Takes in an aligned bam file and searches for anywhere which is marked as an 
	unexpected deletion or insertion and creating two dataframes logging all of them
	It also creates a list of the read depth of every position with an aligned read by first finding
	a max depth which allows all reads to be considered. It also returns a dataframe of the total depth, mean
	depth and number of bases passing filter controls.

	Parameters
	-------------
	bam_file - the file under consideration
	ref - list containing the reference genome used for alignment
	mito_name - reference name of mitochondrial reference in question
	min_base_qual - the minimum quality requirement for a base to be called (Phred scored)
	min_map_qual - the minimum quality requirement for a mapped read to be used (Phred scored)
	min_hf - the minimum heteroplasmy for which a read will be called
	min_depth - the minimum number of reads required at a base for a deletion to be considered
	flank - the required length of the flanking region for indels
	"""
	
	depth_df = pd.DataFrame(columns=['total_depth'])
	mut_df = pd.DataFrame(columns=['POS','REF','ALT','DF_F','DF_R','DF','HF_F','HF_R','HF'])

	Forward_Bases = ['A','C','G','T']
	Reverse_Bases = ['a','c','g','t']

	next_max = 10000
	covDict = {}
    
	while True:
		current_max=next_max
		iter = bam_file.pileup(mito_name,stepper='nofilter',flag_filter=512+256+4, max_depth=current_max)
		depth1=[0]
		for pileupcolumn in iter:
			depth1.append(pileupcolumn.n)
		next_max = current_max+10000
		iter = bam_file.pileup(mito_name,stepper='nofilter',flag_filter=512+256+4, max_depth=next_max)
		depth2=[0]
		for pileupcolumn in iter:
			depth2.append(pileupcolumn.n)
		try:
			if max(depth1)==max(depth2):
				break
		except:
			depther = [0,0,0]
			depth_df.loc[len(depth_df)]=depther
			return covDict, depth_df, mut_df    
    
	passing=0
	iter = bam_file.pileup(mito_name,stepper='samtools',flag_filter=1024+512+256+4,ignore_overlaps=True,
						   min_base_quality=min_base_qual,min_mapping_quality=min_map_qual, max_depth=current_max)

	for pileupcolumn in iter:

		position = pileupcolumn.pos + 1
		sequences = pileupcolumn.get_query_sequences(mark_matches=False,add_indels=True)
	 
		Filtered_Deletes, Filtered_Inserts = filter_indels(pileupcolumn, min_base_qual, flank)

		All_Deletes = [s.upper() for s in Filtered_Deletes]
		Forward_Deletes = [s for s in Filtered_Deletes if s[0] in Forward_Bases]
		Reverse_Deletes = [s for s in Filtered_Deletes if s[0] in Reverse_Bases]

		All_Inserts = [s.upper() for s in Filtered_Inserts]
		Forward_Inserts = [s for s in Filtered_Inserts if s[0] in Forward_Bases]
		Reverse_Inserts = [s for s in Filtered_Inserts if s[0] in Reverse_Bases]

		All_Muts = [s.upper() for s in sequences if s.upper() in Forward_Bases]
		Forward_Muts = [s for s in sequences if s in Forward_Bases]
		Reverse_Muts = [s for s in sequences if s in Reverse_Bases]

		Total_Depth = len(All_Muts)+len(All_Deletes)+len(All_Inserts)
		Forward_Depth = len(Forward_Muts)+len(Forward_Deletes)+len(Forward_Inserts)
		Reverse_Depth =len(Reverse_Muts)+len(Reverse_Deletes)+len(Reverse_Inserts)

		covDict[position] = Total_Depth

		if Total_Depth > min_depth:
			passing += 1

			if All_Muts.count(ref[position-1]) > Total_Depth*(1-min_hf):
				pass

			else:
				insertions = set(All_Inserts)
				
				for number, k in enumerate(insertions):
					if All_Inserts.count(k)>Total_Depth*min_hf:

						refer = k.split('+')[0].upper()
						nexter = ''.join([c for c in k.split('+')[1] if not c.isdigit()])

						try:
							hf_f = Forward_Inserts.count(k)/Forward_Depth
						except:
							hf_f = 0

						try:
							hf_r = Reverse_Inserts.count(k.lower())/Reverse_Depth
						except:
							hf_r = 0

						hf = All_Inserts.count(k)/Total_Depth

						if hf>min_hf:
							mids = [position, refer, refer+nexter, Forward_Depth, Reverse_Depth, Total_Depth, hf_f, hf_r, hf]
							mut_df.loc[len(mut_df)] = mids              


				deletions = set(All_Deletes)

				for number, k in enumerate(deletions):
					if All_Deletes.count(k)>Total_Depth*min_hf:
						
						refer = k.split('-')[0].upper()
						digit = ''.join([c for c in k.split('-')[1] if c.isdigit()])
						nexter = ''.join(ref[(position):(position+int(digit))])

						try:
							hf_f = Forward_Deletes.count(k)/Forward_Depth
						except:
							hf_f = 0

						try:
							hf_r = Reverse_Deletes.count(k.lower())/Reverse_Depth
						except:
							hf_r = 0

						hf = All_Deletes.count(k)/Total_Depth
						
						if hf>min_hf:
							mids = [position, refer+nexter, refer, Forward_Depth, Reverse_Depth, Total_Depth, hf_f, hf_r, hf]
							mut_df.loc[len(mut_df)] = mids

				mutations = set(All_Muts)
				mutations.discard(ref[position-1])

				for k in mutations:
					if All_Muts.count(k) > Total_Depth*min_hf:

						refer = ref[position-1].upper()
						alt = k

						try:
							hf_f = Forward_Muts.count(k)/Forward_Depth
						except:
							hf_f = 0

						try:
							hf_r = Reverse_Muts.count(k.lower())/Reverse_Depth
						except:
							hf_r = 0

						hf = All_Muts.count(k)/Total_Depth                    

						mids = [position, refer, alt, Forward_Depth, Reverse_Depth, Total_Depth, hf_f, hf_r, hf]
						mut_df.loc[len(mut_df)] = mids

	depther = [bam_file.count(mito_name, read_callback = check_read)]
	depth_df.loc[len(depth_df)]=depther
	return covDict, depth_df, mut_df

if __name__ == '__main__':
	import argparse
	from time import time
	import warnings

	warnings.simplefilter("ignore")
	parser = argparse.ArgumentParser(
		description="Compute site read counts for mtDNA given BAM file")
	parser.add_argument('-bam_file_name', default=None, type=str,
		help="Path to bam file needed for variant calling")
	parser.add_argument('-mito_file', default=None, type=str,
		help="Path to fastq file with reference mito genome")    
	parser.add_argument('-cell_id', default=None, type=str,
		help="Unique identifier for of the cell")
	parser.add_argument('-donor_id', default=None, type=str,
		help="Identifier for donor")

	parser.add_argument('-min_base_qual', default=30, type=int,
		help="Minimum base quality for inclusion in counts")
	parser.add_argument('-min_map_qual', default=30, type=int,
		help="Minimum read alignment quality for inclusion in counts")
	parser.add_argument('-flank', default=5, type=int,
		help="Required flanking length for indels")

	parser.add_argument('-min_depth', default=50, type=int,
		help="Minimum read depth to call a heteroplasmy")
	parser.add_argument('-min_hf', default=0.05, type=float,
		help="Minimum heteroplasmy frequency to call a heteroplasmy")
	parser.add_argument('-mito_name', default='chrM', type=str,
		help="Name of the reference sequence used for alignment")

	args = parser.parse_args()

	bam_file_name = args.bam_file_name
	mito_file = args.mito_file
	cell_id = str(args.cell_id)
	donor_id = str(args.donor_id)

	min_base_qual = args.min_base_qual
	min_map_qual = args.min_map_qual
	flank = args.flank
	min_depth = args.min_depth
	min_hf = args.min_hf

	mito_name = args.mito_name

	reffile = pysam.FastaFile(mito_file)
	ref = list(reffile.fetch(mito_name))
	reffile.close()

	bam_file = pysam.AlignmentFile(bam_file_name, "rb")

	covDict, depth_df, mut_df = call_mutations(bam_file, ref, mito_name,
						 min_base_qual, min_map_qual, min_hf, min_depth, flank)
	
	depth_df['cell_id'] = cell_id
	depth_df['donor_id'] = donor_id
	depth_df.to_csv('{}_depth.csv'.format(cell_id),index=False)

	cov_df = pd.DataFrame(index=list(range(1, len(ref)+1)), 
			columns = pd.MultiIndex.from_product([[donor_id],[cell_id]] ,  
						names=['donor_id','cell_id']))
	cov_df[donor_id] = cov_df.index.map(covDict)
	cov_df.to_parquet('{}_coverage.pq'.format(cell_id))
	
		
	mut_df['cell_id'] = cell_id
	mut_df['donor_id'] = donor_id
	mut_df.to_csv('{}_variants.csv'.format(cell_id),index=False)
