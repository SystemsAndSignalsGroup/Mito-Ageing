import numpy as np
import pysam
import os
import re
import pandas as pd
import requests
import pickle

def check_read(read):
             """Check if a read passes mapping quality and has none of the following flags set:
		BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL"""
             return ((read.mapping_quality >= min_map_qual)&((read.flag & (512+256+4)) == 0))

def call_insertion_deletion(bam_file, ref, mito_name, min_base_qual, min_map_qual, min_hf, min_depth):
	"""
	Takes in an aligned bam file and searches for anywhere which is marked as an 
	unexpected deletion or insertion and creating two dataframes logging all of them
	It also creates a list of the read depth of every position with an aligned read by first finding

	a max depth which allows all reads to be considered. It also returns a dataframe of the total depth, mean
	depth and number of bases passing filter controls.

	The flag_filter excludes all reads for which the flags BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL are set.
	Parameters
	-------------
	bam_file - the file under consideration
	ref - list containing the reference genome used for alignment
	mito_name - reference name of mitochondrial reference in question
	min_base_qual - the minimum quality requirement for a base to be called (Phred scored)
	min_map_qual - the minimum quality requirement for a mapped read to be used (Phred scored)
	min_hf - the minimum heteroplasmy for which a read will be called
	min_depth - the minimum number of reads required at a base for a deletion to be considered
	"""
	ins_df = pd.DataFrame(columns=['POS','REF','ALT','DF','HF'])
	del_df = pd.DataFrame(columns=['POS','REF','ALT','DF','HF'])
	depth = {}
	depth_df = pd.DataFrame(columns=['total_depth','mean_depth','bases_passed'])
	next_max = 10000

	while True:
		current_max=next_max
		iter = bam_file.pileup(mito_name,stepper='nofilter',flag_filter=512+256+4, max_depth=current_max)
		depth1=[]
		for pileupcolumn in iter:
			depth1.append(pileupcolumn.n)
		next_max = current_max+10000
		iter = bam_file.pileup(mito_name,stepper='nofilter',flag_filter=512+256+4, max_depth=next_max)
		depth2=[]
		for pileupcolumn in iter:
			depth2.append(pileupcolumn.n)
		try:
			if max(depth1)==max(depth2):
				break
		except:
			depther = [0,0,0]
			depth_df.loc[len(depth_df)]=depther
			return depth, ins_df, del_df, depth_df

	passing=0
	iter = bam_file.pileup(mito_name,stepper='all',flag_filter=512+256+4,ignore_overlaps=False,min_base_quality=min_base_qual,min_mapping_quality=min_map_qual, max_depth=current_max)

	for pileupcolumn in iter:
		position = pileupcolumn.pos + 1
		delete = [s for s in pileupcolumn.get_query_sequences(mark_matches=False,add_indels=True) if "-" in s]
		insert = [s for s in pileupcolumn.get_query_sequences(mark_matches=False,add_indels=True) if "+" in s]
		As = [s for s in pileupcolumn.get_query_sequences(mark_matches=False,add_indels=True) if s.upper() in ['A','C','G','T']]
		depth[position] = len(As)+len(delete)+len(insert)
		if depth[position] > min_depth:
			passing += 1
			if len(insert)>depth[position]*min_hf:
				q = set(insert)
				l = {}
				for k in q:
					if k.upper() in l.keys():
						l[k.upper()] = l[k.upper()] + insert.count(k)
					else:
						l[k.upper()] = insert.count(k)
				for j in l:
					mids = []
					refer = j.split('+')[0]
					nexter = ''.join([c for c in j.split('+')[1] if not c.isdigit()])
					mids.extend([position,refer,refer+nexter,depth[position],l[j]/depth[position]]) # convert from 0-REF to 1-REF
					ins_df.loc[len(ins_df)] = mids

			if len(delete)>depth[position]*min_hf:
				q = set(delete)
				l = {}
				for k in q:
					if k.upper() in l.keys():
						l[k.upper()] = l[k.upper()] + delete.count(k)
					else:
						l[k.upper()] = delete.count(k)
				for j in l:
					mids = []
					refer = j.split('-')[0]
					digit = ''.join([c for c in j.split('-')[1] if c.isdigit()])
					nexter = ''.join(ref[(position):(position+int(digit))])
					mids.extend([position,refer+nexter,refer,depth[position],l[j]/depth[position]]) # to change from 0-REF to 1-REF
					del_df.loc[len(del_df)] = mids

	mean_depth=sum(depth.values())/len(ref)
	depther = [bam_file.count(mito_name),mean_depth,passing]
	depth_df.loc[len(depth_df)]=depther
	ins_df = ins_df[ins_df['HF'] > min_hf]
	del_df = del_df[del_df['HF'] > min_hf]

	return depth, ins_df, del_df, depth_df


def call_snps(bam_file, depth, mito_name, min_base_qual, min_map_qual, min_depth):
	"""
	Takes in aligned bam file and calls all snps subject to desired quality control thresholds

	Parameters
	-------------	
	bam_file - the file under consideration
	depth - dictionary of read depth at every covered position
	mito_name - reference name of chromosome in question
	min_base_qual - the minimum quality requirement for a base to be called (Phred scored)
	min_map_qual - the minimum quality requirement for a mapped read to be used (Phred scored)
	min_hf - the minimum heteroplasmy for which a read will be called
	min_depth - the minimum number of reads required at a base for a deletion to be considered

	"""

	# Now find reads aligning to each base at every position
	t = bam_file.count_coverage(mito_name,read_callback=check_read,quality_threshold=min_base_qual)
	snps_df = pd.DataFrame(t).T.rename(columns={0:'A',1:'C',2:'G',3:'T'})
	snps_df['REF'] = ref
	snps_df.reset_index(inplace=True, drop=False)
	snps_df = snps_df.rename(columns={'index':'POS'})
	snps_df['POS'] = snps_df['POS'] + 1
	snps_df['DF'] = snps_df['POS'].map(depth)
	snps_df = snps_df[snps_df.DF>min_depth].reset_index(drop=True)

	# Look at heteroplasmy of each base if it is not the reference base and join together
	bases = ['A','T','C','G']
	splot = {}
	for base in bases:
		splot[base] = pd.melt(snps_df,id_vars = ['POS'],value_vars=[base])
		splot[base]['REF'] = snps_df['REF']
		splot[base]['DF'] = snps_df['DF']
		splot[base]['HF'] = splot[base]['value']/snps_df['DF']
		splot[base] = splot[base].rename(columns={'variable':'ALT'})
		splot[base] = splot[base][['POS','REF','ALT','DF','HF']]
		splot[base] = splot[base][(splot[base]['REF'] != splot[base]['ALT']) & (splot[base]['REF'] != 'N') & (splot[base]['HF'] > min_hf)]
	snps_df = pd.concat(splot).reset_index(drop=True)

	return snps_df


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
	parser.add_argument('-sample_id', default=None, type=str,
		help="Desired identifier for output file")

	parser.add_argument('-min_base_qual', default=30, type=int,
		help="Minimum base quality for inclusion in counts")
	parser.add_argument('-min_map_qual', default=30, type=int,
		help="Minimum read alignment quality for inclusion in counts")

	parser.add_argument('-min_depth', default=50, type=int,
		help="Minimum read depth to call a heteroplasmy")
	parser.add_argument('-min_hf', default=0.05, type=float,
		help="Minimum heteroplasmy frequency to call a heteroplasmy")
	parser.add_argument('-mito_name', default='chrM', type=str,
		help="Name of the reference sequence used for alignment")

	args = parser.parse_args()

	bam_file_name = args.bam_file_name
	mito_file = args.mito_file
	sample_id = args.sample_id

	min_base_qual = args.min_base_qual
	min_map_qual = args.min_map_qual
	min_depth = args.min_depth
	min_hf = args.min_hf

	mito_name = args.mito_name

	reffile = pysam.FastaFile(mito_file)
	ref = list(reffile.fetch(mito_name))
	reffile.close()

	bam_file = pysam.AlignmentFile(bam_file_name, "rb")

	depth, ins_df, del_df, depth_df = call_insertion_deletion(bam_file, ref, mito_name,
						 min_base_qual, min_map_qual, min_hf, min_depth)
	
	depth_df['sample_id'] = sample_id
	depth_df.to_csv('{}_depth.csv'.format(sample_id),index=False)

	f = open("{}_coverage.pkl".format(sample_id),"wb")
	pickle.dump(depth,f)
	f.close()
	
	if depth_df['bases_passed'][0] != 0:
		snps_df = call_snps(bam_file, depth, mito_name, min_base_qual, min_map_qual, min_depth)

		mut_df = pd.concat([snps_df,ins_df,del_df]).sort_values('POS').reset_index(drop=True)
		mut_df['sample_id'] = sample_id

		mut_df.to_csv('{}_variants.csv'.format(sample_id),index=False)
