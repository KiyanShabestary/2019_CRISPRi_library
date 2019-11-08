#!usr/bin/ncbi-blast-2.6.0+

# map_reads -- Kiyan Shabestary 26.02.18
# This script intends to map trimmed and filtered NGS reads to our sgRNA library. 

import sys
import datetime

# This is the main function. This function goes through all reads and map them back to the reference using an align function.
def map_reads(library,reads_file, duplicate_entry_file):
	perfect_counts = {} #100% match no embiguity		24 + 20 (max 25) + 27(out of 30) = 71  76 max seq data
	read_number = 0

	fh = open(reads_file)
	gh = open(duplicate_entry_file,'w')

	for line in fh.readlines():
		if line[0] == 'C':
			read_number +=1

			line=line.strip()
			if 'N' in line: continue #No embiguities
			
			alignment_results = align(line,library) #We will modify the length of the line with the length of the sgRNA in the loop
			
			#Handle specific case when more than 2 alignments
			for alignment in alignment_results:
				if alignment not in perfect_counts.keys(): perfect_counts[alignment] =1
				else : perfect_counts[alignment] += 1

				if len(alignment_results) >= 2:
					gh.write(str(alignment)+'\t') 

			#To separate double
			if len(alignment_results) >= 2:
				gh.write('\n****\n')

	fh.close()
	gh.close()

	return perfect_counts

# Align to every entry in ths sgRNA library
def align(query_string,library):
	# Convert the sgRNA photospacer string in a fasta file
	alignment_results=[]

	for entry in library.keys():
		if len(library[entry]) >= len(query_string):
			if query_string in library[entry]: alignment_results.append(str(entry))
		else: #Case where some bp non specific to library in sequencing query
			if library[entry] == query_string[0:(len(library[entry]))]: alignment_results.append(str(entry))

	return alignment_results

def print_results(alignment_results, file_to_write):
	fh = open(file_to_write,'w')
	now = datetime.datetime.now()
	fh.write('Library sequencing results - KS %s \n ******************\n Displaying perfect match reads \n ****************** \n'%(str(now)))

	for alignment in alignment_results:
		fh.write('>'+str(alignment)+'\t'+str(alignment_results[alignment])+'\n')

	fh.close()
	return 0

# Read library
def read_library(fasta_file):

	entry = {}
	fh = open(fasta_file)
	for line in fh.readlines():
		if line[0] == '>':
			ID = line.split('|')[0]+'_'+line.split('|')[1]
			entry[ID.split('>')[1]] = ''

		else:
			entry[ID.split('>')[1]] += line.strip()

	fh.close()
	return entry

def main():

	results_library_file = 'results/counts.txt' 				#Output: Counts
	seq_data_file = 'input/test.filtered.fastq'					#Input: Reads in fasta format
	duplicate_ORF_file = 'results/duplicate_sgRNAs.txt' 		#Output: Get a list of reads that align to  multiple sgRNAs
	library = read_library('input/sgRNA_library.txt')			#Input: Library file for referencing

	alignment_results=map_reads(library,seq_data_file,duplicate_ORF_file)
	print_results(alignment_results, results_library_file)
	return 0

main()
#test()
