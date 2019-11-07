# library_design_ncRNA -- Kiyan Shabestary 14.03.16
# This script aims to find two most unique sgRNA per ncRNA to create a library. 

########################################## Secondary modules ############################################
# For simplicity for the reader, all modules are incorporated in the same script

import sys
import datetime

def read_dict(fasta_file):

    entry = {}

    fh = open(fasta_file)
    for line in fh.readlines():
        entry[line.strip()] = ''

    fh.close()
    return entry

# Reads fasta file with gene names and sequences.
def read_fasta(fasta_file):
	n = 1
	entry = {}

	fh = open(fasta_file)
	for line in fh.readlines():
		if line[0] == '>':
			ID = 'Entry' + str(n)
			entry[ID] = ''
			n += 1

		else:
			entry[ID] += line.strip()

	fh.close()
	return entry

# Reads fasta file of whole genome. Used for off-target screening.
def read_genome(fasta_file):

	entry = ''
	fh = open(fasta_file)
	for line in fh.readlines():
		if line[0] == '>': continue
		entry = entry + line.strip()

	fh.close()
	return entry

# Motif to be searched is encoded as a class. 
class Motif:
	
	penalties = []
	motif_parts = 1
	lower_interval = 0
	upper_interval = 0
	second_part_start = 0
	c = 0

	def __init__(self, filename):

		fh = open(filename)

		for line in fh.readlines():
			if line[0] == '#' : fh.readline()
			elif line[0] == '*' : 
				self.second_part_start = self.c
				self.lower_interval= int(line.split('\t')[1].split('-')[0].strip())
				self.upper_interval= int(line.split('\t')[1].split('-')[1].strip())
				self.motif_parts += 1
			elif line[0] == 'A' or 'T' or 'G' or 'C' or 'N':
				penalty = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0}
				self.c += 1
				if len(line.split('\t')[0].strip()) == 2 :
					
					for base in penalty :
						penalty[base] = int(line.split('\t')[1].strip())

					penalty[line.split('\t')[0].strip()[0]] = 0
					penalty[line.split('\t')[0].strip()[1]] = 0
					penalty['N'] = 0
					self.penalties.append(penalty)

				else :

					for base in penalty :
						penalty[base] = int(line.split('\t')[1].strip())

					penalty[line.split('\t')[0].strip()] = 0
					penalty['N'] = 0
					self.penalties.append(penalty)
					
			else : raise ValueError('Unknown symbol ' + line[0])

########################################## Main module #################################################

# The main function uses 4 different loops. 1st loop for every entries in the fasta file. 2nd loop
# go through all possible intervals between the two motif parts. 3rd loop to try different
# starting positions for the sequence and finally 4th loop (while loop) to go through every bases
# of the motif sequence until the calculated penalty is above the deviation. This function returns
# a directory of starting position matches (for every entries) and in brackets for each starting position 
# inside the sequence the length between the 1st and 2nd part of the signal (1st number) and the photospacer
# (2nd number)

# For ncRNAs, since the region is small compare to ORFs, no restrictions concerning sgRNAs position is applied
def find_sequence(deviation, motif, entries):
	#print '#find_sequence'
	solutions = {}
	global total_entries
	global solutions_number
	total_entries = 0
	solutions_number = 0

	for entry in entries:
		solution = {}
		total_entries += 1
		for j in xrange(motif.lower_interval, motif.upper_interval+1, 1):
			k = 0
			for start in range(len(entries[entry])-(len(motif.penalties))-j+1):
				i = 0
				Sum = 0
				#if (start == 500 or start == int(0.75*len(entries[entry]))) : break # We do not want our sgRNA to be far from the start

				while (Sum <= deviation and i < len(motif.penalties)):
					if i < motif.second_part_start :
						Sum += motif.penalties[i][entries[entry][start+i]]
					else :
						Sum += motif.penalties[i][entries[entry][start+i+j]]
					i += 1

				if i == (len(motif.penalties)) and Sum <= deviation and not not_good(entries[entry][k:k+j+len(motif.penalties)]):
					Seq = entries[entry][k:k+j+len(motif.penalties)]
					if start in solution.keys():
						solution[start].extend([j,Seq])
					else:
						solution[start] = [j, Seq] #j is the gap
				k += 1
					
		solutions[entry] = solution
		if solutions[entry]:
			solutions_number += 1

	return solutions

# Makes reverse strand
def reverse_strand(fwd_strand):
	rev_strand = ''
	# find complementary 3'-5'
 	for base in fwd_strand:
 		if base == 'X':
 			rev_strand = rev_strand + 'X'
 		if base == 'A':
 			rev_strand = rev_strand + 'T'
 		if base == 'G':
 			rev_strand = rev_strand + 'C'
 		if base == 'C':
 			rev_strand = rev_strand + 'G'
 		if base == 'T':
 			rev_strand = rev_strand + 'A'
 	# translate to 5'-3'
 	return rev_strand[::-1]

# UPDATE: v3 Now takes the sgRNAs with the two best fitness (position, off-target sites)
# Reminder : 
# results[entry] -> all sgRNAs for a given ORF
# sol -> position of the sgRNA in the gene
# results[entry][sol][0] -> interval between the PAM and the T
# results[entry][sol][1] -> sgRNA sequence

# 1) Make a primary list of unique sgRNA
# 2) a. Find min take it if good
#	 b. If none go for the one with the lowest minScore
# 3) a. Find another one not close (>10 bp) and take it if good
#	 b. If none go for the one with the lowest minScore
def get_two_best(results, genome): 
	#print '#get_two_best'
	best_results = {}
	global entries_number
	entries_number = 0

	for entry in results.keys():
		print '>'+str(entry)
		best_result = {}
		unique_sgRNAs = {}
		previous_solution = 0

		if results[entry]: entries_number += 1
		# Compute the fitness values of the sgRNAs (results[entry])

		# 1) Make list of unique sgRNAs. For mulitple sgRNAs for one start site, only the first is taken (smallest)
		for sol in results[entry].keys():
			if len(find_off_target(results[entry][sol][1], genome)) == 1:
				unique_sgRNAs[sol] = results[entry][sol][1]


		# ***** ROUND 1 *****
		# If it is unique, take this
		if len(unique_sgRNAs) != 0:
			closest_position = min(unique_sgRNAs)
			best_result[closest_position] = results[entry][closest_position][1]
			previous_position = closest_position
			unique_sgRNAs.pop(closest_position, None)
			results[entry].pop(closest_position, None)

		else:#lowest minScore
			scores = {}
			for sol in results[entry].keys():
				#print 'Entered minScore'
				print sol
				minScore = 0
				off_targets = find_off_target(results[entry][sol][1], genome)
				for event in range(len(off_targets)):
					if (off_targets[event][1]==0): 
						minScore += 1
				scores[sol] = minScore
			#Put in a list the ones with minScore
			#print 'Scores'
			#print scores
			if len(scores) != 0:
				candidates = {}
				abs_minScore = min(scores.itervalues())
				for sol in scores.keys():
					if int(scores[sol]) == int(abs_minScore): candidates[sol] = results[entry][sol]
				#Take the one closest
				closest_position = min(candidates)
				best_result[closest_position] = results[entry][closest_position][1]
				previous_position = closest_position
				unique_sgRNAs.pop(closest_position, None)
				results[entry].pop(closest_position, None)

		# ***** ROUND 2 *****

		# If it is unique, take this. Now also integrates selection of sgRNAs not +/- 10bp from previous start
		if len(unique_sgRNAs) != 0:
			#print 'entered if len(unique_sgRNAs) not = 0'
			for target in unique_sgRNAs.keys():
				closest_position = min(unique_sgRNAs)

				if not ((closest_position <= previous_position + 10) and (closest_position >= previous_position - 10)):
					best_result[closest_position] = results[entry][closest_position][1]
					break
				elif len(unique_sgRNAs) >= 2:
					unique_sgRNAs.pop(closest_position, None)
					closest_position = min(unique_sgRNAs)
				else: 
					best_result[closest_position] = results[entry][closest_position][1]
					break

		if len(best_result) == 1:#lowest minScore
			#print 'entered case len(best_result) == 1'
			scores = {}
			for sol in results[entry].keys():
				minScore = 0
				off_targets = find_off_target(results[entry][sol][1], genome)

				for event in range(len(off_targets)):
					if (off_targets[event][1]==0): 
						#print 'Check'
						minScore += 1	
				if not ((sol <= previous_position + 10) and (sol >= previous_position - 10)):
					scores[sol] = minScore

			if len(scores) != 0:
				#Put in a list the ones with minScore
				candidates = {}
				abs_minScore = min(scores.itervalues())
				for sol in scores.keys():
					if int(scores[sol]) == int(abs_minScore): 
						candidates[sol] = results[entry][sol]
				#Take the one closest
				closest_position = min(candidates)
				best_result[closest_position] = results[entry][sol][1]

		best_results[entry] = best_result

	return best_results

# This off_target finder look at mismatches in first 15 and and the remaining part. Now, number of mismatches in the SEED to have off-target is included
# off_targets[2n] = Where it binds
# off_targets[2n+1] = Number of mismatch in SEED
def find_off_target(fwd_target, genome):

	off_targets = []

	#rev_target = reverse_strand(fwd_target)
	rev_genome = reverse_strand(genome)
	#print rev_genome

	#   ******  NGG Search  ******

	for start in range(len(genome)-(len(fwd_target))+1):
		i1 = 0
		Sum1 = 0
		Sum2 = 0

		while (Sum1 <= 1 and Sum2 <= 1 and i1 < len(fwd_target)):
			if (genome[start+i1] != fwd_target[i1] and i1 != 2): # If mismatch occurs at position #3, we do not take that into account
				if i1 == 0 or i1 == 1: break # If mismatch occurs in PAM reject offtarget
				elif i1 >= 15: 
					Sum2 += 1
				else : Sum1 += 1 # SEED region
				# print 'mismatch!'

			# Case no mismatch in SEED
			if i1 == 14 and Sum1 == 0: # As long as no mismatch in the SEED, we can already account for off-targets
				i1 = (len(fwd_target)-1) # will ensure we get out of the loop
			i1 += 1

		if i1 == (len(fwd_target)) and Sum1 <= 1 and Sum2 <=1:# Repechage
			Seq = genome[start:start+len(fwd_target)]
			#sys.stdout.write('Fwd target binds @%s %s\n' %(str(start),str(Seq)))
			off_targets.append([start, Sum1])

	# ** for rev_target **

	#Sum1 should be allowed to be >2 if Sum2 = 0. For computational time sake SUm2 should be computed before. 
	#-> Reverse the genome to keep same order

		i2 = 0
		Sum1 = 0
		Sum2 = 0
		#print '** %s '% str(start)

		while (Sum1 <= 1 and Sum2 <= 1 and i2 < len(fwd_target)): 
			if (rev_genome[start+i2] != fwd_target[i2] and i2 != 2): 
				if i2 == 0 or i2 == 1: break # If mismatch occurs in PAM reject offtarget
				elif i2 >= 15: 
					Sum2 += 1
				else : Sum1 += 1 # SEED region
				#print 'mismatch!'

			# Case no mismatch in SEED
			if i2 == 14 and Sum1 == 0: # As long as no mismatch in the SEED, we can already account for off-targets
				i2 = (len(fwd_target)-1) # will ensure we get out of the loop
			i2 += 1

		if i2 == (len(fwd_target)) and Sum1 <= 1 and Sum2 <=1:
			Seq = rev_genome[start:start+len(fwd_target)]
			#sys.stdout.write('Rev target binds @%s %s\n' %(str(start),str(Seq)))
			off_targets.append([len(genome)-start-1, Sum1])

	#   ***************  NAG Search  **************
	# Now the PAM can be CTN (rev: NAG). 

		fwd_target_CTN = 'CTN'+fwd_target[3:len(fwd_target)]
		# print fwd_target_CTN

		i1 = 0
		Sum1 = 0
		Sum2 = 0
		#print '** %s '% str(start)

		while (Sum1 <= 1 and Sum2 <= 2 and i1 < len(fwd_target_CTN)):
			if (genome[start+i1] != fwd_target_CTN[i1] and i1 != 2): # If mismatch occurs at position #3, we do not take that into account
				if i1 == 0 or i1 == 1: break # If mismatch occurs in PAM reject offtarget
				elif i1 >= 15: 
					Sum2 += 1
				else : Sum1 += 1 # SEED region
				# print 'mismatch!'
			# Case one mismatch in each region, we do not need to continue
			if (Sum1 == 1 and Sum2 == 1) or (Sum1 == 0 and Sum2 == 3): break

			i1 += 1

		if i1 == (len(fwd_target_CTN)) and ((Sum1 <= 1 and Sum2 == 0) or (Sum1 == 0 and Sum2 <= 2)):# Repechage
			Seq = genome[start:start+len(fwd_target_CTN)]
			#sys.stdout.write('Fwd target (CTN) binds @%s %s\n' %(str(start),str(Seq)))
			off_targets.append([start,Sum1])

	# ** for rev_target **

	#Sum1 should be allowed to be >2 if Sum2 = 0. For computational time sake Sum2 should be computed before. 
	#-> Reverse the genome to keep same order

		i2 = 0
		Sum1 = 0
		Sum2 = 0
		#print '** %s '% str(start)

		while (Sum1 <= 1 and Sum2 <= 2 and i2 < len(fwd_target_CTN)): 
			if (rev_genome[start+i2] != fwd_target_CTN[i2] and i2 != 2): 
				if i2 == 0 or i2 == 1: break # If mismatch occurs in PAM reject offtarget
				elif i2 >= 15: 
					Sum2 += 1
				else : Sum1 += 1 # SEED region
				#print 'mismatch!'

			# Case no mismatch in SEED
			if (Sum1 == 1 and Sum2 == 1) or (Sum1 == 0 and Sum2 == 3): break

			i2 += 1

		if ((i2 == len(fwd_target_CTN)) and ((Sum1 <= 1 and Sum2 == 0) or (Sum1 == 0 and Sum2 <= 2))):# Repechage
			Seq = rev_genome[start:start+len(fwd_target_CTN)]
			#sys.stdout.write('Rev target (NAG) binds @%s %s\n' %(str(start),str(Seq)))
			off_targets.append([len(genome)-start-1,Sum1])



	return off_targets



#Tests whether EcoRI, XbaI, SpeI, PstI can cut the sequence, or GGGG, TTTT and 25%<GC content<75% at most
def not_good(sequence):

	EcoRI = 'GAATTC'
	XbaI = 'TCTAGA'
	SpeI = 'ACTAGT'
	PstI = 'CTGCAG'

	EcoRI_2 = 'CTTAAG'
	XbaI_2 = 'AGATCT'
	SpeI_2 = 'TGATCA'
	PstI_2 = 'GACGTC'

	GGGGGG = 'GGGGGG'
	TTTT = 'TTTT'

	GC = 0
	position = 0

	for base in sequence:
		if position >=3:
			if str(base) == 'G' or str(base) == 'C': GC += 1
		position += 1

	GC_content = float(GC)/(len(sequence)-3)

	return ((GC_content <= 0.25) or (GC_content >= 0.75) or (TTTT in sequence[3:len(sequence)]) or (GGGGGG in sequence[3:len(sequence)]) or (EcoRI in sequence[3:len(sequence)]) or (XbaI in sequence[3:len(sequence)]) or (SpeI in sequence[3:len(sequence)]) or (PstI in sequence[3:len(sequence)]) or (EcoRI_2 in sequence[3:len(sequence)]) or (XbaI_2 in sequence[3:len(sequence)]) or (SpeI_2 in sequence[3:len(sequence)]) or (PstI_2 in sequence[3:len(sequence)]))


def main():

	# No parallelization since ncRNAs are small subsets

	# 1) Find all sgRNAs
	now = datetime.datetime.now()
	motif = Motif('motif.txt')
	entries = read_fasta('NC_000911.txt')
	genome = read_genome('Synechocystis_whole_genome.txt')


	key_name = 'results/160704_keys.txt'
	gh = open(key_name, 'w')

	for entry in entries.keys():
		gh.write('>%s\n%s\n' % (entry,entries[entry])) # OUTPUT 2 FASTA

	gh.close()

	print 'screening for a lot of sgRNAs!...'
	results = find_sequence(1, motif, entries)# // at this point


	# 2) Find the two best using a fitness function
	print 'finding the two best...'
	best_results = get_two_best(results, genome) 

	# Writing the file

	library_name = 'results/160705_sgRNA_Library_ncRNA.txt' 

	fh = open(library_name, 'w')

	print 'writing file and screening for off-targets...\nfound %s solutions out of %s total entries' %(solutions_number, total_entries)
	fh.write('sgRNA_Library_III KS July 2016 - Report created on : %s \n - %s/%s solutions found \n' %(str(now), solutions_number, total_entries))

	other_match = {}
	compl_targets = {}

	for entry in best_results.keys():
		if best_results[entry]:
			for sol in best_results[entry].keys():
				fh.write('>%s|%s \n%s\n' % (entry,sol, best_results[entry][sol])) # OUTPUT 2 FASTA

	fh.close()

	then = datetime.datetime.now()
	print 'start time: %s \nend time: %s' % (str(now), str(then))
	print 'script completed'

	return 0

def test():
	subentries = read_dict('test_ORF.txt')

	# 1) Find all sgRNAs
	now = datetime.datetime.now()
	motif = Motif('motif.txt')
	entries = read_fasta('NC000911.txt')
	#entries = read_fasta('test.txt')
	genome = read_genome('Synechocystis_whole_genome.txt')
	#genome = read_genome('test_genome.txt')
	print 'screening for a lot of sgRNAs!...'
	results = find_sequence(1, motif, entries, subentries)

	# 2) Find the two best using a fitness function
	print 'finding the two best...'
	best_results = get_two_best(results, genome) 

	return 0



main()
#test()
