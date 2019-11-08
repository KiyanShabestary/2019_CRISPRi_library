# 23.09.19 - Kiyan Shabestary
# This script aims to filter data obtained after droplet sorting.

# PART 1: Set low reads (<32) to 0, compute relative reads and filter out replicate with > 10-fold differences
# PART 2: Feed part1 data to pandas and select top hits

import pandas as pd

def read_data():

	file_ = 'input/data.txt'
	fh=open(file_)

	samples = ['Lib36_gDNA','Lib36_sorted1','Lib36_sorted2','Lib66_sorted1','Lib66_sorted2','Lib66_gDNA']
	dicts_={}

	total = [0,0,0,0,0,0,0]

	fh.readline()

	for line in fh.readlines():
		dict_={}
		sgRNA = line.split('\t')[0]
		i=1
		for sample in samples:
			if int(line.split('\t')[i].strip()) <= 32: dict_[sample]=0
			else: 
				dict_[sample]=line.split('\t')[i].strip()
				total[i-1]+=int(line.split('\t')[i].strip())
			i+=1

		dicts_[sgRNA]=dict_

	fh.close()

	return dicts_, total

def make_relative_data(data,total):

	samples = ['Lib36_gDNA','Lib36_sorted1','Lib36_sorted2','Lib66_sorted1','Lib66_sorted2','Lib66_gDNA']

	for sgRNA in data.keys():
		i=0
		for sample in samples: 
			data[sgRNA][sample] = float(float(data[sgRNA][sample])/total[i])
			i+=1

	return data

# Criteria for filtering out
# Accepted if: (filter1) All 0.1 < replicate1/replicate2 < 10
def filter_(dict_):

	low_treshold=0.1
	high_treshold=10.0

	return (low_treshold<=ratio(float(dict_['Lib36_sorted1']),float(dict_['Lib36_sorted2']))<=high_treshold)*(low_treshold<=ratio(float(dict_['Lib66_sorted1']),float(dict_['Lib66_sorted2']))<=high_treshold)

# To handle 0/0 types of calculations
def ratio(a,b):

	if (a == 0 or b == 0): ratio = 1
	else: ratio = float(a)/float(b)

	return ratio

# Go through all sgRNAs and sort out noisy data thanks to filter() function and print the data
def filter_data(dicts_):

	filtered_dicts = {}

	for sgRNA in dicts_.keys():
		if filter_(dicts_[sgRNA]):
			filtered_dict = {}
			for sample in dicts_[sgRNA].keys(): 
				filtered_dict[sample]=dicts_[sgRNA][sample]
			filtered_dicts[sgRNA]=filtered_dict

	return filtered_dicts

# Keeping all data until the end since duplicates need to be removed and therefore more accurate when all sampels are present
def write_data(data):

	outfile_name='results/filtered_data.txt'
	fh=open(outfile_name,'w')

	fh.write('sgRNA\tLib36_gDNA\tLib36_sorted1\tLib36_sorted2\tLib66_sorted1\tLib66_sorted2\tLib66_gDNA\n')

	for sgRNA in data.keys():
		fh.write(sgRNA+'\t'+str(data[sgRNA]['Lib36_gDNA'])+'\t'+str(data[sgRNA]['Lib36_sorted1'])+'\t'+str(data[sgRNA]['Lib36_sorted2'])+'\t'+str(data[sgRNA]['Lib66_sorted1'])+'\t'+str(data[sgRNA]['Lib66_sorted2'])+'\t'+str(data[sgRNA]['Lib66_gDNA'])+'\n')

	fh.close()
	return 0

# Get top n fraction for sorted samples only, transfer that to a list and finally use the list to check which sgRNA is here and transform that into a dict
# UPDATE: Get top fraction for sorted/gDNA ratio
def get_top_fractions(n):
	df=pd.read_csv('results/filtered_data.txt',delimiter='\t')
	#print(df.head())
	# Removing duplicates
	df=df.drop_duplicates(subset=['Lib36_gDNA', 'Lib36_sorted1','Lib36_sorted2', 'Lib66_sorted1','Lib66_sorted2', 'Lib66_gDNA'], keep=False)

	# Defining new columns as ration sorted/gDNA
	df['Lib36_sorted1_ratio']=df['Lib36_sorted1']/df['Lib36_gDNA']
	df['Lib36_sorted2_ratio']=df['Lib36_sorted2']/df['Lib36_gDNA']
	df['Lib66_sorted1_ratio']=df['Lib66_sorted1']/df['Lib66_gDNA']
	df['Lib66_sorted2_ratio']=df['Lib66_sorted2']/df['Lib66_gDNA']

	# Get top fractions
	df_Lib36_1=df.sort_values(by=['Lib36_sorted1_ratio'], ascending=False)
	df_Lib36_1=df_Lib36_1.sgRNA.head(n)
	list36_1=df_Lib36_1.tolist()

	df_Lib36_2=df.sort_values(by=['Lib36_sorted2_ratio'], ascending=False)
	df_Lib36_2=df_Lib36_2.sgRNA.head(n)
	list36_2=df_Lib36_2.tolist()

	df_Lib66_1=df.sort_values(by=['Lib66_sorted1_ratio'], ascending=False)
	df_Lib66_1=df_Lib66_1.sgRNA.head(n)
	list66_1=df_Lib66_1.tolist()

	df_Lib66_2=df.sort_values(by=['Lib66_sorted2_ratio'], ascending=False)
	df_Lib66_2=df_Lib66_2.sgRNA.head(n)
	list66_2=df_Lib66_2.tolist()

	all_dicts=make_dicts([df_Lib36_1,df_Lib36_2,df_Lib66_1,df_Lib66_2])

	return all_dicts

# Get a table of dict, go through all sgRNAs and write the file
def write_final_data(data,all_sgRNAs, n):

	outfile_name='results/all_samples_counts_TOP'+str(n)+'.txt'
	fh=open(outfile_name,'w')

	fh.write('\tLib36_sorted1\tLib36_sorted2\tLib66_sorted1\tLib66_sorted2\n')

	for sgRNA in all_sgRNAs.keys():
		if ((sgRNA in data[0].keys()) or (sgRNA in data[1].keys()) or (sgRNA in data[2].keys()) or (sgRNA in data[3].keys())):
			fh.write(sgRNA+'\t'+str(get_value(data[0],sgRNA))+'\t'+str(get_value(data[1],sgRNA))+'\t'+str(get_value(data[2],sgRNA))+'\t'+str(get_value(data[3],sgRNA))+'\n')

	return 0

# Returns list of all sgRNAs
def read_sgRNA_ref_file():

	file_ = '../map_reads/input/sgRNA_library.txt'
	fh=open(file_)

	all_sgRNAs={}

	for line in fh.readlines():
		if line[0] == '>': 
			sgRNA_ID = line.split('|')[0].split('>')[1] + '_' + line.split('|')[1]
			all_sgRNAs[sgRNA_ID] = ''

	fh.close()

	return all_sgRNAs

def get_value(sample_dict,sgRNA):
	if sgRNA in sample_dict.keys(): value=1
	else: value=0
	return value

# Transform list to dicts
def make_dicts(all_lists):

	Lib36_sorted1={}
	Lib36_sorted2={}	
	Lib66_sorted1={}	
	Lib66_sorted2={}	

	for sgRNA in all_lists[0]:
		Lib36_sorted1[sgRNA]=1

	for sgRNA in all_lists[1]:
		Lib36_sorted2[sgRNA]=1

	for sgRNA in all_lists[2]:
		Lib66_sorted1[sgRNA]=1

	for sgRNA in all_lists[3]:
		Lib66_sorted2[sgRNA]=1

	return [Lib36_sorted1,Lib36_sorted2,Lib66_sorted1,Lib66_sorted2]

def main():
	#PART 1: Read data
	data,total=read_data() # Set <32 to 0
	total=[19031710,3748295,2754006,3089167,2943335,21579095] # Total number of reads obtained from pandas when all duplicates except one were removed
	data=make_relative_data(data,total) # Make relative data
	data=filter_data(data) # Filter out replicates with more than 10-fold differences

	write_data(data)

	#PART 2: Pandas treatment to get top fractions

	#Reference
	all_sgRNAs=read_sgRNA_ref_file()
	n_s = [10,100,200,250,500,850,1000]

	for n in n_s:
		all_lists = get_top_fractions(n)

		#write file
		write_final_data(all_lists,all_sgRNAs, n)

	return 0

main()
