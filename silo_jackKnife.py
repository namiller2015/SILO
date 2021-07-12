'''
12/13/2017
Author: Nick Miller
#Silo in python
'''
import sys
import json
import datetime
import shutil, os
def find_CN_changes(training_set, test_set):
	#use the MPD and SD from the training set to check for Copy Number changes in the test set
	#avg the amplicons together to get total gene amplification
	
	for gene in test_set['gene_name']:
		for amplicon_dict in test_set['gene_name'][gene]:
			for amplicon, value in amplicon_dict.items():
				
				MPD = lookup_value(training_set, gene, amplicon, 'Mean_prop_depth')
				SD = lookup_value(training_set, gene, amplicon, 'SD_of_Prop_depth')
				silo_mean_CN = lookup_value(training_set, gene, amplicon, 'silo_mean')
				silo_SD = lookup_value(training_set, gene, amplicon, 'silo_SD')

				upper = MPD + (3 * SD)
				lower = MPD - (3 * SD)
				porp_depths = value[1]['Prop_depths']
				CN_list = []

				for single_porp in porp_depths:
					
					if single_porp == 'FAILED':
						CN_list.append('FAILED')
					else:
						if single_porp < lower or single_porp > upper:
							CN = single_porp/MPD
						else:
							CN = 1.0
					#record this CN value for the test set so we can average across amplicons to get final gene amplification
						CN_list.append(CN)
				value.append({'Copy_number':CN_list})
				value.append({'Mean_prop_depth':MPD})
				value.append({'SD_of_Prop_depth':SD})
				value.append({'silo_mean': silo_mean_CN})
				value.append({'silo_SD': silo_SD})
	
	return test_set


def gene_amplification(test_set):
	
	#determine total gene amplification by taking the mean of the amplicons

	copy_number_dict = {}
	samples_list= []
	for index, sample_name in enumerate(test_set['samples']):
		samples_list.append(sample_name)
		for gene in test_set['gene_name']:

			amp_CN_total = 0
			MPD_total = 0
			SD_total = 0
			silo_mean_total = 0
			silo_SD_total = 0
			for amplicon_dict in test_set['gene_name'][gene]:
				for amplicon, values in amplicon_dict.items():
					
					if values[2]['Copy_number'][index] == 'FAILED':
						pass
					else:

						amp_CN_total += values[2]['Copy_number'][index]
						MPD_total    += values[3]['Mean_prop_depth']
						SD_total     += values[4]['SD_of_Prop_depth']
						silo_mean_total += values[5]['silo_mean']
						silo_SD_total += values[6]['silo_SD']

			amplicon_number = len(test_set['gene_name'][gene])

			if amp_CN_total == 0:
				gene_CN = 'SILO FAILED. RAW READS BELOW 20'
			else:
				gene_CN = round(float(amp_CN_total/amplicon_number), 1)
			#avg_MPD = round(float(MPD_total/amplicon_number), 2)
			#SD_total = round(float(SD_total/amplicon_number), 2)
			
			silo_mean_display = round(float(silo_mean_total/amplicon_number), 2)
			silo_SD_display = round(float(silo_SD_total/amplicon_number), 2)

			if gene not in copy_number_dict:
				copy_number_dict[gene] = [amplicon_number, silo_mean_display, silo_SD_display, gene_CN]
			else:
				copy_number_dict[gene].append(gene_CN)
	
	return samples_list, copy_number_dict


def calc_MPD_SD(RAW_training_set_obj):
	#Calculate the MPD and SD
	for gene in RAW_training_set_obj['gene_name']:
		for amplicon_dict in RAW_training_set_obj['gene_name'][gene]:
			for amplicon, value in amplicon_dict.items():
				
				porp_depths = value[1]['Prop_depths']

				#remove outliers
				cleaned_porp_depths = remove_outliers(10, porp_depths)

				#get MPD and SD from the cleaned porp depths
				MPD = average(cleaned_porp_depths)
				SD = calc_SD(cleaned_porp_depths)
				
				if MPD == 0.0:
					print(gene)
					print(amplicon)
					print(MPD)
					print(porp_depths)
					print(cleaned_porp_depths)

				#add MPD and SD to RAW_training_set_obj and finished_training_set_obj
				value.append({'Mean_prop_depth':MPD})
				value.append({'SD_of_Prop_depth':SD})
				
				CN_list = []
				for single_porp in porp_depths:
					CN = single_porp/MPD
					CN_list.append(CN)
				
				silo_meanCN_cleaned = remove_outliers(5, CN_list)
				#store the silomeanCN which we will avg at the gene level to get the silo_mean to be displayed on the front end
				silo_meanCN = average(silo_meanCN_cleaned)
				value.append({'silo_mean': silo_meanCN})
				
				silo_SD_of_CNs = calc_SD(silo_meanCN_cleaned)
				value.append({'silo_SD': silo_SD_of_CNs})

			
	return RAW_training_set_obj


def calc_porp_depth(data_set_obj):
	#Calculate the Porportional depth for each target in a sample
	
	for gene in data_set_obj['gene_name']:					#iterates over the keys in gene_name dict which is genes
		for amplicon_dict in data_set_obj['gene_name'][gene]:		#iterates over the keys in gene dict which is amplicon dict
			for amplicon, value in amplicon_dict.items():		#iterates over amplicon name and the list of sub amplicon dictionaries 
				porp_depth_list = []
				for index, raw_reads in enumerate(value[0]['Raw_Depths']):	#iterates over the list in raw_depths for each amplicon dict and gives index
					avg_sample_depth = data_set_obj['Ave_sample_depths'][index]		#get the avg_sample_depth using the index position of the raw_depth

					if avg_sample_depth == 'FAILED':
						porp_depth_list.append('FAILED')
					else:
					
						porp_depth = raw_reads/avg_sample_depth						#Calc porporitional depth
						porp_depth_list.append(porp_depth)
			
				value.append({'Prop_depths':porp_depth_list})
	return data_set_obj

def calc_avg_sample_depth(sample_raw_reads_obj, data_set_obj):
	#calculate the avg depth on a sample basis and store in Raw_training_set_obj
	
	avg_sample_depth = []
	for sample_name in data_set_obj['samples']:
		#need to remove outliers 
		cleaned_values = remove_outliers(10, sample_raw_reads_obj[sample_name])
		sample_avg = average(cleaned_values)
		
		if sample_avg == 0 or sample_avg == 0.0:
			sample_avg = 'FAILED'

		avg_sample_depth.append(sample_avg)

	data_set_obj['Ave_sample_depths'] = avg_sample_depth
	return data_set_obj

def remove_outliers(ignPercent, data):
	#remove the outliying percentage from data list based on passed in value
	data = sorted(data)
	n = len(data)

	outliers = int((n * ignPercent) / 100)
	trimmed_data = data[outliers: n-outliers]

	return trimmed_data

def lookup_value(data_obj, gene, amplicon, desired_value):
	#abstracted function to make it easier to get values out of the large data_obj
	gene_dict = data_obj['gene_name'][gene]
	for amplicon_dict_list in gene_dict:
		if amplicon in amplicon_dict_list:
			for dict_lists in amplicon_dict_list[amplicon]:
				if desired_value in dict_lists:
					return dict_lists[desired_value]
	
	print('Could not find:')
	print(gene)
	print(amplicon)
	print(desired_value)
	sys.exit('\n EXIT Lookup value function')


##############################################################
# Fucntions taken from statistics module in Python 3.4
# did not import from module so as to not create dependency and versions problems
def average(values):
	#return the average of the list
	length = len(values)
	if length < 1:
		sys.exit('Error: No values in list')
	return sum(values)/length

def calc_SD(values):
	#return the Standard diviation from a list of values
	#Uses Population SD NEVER sample SD
	n = len(values)
	if n < 2:
		sys.exit('Error: variance requires at least 2 values')
	square_sum = squareSum(values)
	var = square_sum/n

	return var**0.5

def squareSum(values):
	#calulate the sum of square diviations of sequence data
	c = average(values)
	square_sum = sum((x-c)**2 for x in values)

	return square_sum

def pretty_print(samples_list, CN_results, filename, analysis_home_dirs):
	#print out results in an easy to view format
	
	epoch = str(datetime.datetime.today().timestamp()).split(".")[0]
	
	base_file_name = str(filename.split('/')[-1:][0])
	#file_handle = open('/data/software/flypeC/silo/'+str(base_file_name[0])+'_'+str(epoch), 'w')
	#file_handle.write('gene	num_amplicons\tsilo_mean\tsilo_SD\t' + '\t'.join(samples_list) + '\n')
	
	#output_file = '/data/software/flypeC/silo/'+str(base_file_name[0])+'_'+str(epoch)
	
	#iterate through analysis home dirs and write output file to each
	for home_dir in analysis_home_dirs:
		print(home_dir)
		work_accession = os.path.basename(os.path.normpath(str(home_dir)))
		results_file = home_dir+'_SILO_RESULTS_OUTPUT_'+str(epoch)+'.txt'
		results_file_handle = open(results_file, 'w')

		results_file_handle.write('gene\tnum_amplicons\tsilo_mean\tsilo_SD\t' + '\t'.join(samples_list) + '\n')
		
		#write to the file
		for key,value in sorted(CN_results.items()):
			results_file_handle.write(key+'\t'+'\t'.join(map(str,value)) + '\n')
		
		results_file_handle.close()
		output_file = home_dir+'SILO_RESULTS_OUTPUT_'+work_accession+'_'+str(epoch)


	#for key,value in sorted(CN_results.items()):
	#	file_handle.write(key+'\t'+'\t'.join(map(str,value)) + '\n')

	#file_handle.close()

	return output_file


#input and building functions
################################################################
def read_matrixFile(bcmatrix):

	#read in bcmatrix file. This will be created directly from ion torrent. Or if user is providing bam files we will parse bams into bcmatrix format.
	#Prase this file into Base data Object used by both training set and test set
	# also need to create simple_raw_reads obj for vertical iteration to generate the avg_sample_depths
	file_handle = open(bcmatrix, 'r')

	#initalize dicts and lists needed for full data obj
	RAW_data_set_obj = {}
	samples = []
	gene_name = {}
	
	#initalize dict for simple_raw_reads obj
	simple_raw_reads_obj = {}		# sample_raw_reads_obj = {'sample_name': [raw_read_values_by_sample]

	first_line = file_handle.readline()
	first_line = first_line.strip('\n').split('\t')[2:] #slice the first two out

	for sample_name in first_line:
		samples.append(sample_name)
		simple_raw_reads_obj[sample_name] = []

	#iterate over lines in file. First line contains sample names
	for line in file_handle.readlines():
			# each line is a new amplicon so we need to reinit the sub amplicon dicts
			amplicons_dict = {}
			sub_amplicon_dict = {}
			line = line.strip('\n').split('\t')
		
			gene = line[0]
			amplicon = line[1]
			reads = [int(float(x)) for x in  line[2:]]

			sub_amplicon_dict['Raw_Depths'] = reads
			#build amplicon dictionary
			amplicons_dict[amplicon] = [sub_amplicon_dict]
		

			if gene not in gene_name:
				gene_name[gene] = [amplicons_dict]
			else:
				gene_name[gene].append(amplicons_dict)

			#build out simple raw reads obj
			simple_raw_reads_obj = build_simpleRawObj(samples, reads, simple_raw_reads_obj)
	
	#build data obj
	RAW_data_set_obj['samples'] = samples
	RAW_data_set_obj['gene_name'] = gene_name
	return RAW_data_set_obj, simple_raw_reads_obj

def build_simpleRawObj(samples_list, rawReads_list, simpleObj):

	for index, sample_name in enumerate(samples_list):
		simpleObj[sample_name].append(rawReads_list[index])
	
	return simpleObj


def read_bed(bedfile):
	#read in bedfile and parse it
	#Need to get gene name, amplicon name, amplicon region start, amplicon region end

	file_handle = open(bedfile, 'r')

	#skip first line of bed file because its the header line
	lines = file_handle.readlines()[1:]
	
	bed_info = []

	#iterate through
	for line in lines:
		#tab delimited
		line = line.split('\t')

		start = line[1]
		stop = line[2]
		amp = line[3]
		gene = line[7].split('=')

		bed_info.append((start,stop,amp,gene))

	return bed_info


def read_bam(bed_info, bamfiles):
	#bed info is a list of tuples contatining startCoordinate, stopCoordinate, amplicon_name, gene_name
	#function will use the bed info to get the read depths from bam files
	#bamfiles is a list of bam files
	
	#build data object from the bed info and depths from the bam files
	
	#Sam tools command needed
	# samtools depth -r chr1:43814968-43815086 NGS_HotSpot_Control_06-14-2017.275.bam | cut -f3 | paste -sd+ - | bc

	pass


def convert_to_JSON(analysis_home_dirs, data_obj):
	#function converts complex data obj (train or test) into a json object for printing to a file

	for home_dir in analysis_home_dirs:
		file_handle = open(str(home_dir)+'/JSON_SILO_trainingPYTHON.json', 'w')
		file_handle.write(json.dumps(data_obj, indent=4))

	print('training data written to file')

def read_training_JSON():
	
	filename = '/data/software/flypeC/silo/SILO_trainingPYTHON'
	file_handle = open(filename, 'r')
	
	data_obj = json.load(file_handle)
	print('file read from json')
	return data_obj

def add_training_sample(training_obj):
	#function takes in a data obj and adds a single samples worth of data to it.
	pass


def jackKnife(bcmatrix_file, output_dir):

	file_handle = open(bcmatrix_file, 'r')
	first_line = file_handle.readline()
	headers = first_line.strip('\n').split('\t')
	barcodes = first_line.strip('\n').split('\t')[2:] #slice the first two out
	
	try:
		os.mkdir(output_dir+'silo_jackknife_files')
		os.mkdir(output_dir+'silo_output_files')
	except Exception as e:
		raise ("Could not create output dir {}. Exception: {}".format(str(output_dir)+'/silo_jackknife_files', e) )

	silo_output_dir = str(output_dir)+'silo_output_files/'
	output_dir = str(output_dir)+'silo_jackknife_files/'
	for index, barcode in enumerate(headers):
		if index == 0 or index == 1:
			continue
		
		#prep file handles
		test_file_name = str(output_dir)+'/'+str(barcode+'_test.bcmatrix.xls')
		train_file_name = str(output_dir)+'/'+str(barcode+'_train.bcmatrix.xls')
		test_handle = open(test_file_name, 'w')
		train_handle = open(train_file_name, 'w')
		file_handle = open(bcmatrix_file, 'r')
		first_line = file_handle.readline()

		cut_headers = headers[:index]+headers[index+1:]
		test_handle.write('Gene\tTarget\t'+str(headers[index]) + '\n')
		train_handle.write('\t'.join(cut_headers) + '\n')

		for line in file_handle.readlines():
			
			#set flag for 0 values
			flag = False
			split_line = line.strip('\n').split('\t')
			#jackknifing will fail if any coverage value is 0
				# this is because at somepoint that 0 value will be the tested value and will attempt a divide by 0
			for value in split_line[2:]:
				if int(value) == 0:
					flag = True
					continue
			if flag == True:
				continue
			else:
				test_line = split_line[:2] + [split_line[index]]
				test_handle.write('\t'.join(test_line) + '\n')

				train_line = split_line[:index] + split_line[index+1:]
				train_handle.write('\t'.join(train_line) + '\n')

		test_handle.close()
		train_handle.close()
		file_handle.close()
		#call main and pass newly generated bcmatrix files to function
		output_file = silo_main(test_file_name, train_file_name, [silo_output_dir+barcode])

	return 'jack knife completed'


############
# filter functions if a filter file is provided to remove areas of consistant poor coverage
def parse_filter_file(filter_file):
	""" read in and parse filter file into a dictionary """

	file_handle = open(filter_file, 'r')

	filter_dict = {}
	for line in file_handle:
		if line[0] == '#':	continue
		line = line.strip().split('\t')
		filter_dict[line[1]] = line[0]
	
	file_handle.close()
	
	return filter_dict

def create_filtered_bcmatrix(bcmatrix_file, filter_file, analysis_home_dirs):
	""" read in bcmatrix file and check if the target exists in the filter_dict. if it does we want to keep it
		write that line out in the new bcmatrix filterd file.
	"""
	filter_dict = parse_filter_file(filter_file)

	file_handle = open(bcmatrix_file, 'r')
	output_file_fullpath = str(bcmatrix_file)+'.filtered.bcmatrix'
	output_handle = open(output_file_fullpath, 'w')
	headers = file_handle.readline()
	output_handle.write(headers)

	for line in file_handle:
		line = line.strip().split('\t')
		if line[1] in filter_dict:
			if line[0] == filter_dict[line[1]]:
				line_string = '\t'.join(line)
				line_string += '\n'
				output_handle.write(line_string)
			
	file_handle.close()
	output_handle.close()

	# copy the filtered bcmatrix file to every work dir for the run
	for work_dir in analysis_home_dirs:
		work_accession = os.path.basename(os.path.normpath(str(work_dir)))
		try:

			shutil.copyfile(output_file_fullpath, str(work_dir)+'/'+str(work_accession)+'.filtered.bcmatrix')
		except shutil.SameFileError:			pass

	return output_file_fullpath

#MAIN FUNCTION 
###################################################
###################################################

def silo_main(bcmatrix_test, bcmatrix_train,  analysis_home_dirs):

	#read in and build training obj. 
	#train_obj, train_raw_reads_obj = read_matrixFile('/data/software/flypeC/silo/set.TRAINING')
	train_obj, train_raw_reads_obj = read_matrixFile(bcmatrix_train)
	train_obj = calc_avg_sample_depth(train_raw_reads_obj, train_obj)
	train_obj = calc_porp_depth(train_obj)
	train_obj = calc_MPD_SD(train_obj)

	#convert_to_JSON(analysis_home_dirs, train_obj)
	#train_obj = read_training_JSON()

	# read in and build test obj
	print(bcmatrix_test)
	print('build test obj')
	test_obj, test_raw_reads_obj = read_matrixFile(bcmatrix_test)
	print('calc avg sample depth')
	test_obj = calc_avg_sample_depth(test_raw_reads_obj, test_obj)
	print('calc porp depth')
	test_obj = calc_porp_depth(test_obj)

	### determine gene amp in test obj
	print('find CN changes')
	test_obj = find_CN_changes(train_obj, test_obj)
	print('find gene CN amp')
	samples_list, CN_results = gene_amplification(test_obj)
	print('print results to file')
	output_file = pretty_print(samples_list, CN_results, bcmatrix_test, analysis_home_dirs)

	return output_file

