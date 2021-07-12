# Author Nick Miller
# 5/25/2019
# contact at nmiller2@northshore.org
import argparse
import sys
import os
import subprocess
#iterate through the input_dir and combine all the bcmatrix files in this dir into one static training set for SILO
# output that file in the output_directory

def init_parse_args():
	
	desc_str = """This script can be used to create a static training set for SILO from a directory of bcmatrix files \n
			Inputs Required: 
			Directory with bcmatrixs files to be turned into a static training set \n
			Output directory name where static training file will be put
			"""

	parser = argparse.ArgumentParser(description =desc_str, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-b", "--bcmatrixs_files", help="Directory of bcmatrix files to use to create a static training set", required=True)
	parser.add_argument("-o", "--output_directory", help="Directory for Output files", required=True)

	args = parser.parse_args()
	
	bcmatrix_files_list, bcmatrix_files_filename = get_bcmatrix_files(args.bcmatrixs_files)
	combine_bcmatrix_files(bcmatrix_files_list, bcmatrix_files_filename, args.output_directory)

def get_bcmatrix_files(input_dir):
	
	bcmatrix_files = []
	bcmatrix_files_filename = []
	for filename in os.listdir(input_dir):
		if 'bcmatrix' in filename:
			full_file_path = os.path.join(input_dir, filename)
			bcmatrix_files.append(full_file_path)
			bcmatrix_files_filename.append(filename)
	if len(bcmatrix_files) == 0:
		raise Exception('No bcmatrix files were found in {}. Are you sure that is the correct directory?'.format(input_dir) )

	return bcmatrix_files, bcmatrix_files_filename


def combine_bcmatrix_files(bcmatrix_files_list, bcmatrix_files_filename, output_dir_name):
	#combine the bcmatrix files into one
	# slice the first 2 columns out of each file into a "temp" file. 
	# Then stitch all the cleaned temp files into 1
	
	base_path = os.path.abspath(os.path.dirname(__file__))
	#get first file
	first_file = bcmatrix_files_list.pop(0)
	first_file_name = bcmatrix_files_filename.pop(0)
	first_file_fullpath = os.path.join(base_path, first_file)

	temp_dir = create_dirs(output_dir_name, '_temp_files')
	output_dir = create_dirs(output_dir_name, '_static_file')

	cleaned_file_list = []
	#now slice the remaning files and cat everything together
	for index, bcmatrix_file in enumerate(bcmatrix_files_list):
		
		output_filename = temp_dir+bcmatrix_files_filename[index]+'_cleaned'
		slice_cmmd = """cut -f 3- {} > {}""".format(bcmatrix_file, output_filename)
		fp_error = open ("%s/errors.txt" % temp_dir, "w")
		#shell command cut -d " " -f 3- input_filename > output_filename
		fp = subprocess.Popen([slice_cmmd], shell=True, stdout=subprocess.PIPE, stderr=fp_error, universal_newlines=True  )
		fp.communicate()
		fp_error.close()

		cleaned_file_list.append(output_filename)

	#create paste command of cleaned files
	cat_string = "paste {}".format(first_file_fullpath)
	for filename in cleaned_file_list:
		cat_string+= " {}".format(filename)

	cat_string += " > {}".format(output_dir+'static_training_bcmatrix')
	print(cat_string)
	child2_error = open ("%s/errors.txt" % output_dir, "w")
	#child2_output_static_file = open ("{}".format(output_dir+'static_training_bcmatrix') , "w")
	child2 = subprocess.Popen([cat_string], shell=True, stdout=subprocess.PIPE, stderr=child2_error, universal_newlines=True  )
	child2.communicate()
	child2_error.close()
	#child2_output_static_file.close()
		

def create_dirs(output_dir, dir_name):
	
	base_path = os.path.abspath(os.path.dirname(__file__))
	output_dir = os.path.join(base_path, output_dir)
	
	full_dir_output = output_dir+dir_name+'/'
	
	if os.path.isdir(full_dir_output) == False:

		try:
			os.mkdir(full_dir_output)
		except Exception as e:
			print("Could not create output dir {}. Exception: {}".format(str(output_dir)+'/{}'.format(dir_name), e) )
			raise
	
	return full_dir_output

#init_parse_args()





