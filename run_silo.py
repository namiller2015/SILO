# Author Nick Miller
# 5/25/2019
# contact at nmiller2@northshore.org

import subprocess
import argparse
import sys
import os

from silo_jackKnife import silo_main, jackKnife
from create_bcmatrix import bcmatrix
from create_static_training_set import get_bcmatrix_files, combine_bcmatrix_files

def convert_input_boolean(user_input):
	if user_input.lower() in ('yes', 'true', 't', 'y', '1'):
		return True
	elif user_input.lower() in ('no', 'false', 'f', 'n', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')


def init_parse_args():
	

	parser = argparse.ArgumentParser(description ="""Decide if you want to run SILO via static training set or jack knife method\n 
					Input required: \n
					static:	True or False \n 
					True will require a static Training file \n
					False will run SILO via jack knife on 1 depth summary file \n
					depth summary file name to evaluate \n
					static training set file \n
					output directory""" , formatter_class=argparse.RawTextHelpFormatter
					)
	
	parser.add_argument("-c", "--command", help = "Command to execute\n run_silo, create_depth, create_static", required=True)
	
	#arguments for running SILO in either jackknife or static mode
	parser.add_argument("-s", "--static", help = "Enter True or False to run silo in Static mode or jack knife mode", required=False)
	parser.add_argument("-tf", "--depth_testfile", help="filename of depth summary file to run SILO against", required=False)
	parser.add_argument("-ts", "--training_set", help="filename of static training set", required=False)
	parser.add_argument("-o", "--output_directory", help="Directory for Output files", required=True)
	
	#arguments for creating a bcmatrix file from a directory of BAM files
	parser.add_argument("-bed", "--bedfile", help = "bed file input", required=False)
	parser.add_argument("-bam", "--bams", help = "Directory of BAM and BAI files", required=False)
	
	#arguments for combining multiple bcmatrix files into a static training set
	parser.add_argument("-df", "--depth_files", help="Directory of depth files to use to create a static training set", required=False)

	args = parser.parse_args()
	

	if args.command == 'run_silo':
		run_silo(args)
	elif args.command == 'create_depth':
		create_bcmatrix_file(args)
	elif args.command == 'create_static':
		create_static_training_set(args)
	else:
		raise argParse.ArgumentTypeError('Command not recognized'.format(args.command) )
	
def create_bcmatrix_file(args):

	bcmatrix_obj = bcmatrix(args.bedfile, args.bams, args.output_directory)
	bcmatrix_obj.make_output_dir()

	output_file = bcmatrix_obj.create_bcmatrix()
	print(output_file)

	ready_to_run = """depth summary file generated. You are now ready to run SILO. Please execute the run_silo.py script"""
	print(ready_to_run)


def create_static_training_set(args):
	
	bcmatrix_files_list, bcmatrix_files_filename = get_bcmatrix_files(args.depth_files)
	combine_bcmatrix_files(bcmatrix_files_list, bcmatrix_files_filename, args.output_directory)


def run_silo(args):


	user_input_bool = convert_input_boolean(args.static)
	output_dir = create_dirs(args.output_directory)

	base_path = os.path.abspath(os.path.dirname(__file__))
	bcmatrix_file = os.path.join(base_path, args.depth_testfile)

	if user_input_bool == True and args.depth_testfile is None:
		parser.error("You selected to run SILO in static mode. Static command: {}. but no static_training file was found. static_set command {}".format(args.depth_testfile, args.depth_testfile) )

	if user_input_bool == True:
		#Running Static SILO
		
		#check if static training set input file exists
		full_static_file_path = os.path.join(base_path, args.depth_testfile)
		file_present = os.path.isfile(full_static_file_path)
		if file_present == False:
			parser.error("Static training set file could not be found. \n static_set file name input {}. Full static file path found {} \n Do you have the correct file path to your static file?".format(args.depth_testfile, full_static_file_path) )

		print('Runnning Static silo with depth summary file {}, static file {}, and output {}'.format(args.depth_testfile, args.depth_testfile, output_dir) )
		output_dir_list = [output_dir]
		silo_main(bcmatrix_file, full_static_file_path, output_dir_list)
		print('SILO complete')
	elif user_input_bool == False:
		#Running jack knife SILO
		print('Runnning SILO jack knife  with depth summary {} and output {}'.format(bcmatrix_file, output_dir) )
		jackKnife(bcmatrix_file, output_dir)
		print('SILO Jack Knife complete')
	else:
		raise argParse.ArgumentTypeError('Boolean value expected for static command. You input {}'.format(args.depth_testfile) )

def create_dirs(output_dir):
	
	base_path = os.path.abspath(os.path.dirname(__file__))
	output_dir = os.path.join(base_path, output_dir)
	
	full_dir_output = output_dir+'/'
	
	if os.path.isdir(full_dir_output) == False:

		try:
			os.mkdir(full_dir_output)
		except exception as e:
			print("could not create output dir {}. exception: {}".format(str(output_dir), e) )
			raise
	
	return full_dir_output




init_parse_args()


