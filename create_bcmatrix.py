#### script to take BAM files and generate a bcmatrix file as input for running SILO

# author: Nick Miller
# 3/14/2019

# check if bams need bai files generated

import subprocess
import argparse
import sys
import os


class bcmatrix:

	def __init__(self, bed, bam_input_dir, output_dir):
		
		self.bed = bed
		self.bam_input_dir = bam_input_dir
		
		#create absolute path to output dir
		base_path = os.path.abspath(os.path.dirname(__file__))
		self.output_dir = os.path.join(base_path, output_dir)
		
	def get_list_of_bam_files(self):
		# use given input_directory from user to gather list of bam files
		bam_list = []
		barcode_names_list = []
		for filename in os.listdir(self.bam_input_dir):
			if '.bam' in filename and '.bai' not in filename:
				full_file_path = os.path.join(self.bam_input_dir, filename)
				bam_list.append(full_file_path)
				barcode_names_list.append( filename.split('.bam')[0] )
		if len(bam_list) == 0:
			raise Exception('No bam files were found in {}. Are you sure that is the correct directory?'.format(bam_input_dir) )

		return bam_list, barcode_names_list

	def run_bedtools(self):
		# use bedtools multi cov command to generate coverage values from the bam files
		#bedtools multicov -bams IonXpress_001_rawlib.bam IonXpress_002_rawlib.bam -bed plugin_out/variantCaller_out.1169/CHP2.20131001.designed.bed > output_multicov_rawlib.txt
		
		bed = self.clean_bed_file()

		bam_files_list, barcode_names_list = self.get_list_of_bam_files()
		bam_files = ' '.join(bam_files_list)
		
		output_file = os.path.join(self.output_dir, "bedtools_output.txt")
		bed_tools_command = "bedtools multicov -bams %s -bed %s > %s" % (bam_files, bed, output_file)
		
		print('Running bedtools multicov on bam files')
		fp_error = open ("%s/errors.txt" % self.output_dir, "a")
		fp = subprocess.Popen([bed_tools_command], \
				shell=True, stdout=subprocess.PIPE, stderr=fp_error, universal_newlines=True  )
		fp.communicate()
		fp_error.close()
		print('Finished running bedtools multicov on bam files')
		return output_file, barcode_names_list
	
	def clean_bed_file(self):
		# use bedtools command to make sure bedfile is properly formed
		#sort <BED> | uniq > BED.noDuplicates
		
		print('Checking Bed file {} to make sure there are no duplicates'.format(self.bed) )
		output_file = os.path.join(self.output_dir, str(self.bed)+'_cleaned.bed')
		fp_error = open ("%s/errors.txt" % self.output_dir, "w")
		fp = subprocess.Popen(["sort %s | uniq > %s" %(self.bed, output_file  )], shell=True, stdout=subprocess.PIPE, stderr=fp_error, universal_newlines=True )
		fp.communicate()
		fp_error.close()
		print('Finished checking bed file')
		return output_file


	def create_bcmatrix(self):
		
		bedtools_output, barcode_list = self.run_bedtools()
		file_handle = self.file_handle = open(bedtools_output, 'r')

		bcmatrix_file = os.path.join(self.output_dir, 'bcmatrix')
		bcmatrix_handle = open(bcmatrix_file, 'w')
		bcmatrix_handle.write("""Gene\tTarget\t"""+'\t'.join(barcode_list)+'\n')
		print('Creating Bcmatrix file from multicov output')
		for line in self.file_handle:

			line = line.strip().split('\t')
			
			gene = line[7].split('=')[1]
			target = line[3]
			covg_depths = line[8:]

			bcmatrix_handle.write(gene+'\t'+target+'\t'+'\t'.join(covg_depths)+'\n')


		self.file_handle.close()
		bcmatrix_handle.close()
		print('Finished Creating Bcmatrix file {}. Ready to run silo'.format(bcmatrix_file) )

		return bcmatrix_file

	def make_output_dir(self):
		# check if output dir exists. if not create

		if not os.path.exists(self.output_dir):
			os.makedirs(self.output_dir)



def parse_args():
	parser = argparse.ArgumentParser(description ="""This program produces a bcmatrix file from BAM files for use in the SILO algorithm\n
			Inputs requried are:
				1. bed file
				2. Directory of BAM and BAI files
				3. Name of output directory (does not need to exist)""")

	parser.add_argument("-b", "--bedfile", help = "bed file input", required=True)
	parser.add_argument("-bm", "--bams", help = "Directory of BAM and BAI files", required=True)
	parser.add_argument("-o", "--output", help = "Name of output directory. Will be created if it does not exist already ", required=True)
	
	args = parser.parse_args()
	
	bcmatrix_obj = bcmatrix(args.bedfile, args.bams, args.output)
	bcmatrix_obj.make_output_dir()

	output_file = bcmatrix_obj.create_bcmatrix()
	print(output_file)

	ready_to_run = """bcmatrix file generated. You are now ready to run SILO. Please execute the run_silo.py script"""
	print(ready_to_run)


#parse_args()



