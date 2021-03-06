# SILO
#Silo README File
#Author: Nick Miller
#email: namiller2015@gmail.com
#05/07/2021

Silo is an algorithm for determining fold change in gene expression. Please see manuscript for more details.
Requirements:
Python 36
Bedtools
Execute silo in a Python 36 environment.

Commands available
	create_depth, create_static, run_silo,

create_depth:
	Create a depth summary file from a directory of BAM and BAI files using bedtools multicov 
	Output is a depth summary file that can be used as input for the silo algorithm
	
	Requirements:
		-bedfile: bed file correctly formatted
		-bams: directory of bam and bai files
		-output_directory: output directory name (will be created for you)
		-bedtools: Bedtools program should be installed and in your path. 
		
		Bed file Format
		Chrom	start		stop		target name	score	strand		gene name
		chr1	43814968	43815086	CHP2_MPL_1	0	+	.	GENE_ID=MPL


Create_static:
	Merge multiple depth summary files into one static training set
	
	Requirements:
		-depth_files: directory of depth files to use to create a static training set
		-output_directory: output directory name (will be created for you)

Run_silo:
	Run the Silo Algorithm in either static or jackknife mode
	Requirements:
		-static: enter True to run in Static mode using a training set
			Enter False to run in Jackknife mode using only 1 depth summary file
		
		-depth_testfile: filename to run SILO against
		-training_set: static training set required if running in static mode
		-output_directory: output directory name (will be created for you)


