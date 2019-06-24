# Manisha Munasinghe - Last Updated: 6/24/19
# Script containing functions used for:
#   - QC, Mapping, and Read Counting of RNASeq Libraries
# For More Details See:
#   - qc.py 
#   - 

import os
import subprocess
import numpy as np 
import pandas as pd 

#############################################################FUNCTIONS########################################################

def unzip_files(directory_in):
	for filename in os.listdir(directory_in):
		if filename.endswith('.gz'):
			print('Unzipping ' + filename)
			os.system('gunzip ' + directory_in + filename)


def fastqc_run(directory_in, directory_out):
	# This function runs FastQC on a file (filename) in a given directory (directory_in)
	# and outputs the FastQC report into the requested directory (direct	
	if not os.path.exists(directory_out):
		os.makedirs(directory_out)
		for filename in os.listdir(directory_in):
			if filename.endswith('fastq'):
				print(filename)
				os.system('/programs/FastQC-0.11.5/fastqc -q -o ' + directory_out + ' ' + directory_in + filename)
	#os.system = executes the command within in it in the shell
	#We are telling it to run FastQC to assess the library quality of the samples
	#Good Manual explanation of FastQC: https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf

def obtain_read_counts(directory_raw,directory_trim):
	### This function reads through  FastQC reports to determine the number of reads in the file
	read_dict = {}

	for filename in os.listdir(directory_raw):
		#Go through every file in the directory storing raw reads
		if filename.split('_')[-1]=='fastqc.zip':
			#If the file ends in zip, we know it contains read count information
			
			tag = filename.split('_')[0]+'_'+filename.split('_')[-2]+'_'+filename.split('_')[-1]
			#what's the unique identifier for this file/sample, tag will look something like 85-4_R1_fastqc
			
			inp_1 ='unzip -c ' + directory_raw + filename + ' ' +filename.split('.')[0]+'/fastqc_data.txt | sed -n 9p'		
			out_1 = subprocess.check_output(inp_1,shell=True)
			#We submit the input command (inp_1) to the shell and store the output (out_1)
			#This is to obtain the number of raw reads/sample

			inp_2 ='unzip -c ' + directory_trim+ 'trimmed_'+ tag + ' ' +'trimmed_'+tag.split('.')[0]+'/fastqc_data.txt | sed -n 9p'
			out_2 = subprocess.check_output(inp_2,shell=True)
			#We submit the input command (inp_2) to the shell and store the output (out_2)
			#This is to obtain the number of trimmed reads/sample
			
			key = tag.split('.')[0]
			value_1 = int(out_1.split()[-1])
			value_2 = int(out_2.split()[-1])
			read_dict[key]=[value_1,value_2]
			#We store the information we just obtained in a dictionary
			#where each key:value is as follows
			#key = sample name
			#value = [raw read count, trimmed read count]
	

	#Once we work through all of our files/samples, we translate all of the data in a dictionary
	#into a dataframe, which will allow us to output it as a text file for processing elsewhere
	rc_dataframe = pd.DataFrame(index=range(0,28),columns = range(0,3))
	rc_dataframe.columns = ['Sample','Raw','Trimmed']
	row = 0
	for key in read_dict:
		sample_dict = read_dict[key]
		rc_dataframe.loc[[row],0:3] = key.split('_')[0], sample_dict[0],sample_dict[1]
		row +=1
	rc_dataframe = rc_dataframe.sort_values(by='Sample')
	return(rc_dataframe)

def trim_reads(directory_in,directory_out):
	if not os.path.exists(directory_out):
		os.makedirs(directory_out)
		for filename in os.listdir(directory_in):
			os.system('java -jar /programs/trimmomatic/trimmomatic-0.36.jar SE ' + directory_in + filename + " " + directory_out + "trimmed_" +filename.split('_')[0] + '_' + filename.split('_')[-1] + ' LEADING:3 TRAILING:3 HEADCROP:10 SLIDINGWINDOW:4:20 MINLEN:20')


#STAR Manual Page: http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STAR.sandbox/doc/STARmanual.pdf
def star_gen_index(directory_out):
	if not os.path.exists(directory_out):
		os.makedirs(directory_out)
		os.system('/programs/STAR-2.5.2b/bin/Linux_x86_64/STAR --runMode genomeGenerate --runThreadN 8 --genomeDir ' + directory_out + ' --genomeFastaFiles ./ref_seq/dmel-all-chromosome-r6.18.fasta --sjdbGTFfile ./ref_seq/updated_dmel_r6.18_gffread.gtf')

def align_reads(directory_in,directory_out):
	if not os.path.exists(directory_out):
		os.makedirs(directory_out)
		print('Running STAR on Trimmed Files')
		for filename in os.listdir(directory_in):
			os.system('/programs/STAR-2.5.2b/bin/Linux_x86_64/STAR --genomeDir ./ref_seq/STARindex/ --readFilesIn ' + directory_in + filename + ' --runThreadN 12 --outFileNamePrefix ' + directory_out + filename.split('.')[0] + ' --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts')

#HTSEQ-Count Manual Page: http://htseq.readthedocs.io/en/master/count.html
def get_counts(directory_in,directory_out):
	if not os.path.exists(directory_out):
		os.makedirs(directory_out)
		for filename in os.listdir(directory_in):
			if filename.endswith('.bam'):
				os.system('/programs/HTSeq-0.11.2/bin/htseq-count -f bam -i gene_id ' + directory_in + filename + ' ./ref_seq/updated_dmel_r6.18_gffread.gtf > ' + directory_out + filename.split("Aligned")[0] + '_output_basename.counts')

###################################### FUNCTIONS #############################################################################
