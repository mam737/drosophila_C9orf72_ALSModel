# Manisha Munasinghe - Last Updated: 6/24/19
# Script used for:
#   - QC, Mapping, and Read Counting of RNASeq Libraries
# For More Details See: 
#   - qc_functions.py
#   - 


import os
import subprocess
import qc_functions as qcf
import numpy as np 
import pandas as pd 


##############################################################################################################################
#####################################  RELEVANT DIRECTORES ###################################################################
##############################################################################################################################

os.chdir('/local/workdir/arh223/brc_data/')
#These are the directories where we want to store things
directory_1 = './raw_files/'
directory_2 = './raw_fastqc/'
directory_3 = './trimmomatic_files/'
directory_4 = './trimmomatic_fastqc/'
directory_5 = './ref_seq/'
directory_6 = './ref_seq/STARindex/'
directory_7 = './aligned_bam/'
directory_8 = './htseq_count_bam/'
##############################################################################################################################
######################  RELEVANT DIRECTORES ##################################################################################
##############################################################################################################################


######################  QUALITY CONTROL ######################################################################################

#Unzip raw fasta files
qcf.unzip_files(directory_1)

#Run FASTQC on raw data
qcf.fastqc_run(directory_1,directory_2)

#Trim Reads
qcf.trim_reads(directory_1,directory_3)

#Run FASTQC on trimmed data
qcf.fastqc_run(directory_3,directory_4)

#This is to generate a text file saving the read count information for raw and after trimming
if not os.path.exists('./a_brc_readcount.txt'):
	rc_dataframe = 	qcf.obtain_read_counts(directory_2,directory_4)
	rc_dataframe.to_csv('./a_brc_readcount.txt', header=True, index=False, sep='\t', mode='a')

######################  QUALITY CONTROL ######################################################################################


######################  MAPPING ##############################################################################################

#Generate STAR Index
qcf.star_gen_index(directory_6)

#Map reads
qcf.align_reads(directory_3,directory_7)

qcf.bam_index(directory_7)
######################  MAPPING ##############################################################################################

######################  READ COUNT ###########################################################################################

qcf.get_counts(directory_7,directory_8)


