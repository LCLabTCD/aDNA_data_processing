"""
HapSNP.py

Lara Cassidy

v1 05.08.2018 - HapSNP.py created from earlier pileup_filter.py by Eppie Jones

v2 28.03.2019 - removed triallelic filtering. Instead use allelic info in the ref bim to filter incorrect/unexpected bases.

v5 25.01.2021 - change GATK pileup output to verbose, including read length and mapping quality info, which can then be filtered

v9 21.02.2022 - multiple types of SNPCAP data (1240k, 390k etc.) now accounted for.

USAGE

usage: HapSNP_v5.py [-h] -data DATASET -I LIST_OF_BAMS [-R REF] [-p PARALLEL]
                    [-o OUTPUT_MERGE_FILE] [-mq MAPPING_QUAL] [-bp BP_LENGTH]
                    [-mem JAVA_MEM] [-gatk PATH_GATK]
DESCRIPTION

This program takes an input list of indexed BAM files and calls pseudo-haploid genotypes in these for a dataset of choice, outputting PLINK format files.
The following steps are performed

1. GATK pileup (verbose output)
2. Base quality 30 pileup filtering
3. Mapping quality and read length pileup filtering [defaults: 25, 34]
4. Removal of base calls not present in the reference bim. This is why it's important to have all allelic info present in the bim!
5. Random selection of one base call to create the pseudo-haploid genotype.

OUTPUT

Four intermediate files are generated and placed in an output directory ending in ".pifilter"
1. sample.dataset.bq30.pileup - pileup file filtered by base quality 30
2. sample.dataset.bq30.mq[mq].[bp]bp.pileup - pileup file filtered by user-specified mapping qualit and read length
2. sample.dataset.bq30.homozyg.pileup - Single base calls randomly selected from the above pileup file
4. sample.dataset.bq30.homozyg.potential_damage.snp_list - A list of SNPs in which damage could account for the final genotype (C/T and G/A SNPs, where sample is called as AA or TT) 

Plink Files are also outputted into the main working directory

sample.dataset.bq30.mq[mq].[bp]bp.homozyg.bed
sample.dataset.bq30.mq[mq].[bp]bp.homozyg.bim
sample.dataset.bq30.mq[mq].[bp]bp.homozyg.fam

Example dataset input string

/path/to/dataset

where the folder /path/to/ contains 

/path/to/dataset.pileup.bed
/path/to/dataset.pileup.bim

The pileup.bed file should look like this and include all positions you want called.

CHROM POS-1 POS
4 4999 5000

The pileup.bim file should contain the same positions but with extra info. Make sure you have biallelic info for all sites. cM distances don't matter and can be set to 0.

CHROM	SNPID	cM_Pos	POS	REF	ALT
4	SNP1	0	5000	A	T

Make sure your pileup file contigs match the reference genome you are using!


"""
#Import Statements

from __future__ import division
import os.path
import sys
import os
import subprocess
from subprocess import call
from subprocess import check_output
import errno
import csv
import gzip
from joblib import Parallel, delayed
import sys
import os
import random
import argparse
import signal
import pysam

## error log for gunzipping in python
def default_sigpipe():
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

##check directory 

def ensure_directory_exists(dir) :
    try:
        os.mkdir(dir)
    except OSError as exception :
        if exception.errno != errno.EEXIST :
            raise

## get current directory. final merged plink is labelled using the cwd name
cwd = os.getcwd()

## use parser module
parser = argparse.ArgumentParser()

## add arguments to accept
## required arguments
parser.add_argument("-data","--dataset",required=True,help="Path to base label for reference dataset (*pileup.bed and *pileup.bim files)")
parser.add_argument("-I","--list_of_bams",required=True,help="List of BAMS for which to call pseudohaploid genotypes")

## optional arguments
parser.add_argument("-R","--ref",help="Reference genome [hs37d5]",default="/Reference_Genomes/hs37d5/hs37d5.fa")
parser.add_argument("-p","--parallel",help="Number of samples to run in parallel - includes bwa samse/sampe but not aln [10]",type=int,default=10)
parser.add_argument("-o","--output_merge_file",help="Final merged PLINK file name [directory name]",default=cwd.split('/')[-1])
parser.add_argument("-mq","--mapping_qual",help="Mapping quality [25]",type=int,default=25)
parser.add_argument("-minbp","--min_bp_length",help="Minimum length filter [34]",type=int,default=34)
parser.add_argument("-maxbp","--max_bp_length",help="Maximum length filter [off]",type=int,default=-1)
parser.add_argument("-mem","--java_mem",help="Java memory for GATK [25]",type=int,default=25)
parser.add_argument("-gatk","--path_gatk",help="Path to GATK java file",default="/Software/GenomeAnalysisTK.jar")

## parse comand line arguments
args = parser.parse_args()
path_to_dataset = args.dataset
sample_list =  open(str(args.list_of_bams), 'r')

path_to_genome = args.ref
parallel_jobs = args.parallel
merge_label = args.output_merge_file
MQ = str(args.mapping_qual)
MINBP = str(args.min_bp_length)
MAXBP = str(args.max_bp_length)
Java_Mem = str(args.java_mem)
path_to_gatk = args.path_gatk

## get sample and labelling information
sample_IDs = ['.'.join(line.split('.')[0:-1]) for line in sample_list]
dataset_name= path_to_dataset.split('/')[-1].split('.')[0]
genome_name = path_to_genome.split('/')[-1].split('.')[0].split('_')[0]
base_name = '.' + dataset_name + '.' + genome_name + '.'

if int(MAXBP) < 0:
	infix= base_name + 'bq30.mq' + MQ + '.' + MINBP + 'bp'
else:
	infix= base_name + 'bq30.mq' + MQ + '.' + MINBP +  '-' + MAXBP + 'bp'


print sample_IDs

#Definitions

def pileup(full_ID):
	sample= full_ID.split('/')[-1].split('.')[0]
	try:
		file = open(sample + '.' + dataset_name + '.' + genome_name +  '.pileup' , 'r')
		print "Pileup already performed. Moving on..."
	except IOError:
                call('java -Xmx' + Java_Mem + 'g -jar ' + path_to_gatk + ' -T Pileup -R ' + path_to_genome + ' -I ' + full_ID + '.bam -o ' + sample + base_name + 'pileup -verbose -L ' + path_to_dataset + '.pileup.bed >' + sample + base_name + 'pileup.log 2>'  + sample + base_name + 'pileup.err.log', shell=True) 

def bq_filter(full_ID):
	sample= full_ID.split('/')[-1].split('.')[0]
	try:
		file =open(sample + base_name + 'pifilter/' + sample + base_name + 'bq30.pileup' , 'r')
		print "Pileup base quality filtering (30) already performed. Moving on..."

	except IOError:
		try :
			file =open(sample + base_name + 'bq30.pileup' , 'r')
			print "Pileup base quality filtering (30) already performed. Moving on..."
		except IOError:
			#Input files and strings
			input_file = open(sample + base_name + 'pileup','r')

			# for input pileup file, ignore header
			quality_filtered_file = open(sample + base_name + 'bq30.pileup', 'a')
			for line in input_file:
				if "REDUCE" in line:
					pass

			# Filter pileup file by base quality 30
				else:
					each = line.strip()
					each = each.split()
					my_SNP = each[3]
					quality = each[4]
					info = each[6].strip().split(',')
					bps = [read.split('@')[2] for read in info]
					mqs = [read.split('@')[3] for read in info]
					quality = quality.strip()
	
					counter = 0
					quality_list = []
					mySNP_list = []
					bp_list = []
					mq_list = []

					# Cycle through each value in quality string. Quality offset is 33.
					for i in quality[0:len(quality)]:
						offset_quality = ord(i)-33
	
						# If the quality of the base is greater than or equal to 30, store this base and the corresponding quality in a list (well 2 lists - one for bases, one for qualities.)
						if offset_quality >= 30:
							quality_list.append(i)
							mySNP_list.append(my_SNP[counter])
							bp_list.append(bps[counter])
							mq_list.append(mqs[counter])
						counter = counter + 1
	
					# Output high quality base calls for this position to file
					if len(quality_list) > 0:
						header_string = '\t'.join(each[0:3])
						mySNP_string = ''.join(mySNP_list)
						quality_string = ''.join(quality_list)
						bp_string = ','.join(bp_list)
						mq_string = ','.join(mq_list)
						quality_filtered_file.write(header_string + '\t' + mySNP_string + '\t' + quality_string + '\t' + bp_string + '\t' + mq_string + '\n')

			quality_filtered_file.close()


def final_filter(full_ID):
	
	sample= full_ID.split('/')[-1].split('.')[0]
	try:
		file =open(sample + base_name + 'pifilter/' + sample + infix + '.pileup' , 'r')
		print "Pileup mapping quality and length filtering already performed. Moving on..."

	except IOError:
		# Open quality filtered file for reading
		try :
			file = open( sample + infix + '.pileup', 'r')
			print "Pileup mapping quality and length filtering already performed. Moving on..."
		except IOError:
			try :
				quality_filtered_file = open(sample + base_name + 'bq30.pileup', 'r')
			except IOError:
				quality_filtered_file = open(sample + base_name + 'pifilter/' + sample + base_name + 'bq30.pileup', 'r')
			# Open output mapping quality and length filtered file for writing
			output = open(sample + infix + '.pileup', 'a')
       	                # Filter pileup file by mapping quality and length
          	      	for line in quality_filtered_file:
				each = line.strip()
				each = each.split('\t')
				chromosome = each[0]
				pos = each[1]
				my_SNP = each[3]
				quality = each[4]
				bp = each[5].split(',')
				mq = each[6].split(',')
				counter = 0
				mySNP_list = []
				quality_list = []
				bp_list = []
				mq_list = []

				# Cycle through each value in mq and bp string
				for i in bp[0:len(bp)]:
					if int(MAXBP) < 0 :
						if int(i) >= int(MINBP) and int(mq[counter]) >= int(MQ) :
							bp_list.append(i)
							mq_list.append(mq[counter])
							mySNP_list.append(my_SNP[counter])
							quality_list.append(quality[counter])
					elif int(i) >= int(MINBP) and int(i) <= int(MAXBP) and int(mq[counter]) >= int(MQ) :
						bp_list.append(i)
                                                mq_list.append(mq[counter])
                                                mySNP_list.append(my_SNP[counter])
                                                quality_list.append(quality[counter])

					counter = counter + 1
		
                	        if len(quality_list) > 0:
                        		header_string = '\t'.join(each[0:3])
	                               	mySNP_string = ''.join(mySNP_list)
	                                quality_string = ''.join(quality_list)
	                                bp_string = ','.join(bp_list)
	                                mq_string = ','.join(mq_list)
	                                output.write(header_string + '\t' + mySNP_string + '\t' + quality_string + '\t' + bp_string + '\t' + mq_string + '\n')

	                output.close()

def pseudohap(full_ID):
	sample= full_ID.split('/')[-1].split('.')[0]
	try:
		file =open(sample + infix + '.homozyg.map' , 'r')
		print "Pseudohaploidisation already performed. Moving on..."

	except IOError:

		#IDs
		sample= full_ID.split('/')[-1].split('.')[0]
		name = sample.split('-')[0]

		#Get data type info for plink ID
		if os.path.isfile(full_ID + ".bam"):
			read_groups =  subprocess.check_output("samtools view -h " + full_ID + ".bam | head -10000 |  awk '{if ($1 ~ /@RG/) print $0}'  ",shell=True, preexec_fn=default_sigpipe)
			read_groups = read_groups.replace('\n','\t')
			descriptions = []
			ids = []
			for i in read_groups.split("\t") :
				if 'DS:' in i :
					descriptions.append(":".join(i.split(":")[1:]))
				if 'ID:' in i :
					ids.append(":".join(i.split(":")[1:]))
			seq_types = []
			UDG_status = []
			lib_types = []
			if len(descriptions) < 1 :
				data_type = "NK"
			else :
				for entry in descriptions :
					info = entry.split('_')
					if len(info) < 4 :
						data_type = "NK"
					else :
						seq_types.append(info[0])
						UDG_status.append(info[1])
						lib_types.append(info[2])
				if len(set(seq_types)) == 1 :
					seq_type = seq_types[0]
				elif all("SNPCAP" in type for type in seq_types) :
					seq_type = "SNPCAP"
				else :	
					seq_type = "MIX"
				if len(set(lib_types)) == 1 :
					lib_type = lib_types[0]
				else :
					lib_type = "MIX"
				if len(set(UDG_status)) == 1 :
                        	        UDG = UDG_status[0]

				#If more than one type of UDG data is present in the file and some of it is non-UDG calculate the amount of non-UDG reads
				elif "nUDG" not in UDG_status :
					UDG = "UDG"
	
				else :
					nUDG_idxs = [n for (n, e) in enumerate(UDG_status) if e == 'nUDG']
					read_counts = []
					for i in nUDG_idxs:
						id = ids[i]
						count = subprocess.check_output("samtools view  " + full_ID + ".bam | head -100000 | grep " + id + " | wc -l",shell=True, preexec_fn=default_sigpipe)
						read_counts.append(int(count))
					nUDG_pc = int(sum(read_counts)/1000)
					if nUDG_pc < 10 :
						UDG = "UDG"
					elif nUDG_pc > 90 :
						UDG = "nUDG"
					else :
						UDG = str(100-nUDG_pc)+"pcUDG"
	
				data_type = seq_type+'_'+UDG+'_'+lib_type
		else : 
			data_type= subprocess.check_output("head -1 " + sample + ".*fam | awk '{print $2}' " ,shell=True)

		#Dictionaries
		SNP_dictionary = {}
		allele_dictionary = {}
		#Define Variables
		SNP_list=[]
		damage_list = []

		# Specify output files and folders
		output_folder = sample + base_name + 'pifilter'
		homozyg_file =open(sample + infix + '.homozyg.pileup' , 'a')
		triallelic_file=open(sample + infix +'.triallelic_sites.txt' , 'a')
		map_file =open(sample + infix + '.homozyg.map' , 'a')
		ped_file =open(sample + infix + '.homozyg.ped' , 'a')

		#For each line in bim file, join chromosome and SNP position (separated by a _) and create dictionary (called SNP_dictionary) relating this string to the SNP ID.
		#For each line in bim file, join chromosome and SNP position (separated by a _) and create dictionary (called allele_dictionary) relating this string to the possible alleles at this locus.

		input_bim = open(path_to_dataset + '.pileup.bim', 'r')

		for lines in input_bim:
			lines = lines.strip()
			col = lines.split('\t')
			chrom = col[0]
			SNP_ID = col[1]
			SNP_pos = col[3]
			ref = col[4]
			alt = col[5]
			SNP_info = chrom+'_'+SNP_pos # don't print "chr"
			alleles = ref + alt
			SNP_dictionary.update({SNP_info:SNP_ID})
			allele_dictionary.update({SNP_info:alleles})


# Now take this quality filtered file and randomly select one allele per SNP, making the dataset pseudo-homozygous.
# For each line in filtered file...
		try :
			quality_filtered_file =open(sample + infix + '.pileup' , 'r')
		except IOError:
                        quality_filtered_file = open(sample + base_name + 'pifilter/' + sample + infix +  '.pileup', 'r')		
		for line in quality_filtered_file:
				each = line.strip()
				each = each.split('\t')
				chromosome = each[0]
				chromosome_number =chromosome.replace("chr","")
				pos = each[1]
				my_SNP = each[3]
				quality = each[4]
				bp = each[5].split(',')
				mq = each[6].split(',')
				counter = 0
				quality_list = []
				mySNP_list = []
				mq_list = []
				bp_list = []
				SNP = chromosome+'_'+pos

				for i in my_SNP:
					if i in allele_dictionary[SNP]:
                                                quality_list.append(quality[counter])
                                                mySNP_list.append(i)
                        			mq_list.append(mq[counter])
                        			bp_list.append(bp[counter])
					else :
						triallelic_file.write(line)
                                        counter = counter + 1

				if len(mySNP_list) >= 1 :
					random_number= random.randint(0,len(mySNP_list)-1)
					one_allele = mySNP_list[random_number]
					corresponding_quality = quality_list[random_number]
					corresponding_bp = bp_list[random_number]
					corresponding_mq = mq_list[random_number]
					homozyg_file.write(each[0]+'\t'+each[1]+'\t'+each[2]+'\t'+one_allele+'\t'+corresponding_quality+'\t'+corresponding_bp+'\t'+corresponding_mq+'\n') 
					# Append SNP twice to SNP_list (to make SNP homozygous)
					SNP_list.append(one_allele+' '+one_allele)
					# Add SNP to map file
					map_file.write(chromosome_number+'\t'+SNP_dictionary[SNP]+'\t'+'0'+'\t'+pos+'\n')


		# Make ped file
		SNP_string = ' '.join(SNP_list)
		ped_file.write("ANC"+' '+name+'_'+data_type+' 0 0 0 -9 '+SNP_string)
		map_file.close()
		ped_file.close()
		homozyg_file.close()
		triallelic_file.close()


# Make list of SNPs compatible with being damage (i.e. my call is T where in the modern reference dataset there is a C, my call is an A where in the modern reference dataset there is a G)

		homozyg_file =open(sample + infix + '.homozyg.pileup' , 'r')
		damage_file=open(sample + infix + '.homozyg.damage_sites.txt' , 'a')

		for line in homozyg_file:
			each = line.strip()
			each = each.split('\t')
			chromosome = each[0]
			pos = each[1]
			my_SNP = each[3]
			SNP = chromosome+'_'+pos

			if my_SNP == 'T' and 'C' in allele_dictionary[SNP]:
				damage_file.write(SNP_dictionary[SNP] + '\n')
			if my_SNP == 'A' and 'G' in allele_dictionary[SNP]:
				damage_file.write(SNP_dictionary[SNP] + '\n')

		damage_file.close()

# Move all files you've created to folder specified on command line
		ensure_directory_exists(output_folder)
		cmd = "mv %s*bq30*pileup %s*bq30*txt  %s/"%(sample + '.' + dataset_name + '.' + genome_name, sample + '.' + dataset_name + '.' + genome_name, output_folder)
		os.system(cmd)

def plink_process(full_ID):
	sample= full_ID.split('/')[-1].split('.')[0]
	try:
                file =open(sample + infix + '.homozyg.bed', 'r')
		print "Binary plink files made! Moving on..."
	except IOError:
		call('plink1.9 --file ' + sample + infix + '.homozyg --make-bed --out ' + sample + infix + '.homozyg', shell=True)

Parallel(n_jobs=parallel_jobs)(delayed(pileup)(full_ID) for full_ID in sample_IDs)
Parallel(n_jobs=parallel_jobs)(delayed(bq_filter)(full_ID) for full_ID in sample_IDs)
Parallel(n_jobs=parallel_jobs)(delayed(final_filter)(full_ID) for full_ID in sample_IDs)
Parallel(n_jobs=parallel_jobs)(delayed(pseudohap)(full_ID) for full_ID in sample_IDs)
Parallel(n_jobs=parallel_jobs)(delayed(plink_process)(full_ID) for full_ID in sample_IDs)

#Final Merge

try:
	file=open(merge_label+'list','r')
	print "Merge list already made. Assuming merge complete. Delete merge list to run plink merge."
except IOError:
	#Make merge list for plink merge
	merge_list=open(merge_label+infix+'.list','a')
	for full_ID in sample_IDs:
		sample= full_ID.split('/')[-1].split('.')[0]
		merge_list.write(sample + infix + '.homozyg.bed ' + sample + infix + '.homozyg.bim ' + sample + infix + '.homozyg.fam' + '\n')
	merge_list.close()
	call('plink1.9 --merge-list ' + merge_label+infix+'.list  --make-bed --out ' + merge_label + infix + '.homozyg', shell=True)

