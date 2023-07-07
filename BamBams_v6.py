"""

BamBams.py : Final BAM Merging, Processing and Basic Statistics

07.06.20

Lara Cassidy

usage: BamBams_v6.py [-h] [-mq MAPPING_QUAL] [-bp READ_BP_LENGTH]
                     [-clip SOFTCLIP] [-R REF] [-mem JAVA_MEM] [-p PARALLEL]
                     [-op OP_DUP] [-i BAMLIST] [-o OUTDIR]


v6. updated steps for processing non-TCD files.

===========
DESCRIPTION
===========


This program is designed to merge BAM files from across different sequencing runs following the 
TCD aDNA labelling system. It can also be used on single sample BAMS that do not follow this system.

The input is a list of BAM files that have been read-grouped, merged to a within sequencing run sample level and had duplicates removed. If you have used seqPro.py you'll find these here:

ls $seq_run/alignments/$ref_genome/bam_files/*trimmed$bp$q.*duprm.bm

You can include BAMS from as many sequencing runs as you want in the list - just make sure to give the full paths!

NOTE: BamBams.py expects all files to have been processed in the same way with the same suffix $trimming_infix.sorted.grouped.duprm.bam. If your files do not end in this fashion, please use the -tcd N switch and be sure to add appropriate reads groups and sort reads prior to running the program.

If the program detects separate files containing the same amplified library it will remove duplicates again. If you are not using TCD labelling system be sure to use the -tcd N option, which forces duplicate removal.

The program then carries out indel realignment using GATK. It will also apply mapping quality, read length and softclipping filters, which can be switched on or off by the user.

Finally, it will calculate coverage using qualimap.

IMPORTANT: You must give an output directory to write the merges to!


==================
SETUP INSTRUCTIONS
==================

You must have GATK and qualimap software on your server, as well as the softclip.py program


"""

#Import Statements

from __future__ import division
import sys
import os
import subprocess
from subprocess import call
from subprocess import check_output
import errno
import csv
import gzip
from joblib import Parallel, delayed
import fileinput
import argparse
from datetime import date
import signal 
import pysam
from contextlib import contextmanager
from io import BytesIO as StringIO

###set up stdout_redirector function for filter step

@contextmanager

def stdout_redirector(stream):

    old_stdout = sys.stdout

    sys.stdout = stream

    try:

        yield

    finally:

        sys.stdout = old_stdout


## Error log for gunzipping in python
def default_sigpipe():
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)


#Take command line arguments
parser = argparse.ArgumentParser()

parser.add_argument("-mq","--mapping_qual",help="Mapping quality [20]",type=int,default=20)
parser.add_argument("-bp","--read_bp_length",help="Read length filter [25]",type=int,default=25)
parser.add_argument("-clip","--softclip",help="Reduce base qualities at the ends of reads [default 2bp] to a phred score of 2 ",type=int,default=2)
parser.add_argument("-R","--ref",help="Reference genome [hs37d5]",default="/Reference_Genomes/hs37d5/hs37d5.fa")
parser.add_argument("-mem","--java_mem",help="Java memory for merging [25]",type=int,default=25)
parser.add_argument("-p","--parallel",help="Number of samples to run in parallel - includes bwa samse/sampe but not aln [10]",type=int,default=10)
parser.add_argument("-op","--op_dup",help="optical duplicate distance (default: 100). Please change to 2500 manually if duplicates from a NovaSeq, HiSeqX or HiSeq4000 are being removed",type=int,default=100)
parser.add_argument("-tcd","--tcd_labels",help="Does data follow the TCD labelling system? [Y]",default="Y")
parser.add_argument("-gatk","--gatk", help="Path to gatk jarfile", default = "/Software/GenomeAnalysisTK.jar")
parser.add_argument("-picard","--picard", help="Path to picard jarfile", default = "/Software/picard.jar")


#REQUIRED NO DEFALUTS
parser.add_argument("-i","--bamlist",help="list of bams",type=str)
parser.add_argument("-o","--outdir",help="output directory",type=str)

#Parse Arguments
args = parser.parse_args()
op_dup = str(args.op_dup)
parallel_jobs = args.parallel
MQ = str(args.mapping_qual)
Java_Mem = str(args.java_mem)
bp = str(args.read_bp_length)
clip = str(args.softclip)
tcd = str(args.tcd_labels)
bamlist  =  open(args.bamlist, 'r')
outdir = str(args.outdir)

#Get Ref Genome Path, Name and Length
path_to_genome = args.ref
Genome_Name = path_to_genome.split('/')[-1].split('.')[0]
ref_check = "/".join(path_to_genome.split('/')[0:-1])
Genome_Length = int(subprocess.check_output("cat " + path_to_genome + "*fai | awk '{SUM+=$2}END{print SUM}'",shell=True))

#Set paths to qualimap, picard and GATK
path_to_qualimap = str('.'.join(check_output('locate qualimap.jar | head -1', shell=True).rstrip('\n').split('.')[0:-1]))
path_to_gatk = str(args.gatk)
path_to_picard = str(args.picard)


#Set definition to create required directories
def ensure_directory_exists(dir) :
    try:
        os.mkdir(dir)
    except OSError as exception :
        if exception.errno != errno.EEXIST :
            raise

#Create output directory if it doesn't exist
ensure_directory_exists(outdir)

#Get list of samples and bam files
all_files = [line.rstrip() for line in bamlist]
samples = sorted(set([line.split("/")[-1].split('-')[0].split('.')[0] for line in all_files]))
infix =  "." + all_files[0].split("/")[-1].split('.')[1]

#####################################################################################
######################## BAM  MERGES AND DUPLICATE REMOVAL ##########################
#####################################################################################

#Store list of final files for later steps
bam_merges = []

#The below definition is run iteratively through samples.

def merge_duprm(merge_label,file_string):

	picard_merge= "java -Xmx" + Java_Mem + 'g -jar  ' + path_to_picard + ' MergeSamFiles '
	picard_duprm="java -Xmx" + Java_Mem + 'g -jar  ' + path_to_picard + ' MarkDuplicates OPTICAL_DUPLICATE_PIXEL_DISTANCE=' + op_dup + ' REMOVE_DUPLICATES=TRUE  I='

	try:
		file = open(outdir + "/" + merge_label + infix + '.sorted.grouped.duprm.bam.bai', 'r')
		print merge_label + " already merged and deduped if needed. Moving on..."
        except IOError:
		if duprm == "N"	:
			print "No duplicate removal required for " + merge_label
                	merge_command = picard_merge + ' I=' + file_string + ' O=' + outdir + "/" + merge_label + infix + '.sorted.grouped.duprm.bam 2>' + outdir + "/" + merge_label + infix + '.sorted.grouped.duprm.merge.log'
                	call(merge_command, shell=True)
			call('samtools index ' + outdir + "/" + merge_label + infix + '.sorted.grouped.duprm.bam',shell=True)
		else :
			try: 
				file = open(outdir + "/" + merge_label + infix + '.sorted.grouped.bam', 'r')
				print merge_label + " already merged. Deduping..."
				call(picard_duprm + outdir + "/" + merge_label + infix + '.sorted.grouped.bam O=' + outdir + "/" + merge_label  + infix + '.sorted.grouped.duprm.bam M=' + outdir + "/" + merge_label + infix + '.sorted.grouped.duprm.stats 2>' + outdir + "/" + merge_label + infix + '.sorted.grouped.duprm.log',shell = True)
				call('samtools index ' + outdir + "/" + merge_label + infix + '.sorted.grouped.duprm.bam',shell=True)
			except IOError:
	                	merge_command = picard_merge + ' I=' + file_string + ' O=' + outdir + "/" + merge_label + infix + '.sorted.grouped.bam 2>' + outdir + "/" + merge_label + infix + '.sorted.grouped.merge.log'
        	        	call(merge_command, shell=True)
				call(picard_duprm + outdir + "/" + merge_label + infix + '.sorted.grouped.bam O=' + outdir + "/" + merge_label  + infix + '.sorted.grouped.duprm.bam M=' + outdir + "/" + merge_label + infix + '.sorted.grouped.duprm.stats 2>' + outdir + "/" + merge_label + infix + '.sorted.grouped.duprm.log',shell = True)
				call('samtools index ' + outdir + "/" + merge_label + infix + '.sorted.grouped.duprm.bam',shell=True)


for i in samples :
	print i
#Define new sample files from your bam_list and old sample files that might be present in the output directoy
        sample_paths = [file for file in all_files if ("/" + i + "-") in file or ("/" + i + ".") in file] 
#Set up list of files for merging with picard
        Merge_Connector = ' I='
        file_string = Merge_Connector.join(sample_paths)

#Set up merge labels for output file
	label_lengths = [len((file.split("/")[-1].split(".")[0]).split("-")) for file in sample_paths]
	shortest_label = min(label_lengths)

	if shortest_label < 2 :
		merge_label= i
	elif shortest_label == 2 :
		BP_Tubes = sorted(set(['-' + line.split("/")[-1].split('-')[1].split(".")[0]  for line in sample_paths]))
		if len(BP_Tubes) == 1 :
			merge_label= i + str(BP_Tubes[0])
		else :
			merge_label= i
	elif shortest_label == 3 :
		BP_Tubes = sorted(set(['-' + line.split("/")[-1].split('-')[1] for line in sample_paths]))
		Extractions = sorted(set(['-' + line.split("/")[-1].split('-')[2].split(".")[0]  for line in sample_paths]))
		if len(BP_Tubes) > 1 :
                        merge_label= i 
		elif len(Extractions) > 1:
			merge_label= i + str(BP_Tubes[0])
		else :
			merge_label= i + str(BP_Tubes[0]) + str(Extractions[0])
        elif shortest_label == 4 :
		BP_Tubes = sorted(set(['-' + line.split("/")[-1].split('-')[1] for line in sample_paths]))
		Extractions = sorted(set(['-' + line.split("/")[-1].split('-')[2] for line in sample_paths]))
		Libraries = sorted(set(['-' + line.split("/")[-1].split('-')[3].split(".")[0]  for line in sample_paths]))
		if len(BP_Tubes) > 1 :
                        merge_label= i 
		elif len(Extractions) > 1:
			merge_label= i + str(BP_Tubes[0])
		elif len(Libraries) > 1 :	
			merge_label= i  + str(BP_Tubes[0]) + str(Extractions[0])
		else : 
			merge_label= i + str(BP_Tubes[0]) + str(Extractions[0])  + str(Libraries[0])

        elif shortest_label == 5 :
                BP_Tubes = sorted(set(['-' + line.split("/")[-1].split('-')[1] for line in sample_paths]))
                Extractions = sorted(set(['-' + line.split("/")[-1].split('-')[2] for line in sample_paths]))
                Libraries = sorted(set(['-' + line.split("/")[-1].split('-')[3] for line in sample_paths]))
                seq_run = sorted(set(['-' + line.split("/")[-1].split('-')[4].split(".")[0] for line in sample_paths]))
                if len(BP_Tubes) > 1 :
                        merge_label= i 
                elif len(Extractions) > 1:
                        merge_label= i + str(BP_Tubes[0])
                elif len(Libraries) > 1 :       
                        merge_label= i + str(BP_Tubes[0]) + str(Extractions[0])
                elif len(seq_run) > 1 : 
                        merge_label= i + str(BP_Tubes[0]) + str(Extractions[0]) + str(Libraries[0])
		else : 
			merge_label= i + str(BP_Tubes[0]) + str(Extractions[0])  + str(Libraries[0]) + str(seq_run[0])
	elif shortest_label == 6 :
                BP_Tubes = sorted(set(['-' + line.split("/")[-1].split('-')[1] for line in sample_paths]))
                Extractions = sorted(set(['-' + line.split("/")[-1].split('-')[2] for line in sample_paths]))
                Libraries = sorted(set(['-' + line.split("/")[-1].split('-')[3] for line in sample_paths]))
		PCRs = sorted(set(['-' + '-'.join(line.split("/")[-1].split('-')[4:6]).split(".")[0] for line in sample_paths]))
                if len(BP_Tubes) > 1 :
                        merge_label= i 
                elif len(Extractions) > 1:
                        merge_label= i + str(BP_Tubes[0])
                elif len(Libraries) > 1 :       
                        merge_label= i  + str(BP_Tubes[0]) + str(Extractions[0])
                elif len(PCRs) > 1 : 
                        merge_label= i + str(BP_Tubes[0]) + str(Extractions[0])  + str(Libraries[0])
                else : 
                        merge_label= i + str(BP_Tubes[0]) + str(Extractions[0])   + str(Libraries[0]) + str(PCRs[0])
	elif shortest_label > 6 : 
        	BP_Tubes = sorted(set(['-' + line.split("/")[-1].split('-')[1] for line in sample_paths]))
                Extractions = sorted(set(['-' + line.split("/")[-1].split('-')[2] for line in sample_paths]))
                Libraries = sorted(set(['-' + line.split("/")[-1].split('-')[3] for line in sample_paths]))
                PCRs = sorted(set(['-' + '-'.join(line.split("/")[-1].split('-')[4:6]) for line in sample_paths]))
		seq_run = sorted(set(['-' + '-'.join(line.split("/")[-1].split('-')[6:7]).split(".")[0] for line in sample_paths]))
                if len(BP_Tubes) > 1 :
                        merge_label= i 
                elif len(Extractions) > 1:
                        merge_label= i + str(BP_Tubes[0])
                elif len(Libraries) > 1 :       
                        merge_label= i + str(BP_Tubes[0])  + str(Extractions[0])
                elif len(PCRs) > 1 : 
                        merge_label= i  + str(BP_Tubes[0])  + str(Extractions[0])   + str(Libraries[0])
		elif len(seq_run) > 1:
			merge_label= i + str(BP_Tubes[0]) + str(Extractions[0])  + str(Libraries[0]) + str(PCRs[0])
                else : 
                        merge_label= i + str(BP_Tubes[0]) + str(Extractions[0])  + str(Libraries[0]) + str(PCRs[0]) + str(seq_run[0])


### Figure out whether duplicates need to be removed
	full_libs = []
	for i in sample_paths :
		file_libs = []
		proc = subprocess.Popen('samtools view -h ' + i + " | head -999 | awk '{if ($1 ~ /@RG/) print $0}'",shell=True, stdout=subprocess.PIPE, preexec_fn=default_sigpipe)
		for line in iter(proc.stdout.readline, ''):
			rgs_split = line.split("\t")
			for i in rgs_split :
				if "LB:" in i :
					file_libs.append((i.replace("LB:","")))
		file_libs  = set(file_libs)
		full_libs.extend(file_libs)

	print full_libs
	if tcd is "N" :
		duprm = "Y"
		infix= ".trimmedx"
	elif len(full_libs) == len(set(full_libs)) :
		duprm = "N"
	else :
		duprm = "Y"

####RUN THE MERGE #####

	merge_duprm(merge_label,file_string)
	bam_merges.append(merge_label + infix + '.sorted.grouped.duprm')

################################################
###########  INDELS AND FILTERS    #############
################################################

def indel(full_ID):
	try : 
		file = open( full_ID + '.indel.bai', 'r')
		print "Indels already realigned. Moving on..."
	except IOError:
		try :
			file = open( full_ID + '.intervals', 'r')
			print "Interval list already created. Moving on..."
			call('java -Xmx' + Java_Mem + 'g -jar ' + path_to_gatk + ' -T IndelRealigner -R ' + path_to_genome + ' -I ' + full_ID + '.bam -targetIntervals ' + full_ID + '.intervals -o ' + full_ID + '.indel.bam', shell=True)
		except IOError:
			call('java -Xmx' + Java_Mem + 'g -jar ' + path_to_gatk + ' -T RealignerTargetCreator -R ' + path_to_genome + ' -I ' + full_ID + '.bam -o ' + full_ID + '.intervals', shell=True)
			call('java -Xmx' + Java_Mem + 'g -jar ' + path_to_gatk + ' -T IndelRealigner -R ' + path_to_genome + ' -I ' + full_ID + '.bam -targetIntervals ' + full_ID + '.intervals -o ' + full_ID + '.indel.bam', shell=True)
			




def filter(full_ID):
	if MQ > 0 :
		file1 = full_ID + '.indel.mq' + MQ 
	else :
		file1 = full_ID + '.indel'
	if clip > 0 :
		file2 = file1 + '.clip' + clip + 'bp'
	else : 
		file2 = file1
        if bp > 0 :
                file3 = file2 + '.' + bp + 'bp'   
        else : 
                file3 = file2
	try:
        	file = open(file3 + '.bam', 'r')
                print "Filtering for " + full_ID + " already complete!"
        except IOError:	
		samfile = pysam.AlignmentFile(full_ID + ".indel.bam", "rb")			
		print "Filtering " + full_ID + "..."
		output = pysam.AlignmentFile(file3 +  '.bam', "wb", template=samfile) 

		for read in samfile :
			if read.mapping_quality >= int(MQ) and read.query_length >= int(bp) :
				clipped = "#" * int(clip)
				internal_pos =  read.query_length - int(clip)
				internal_string = str(read.qual[int(clip):internal_pos])
				soft_clipped_quality_string =  clipped + internal_string + clipped
				read.qual = soft_clipped_quality_string
				output.write(read)
		output.close()
		call('samtools index ' + file3  + '.bam',shell=True)

#Processing

os.chdir(outdir)
cwd = os.getcwd()
Parallel(n_jobs=1)(delayed(indel)(full_ID) for full_ID in bam_merges)
Parallel(n_jobs=parallel_jobs)(delayed(filter)(full_ID) for full_ID in bam_merges)

#New sample list
final_files = []
for i in bam_merges :
	new = [file for file in os.listdir(cwd) if file.startswith(i) and file.endswith(".bam.bai")][0]
	final_files.append(new)

print final_files

"""
#Coverage
for i in sample_IDs:
        if i.endswith("rmdup") :
                file = i + ".indel.softclip.bam" 
        else :
                file = i + ".rmdup.indel.softclip.bam"
	call(path_to_qualimap + ' bamqc -bam '  + file + ' -outdir ' + file + '_qualimap', shell=True)
"""

print "Bam Complete finished! Remember to back up your final files and delete any unneeded intermediate files!"


