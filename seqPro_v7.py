#!/bin/bash

"""
seqPro_v1.py

Lara Cassidy, Iseult Jackson and Emily Breslin

usage: seqPro.py [-h] [-R REF] [-n THREADS] [-p PARALLEL] [-o STATS]
                    [-mq MAPPING_QUAL] [-trimbp TRIM_BP_LENGTH]
                    [-mem JAVA_MEM] [-trimq TRIM_BP_QUAL] [-data DATA_TYPE]
                    [-index INDEX_CHECK] [-fish INDEX_FISH]
                    [-end SEQUENCE_STRATEGY] [-lib LIBRARY] [-RGCN CENTER]
                    [-RGPM PLATFORM]


v1 - 04/05/20 (lara edits, merged novaSeq.py and miSeqPro.py scripts)
v2 -  20/05/20 (iseult + lara edits, new coverage stats, user manual input for sequencer ID if needed)
bug alert: if you set parallel job number get_stats doesn't ouput the final file's stats. Haven't figured out why
v3 - 08/06/20 (lara edit - fix EC calculation for PE. Added a way to check if reads IDs are from older version of Illumina)
09/06/20 (bug fix: getting stuck at the mapDamage step)
08/03/21 - bug fix: correct mq read count for files where sample ids begin with same string (e.g. CpM3 and CpM30)
16/03/21 - fix p5 indices in dictionary to reverse complement for index fishing
v4 - indexing fishing revamp
New p5 and p7 indices have been added
v7 - 9/12/21: Add option to keep unmapped reads for metagenomic analysis; changed default -end to PE
===========
DESCRIPTION
===========

This script takes raw FASTQ files (PE or SE), trims them for adapters and base quality, aligns them to a chosen reference genome, merges to sample level and applies PCR duplicate and mapping quality filters. 

A statistics file with endogenous contents and median read lengths will be generated. mapDamage will be run on all libraries for each extract without UDG treatement.

If your FASTQs have index information, calculate the percentage of mismatches and generate index fished files (zero mismatches) with the -index and -fish switches (option: Y)

==================
SETUP INSTRUCTIONS
==================

Run this script in a directory labelled "raw_data" containing raw FASTQ files, within a sequencing run specific directory.

/SEQUENCE_RUN_DATE/raw_data/

FASTQC files must be labelled as below, or the -lib switch turned off (option: N). Samples must have 7 or more '-' delimited columns and controls 6 or less. The lane and read columns are optional.

If -lib switch is off, file merging and mapDamage will not be performed.

Sample examples:
    - Paired end: ID-BP-EXT-LIB-P7_P5-PCR-RUN-LANE_RX.fastq.gz
    - Single end: ID-BP-EXT-LIB-P7-PCR-RUN.fastq.gz

Control examples: 
    - H2O/AIR/EXT: CTRID-EXT-LIB-P7-PCR-RUN.fastq.gz
    - LIB: CTRID-P7_P5-PCR1-RUN-LANE_RX.fastq.gz
    - PCR: CTRID-P7-RUN.fastq.gz 

This script assumes single end (SE) data from a MiSeq unless the -end switch is set to paired end (option: PE). 

In paired end mode the NovaSeq platform is assumed, AdapterRemoval replaces cutadapt and bwa sampe is used for read pairs. 


=====
STEPS
=====

Files run in parallel for all steps except bwa aln and Picard programs. Set the number of samples to run in parallel with the -p switch.

1.  Gzip .fastq and .fq files if needed
2.  Alter read identifiers if -fish switch on (-fish Y)
3.  Adapter trimming with cutadapt or Adapter removal. Illumina TruSeq DNA adapters assumed. Minimum adapter overlap is set to 1. Base quality trimming and length filter set by user (-trimbp, -trimq).
4.  BWA alignment to chosen reference genome (-R /Path/to/ref/genome.fastq.gz) with parameters -o 2 -n 0.01 -l 165000. Threads can be set by user (-n).
5.  BWA samse or sampe conversion to BAM file. Reads are sorted (samtools sort) by chromosomal positions and read groups added. See more on read groups below.
6.  BAMs are merged to sample level using Picard (MergeSamFiles). This is done following file name unless the -lib switch is turned off (-lib N), in which case no merge is performed. Java memory set with -mem.
7.  PCR and optical duplicates are removed from BAMS (merged or otherwise) using Picard (MarkDuplicates). Java memory set with -mem. 
8.  Mapping quality filter is applied with samtools if set (-mq).
9.  Raw indexed fished (no mismatches) BAMs are generated if -fish switch is turned on. These can be used for merges with other sequencing runs (e.g. HiSeq)
10. mapDamage is run on extracts that have not undergone UDG treatment (-lib switch must be on).
11. Statistics are generated an placed in output file (-o).

===========
READ GROUPS
===========

Switch -lib switch off unless your file name is formatted correctly.

Fixed read groups:
Illumina sequencing technology assumed. Platform unit and read group ID are taken from the FASTQ read identifiers (@ line) and the file ID.

RGPL=ILLUMINA
RGPU=FLOWCELL:LANE (taken from FASTQ file)
RGID=SEQUENCER:RUN:FLOWCELL:LANE:FILE_ID

Automatic entry read groups:
RGSM and RGLB are taken from the file name unless lib switch is turned off (-lib N). If turned off, RGLB is set to "NK" and the sample take from the first column ('-' or '_' delimited) of the file ID.

RGSM=ID
RGLB=EXT-LIB-P7-PCR

Manual entry read groups:
Sequencing center and platform can be set manually (-RGCN, -RGPM). RGPM defaults to MiSeq for single end (-end SE) and NovaSeq for paired end (-end PE). RGCN defaults to TrinSeq_Dublin_Ireland.

RGCN=TrinSeq_Dublin_Ireland
RGPM=MiSeq

RGDS (description) includes the data_type (-data), the UDG status (nUDG, pUDG or UDG) and library type (dsLIB or ssLIB) taken from file name, the type of sequencing (-end) and the number of cycles (measured using the length of the first raw FASTQ read).

Examples:
  RGDS=WGS_nUDG_DS_PE75
  RGDS=MTCAP_UDG_SS_SE50

"""

#Import statements

from __future__ import division
import sys
import os
import subprocess
from subprocess import call
from subprocess import check_output
import errno
import csv
from joblib import Parallel, delayed
import multiprocessing
from multiprocessing import Queue
import argparse
import signal
import operator
from operator import itemgetter

#current working directory- where raw reads are stored
cwd = os.getcwd()

## use parser module
parser = argparse.ArgumentParser()

## Error log for gunzipping in python
def default_sigpipe():
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

## add arguments to accept
parser.add_argument("-R","--ref",help="Reference genome [hs37d5]",default="/Reference_Genomes/hs37d5/hs37d5.fa.gz")
parser.add_argument("-n","--threads",help="Number of threads for bwa [10]",type=int,default=10)
parser.add_argument("-p","--parallel",help="Number of samples to run in parallel - includes bwa samse/sampe but not aln [10]",type=int,default=10)
parser.add_argument("-o","--stats",help="Output statistics file name [seqPro_stats]",default="seqPro_stats")
parser.add_argument("-mq","--mapping_qual",help="Mapping quality [20]",type=int,default=0)
parser.add_argument("-trimbp","--trim_bp_length",help="Length filter pre-alignment [25]",type=int,default=25)
parser.add_argument("-mem","--java_mem",help="Java memory for merging [25]",type=int,default=25)
parser.add_argument("-trimq","--trim_bp_qual",help="Base pair quality trimming threshold for cutadapt [25]",type=int,default=25)
parser.add_argument("-data","--data_type",help="Is your data shotgun sequence (default: WGS) or capture? Input should be single string with no spaces. Common values - WGS, SNPCAP, YCAP, MTCAP, WGCAP, EXCAP",default="WGS")
parser.add_argument("-index","--index_check",help="Outputs meyer barcode in stats file. Default is on [Y]",default="Y")
parser.add_argument("-fish","--index_fish",help="Outputs percentage raw and aligned index mismatches. Removes whitespace from raw FASTQ read identifiers. Creates extra BAM files with no mismatches. Default is off [N]",default="N")
parser.add_argument("-end","--sequence_strategy",help="Single end (SE) or pair end (PE) sequencing? [PE]",default="PE")
parser.add_argument("-md","--map_damage",help="Run mapDamage stats (no bam recalibration)? [N]",default="N")
parser.add_argument("-un","--unmap",help="Output unaligned reads into BAM file",default="Y")
##Read groups
parser.add_argument("-lib","--library",help="Take sample and library ID info from the file name for read groups? [Y]",default="Y")
parser.add_argument("-RGCN","--center",help="Sequencing center/institute? [TrinSeq_Dublin_Ireland]",default="TrinSeq_Dublin_Ireland")
parser.add_argument("-RGPM","--platform",help="Sequencing Platform. This overwrites the defaults set by the PE/SE switch (-end) - NovaSeq for PE and MiSeq for SE ",default="none")
parser.add_argument("-seqid","--seqid",help="Sequencer unique ID. This is automatically taken from FASTQ read ids unless specificed",default="none")
parser.add_argument("-picard","--picard", help="Path to picard jarfile.", default="/Software/picard.jar")
## parse comand line arguments
args = parser.parse_args()
user_seq_id = args.seqid
parallel_jobs = args.parallel
n_threads = str(args.threads)
output = str(args.stats)
MQ = str(args.mapping_qual)
Java_Mem = str(args.java_mem)
trimbp = str(args.trim_bp_length)
trimq = str(args.trim_bp_qual)
index_check = str(args.index_check)
fish = str(args.index_fish)
library = str(args.library)
md = str(args.map_damage)
un = str(args.unmap)
path_to_picard = str(args.picard)
#Manual read group arguments
center = str(args.center)
platform = str(args.platform)
data_type = str(args.data_type)
end = str(args.sequence_strategy)

# Set default platform information if needed and optical duplicate distance (need to add option of other platforms)
if platform == "none" and end == "PE" :
	platform="NovaSeq"
	op_dup = "2500"
elif platform == "none" and end == "SE":
	platform="MiSeq"
	op_dup = "100"
#Other patterned flowcell platforms that require 2500 pixel distance
elif platform == "HiSeqX" or platform == "HiSeq4000" :
	op_dup= "2500"
else:
	op_dup="100"


#Set Alignment Genome : user inputs the path to the genome. From this we extract the name of the genome.
#We will use this name to create a directory in the Alignments folder for our output aligned files.
Path_to_genome = args.ref
Genome_Name = Path_to_genome.split('/')[-1].split('.')[0]
ref_check = "/".join(Path_to_genome.split('/')[0:-1])
Genome_Length = float(subprocess.check_output("cat " + Path_to_genome + ".fai | awk '{SUM+=$2}END{print SUM}'",shell=True))

def gzip(file):
        if file.endswith("fastq") or file.endswith("fq") :
                print "gzipping " + file
                call('gzip '+ file,shell=True)
        else :
                print file + " already gzipped. Moving on!"

#Run on all potential FASTQ files in the directory
print "Gzipping any unzipped files..."
Parallel(n_jobs=parallel_jobs)(delayed(gzip)(file) for file in os.listdir(cwd)) 


#Data description (RGDS) is based on the data type (e.g. WGS), sequencing strategy (PE or SE) and number of cycles. We also include library protocol information (e.g. UDG, ssLIB) on a sample by sample basis below
#Check the number of bps sequenced per read.
cycles =  str(subprocess.check_output("gunzip -c *q.gz |  head -2 | tail -n1 | awk '{print length($1)}'",shell=True, preexec_fn=default_sigpipe).rstrip('\n'))
seq_type = end + cycles

#Dictionary of P5 indexes:

#p5_dict = {"indexing1": "CCTGCGA", "indexing2": "TGCAGAG", "indexing3": "ACCTAGG", "indexing4": "TTGATCC", "indexing5": "ATCTTGC", "indexing6": "TCTCCAT", "indexing7": "CATCGAG", "indexing8": "TTCGAGC", "indexing9": "AGTTGGT", "indexing10": "GTACCGG", "indexing11": "CGGAGTT", "indexing12": "ACTTCAA", "indexing13": "TGATAGT"}

p5_dict = {"indexing1": "TCGCAGG", "indexing2": "CTCTGCA", "indexing3": "CCTAGGT", "indexing4": "GGATCAA", "indexing5": "GCAAGAT", "indexing6": "ATGGAGA", "indexing7": "CTCGATG", "indexing8": "GCTCGAA", "indexing9": "ACCAACT", "indexing10": "CCGGTAC", "indexing11": "AACTCCG", "indexing12": "TTGAAGT", "indexing13": "ACTATCA","indexing14":"TTGGATC","indexing15":"CGACCTG","indexing16":"TAATGCG","indexing19":"GAATCTC","indexing20":"CATGCTC","indexing21":"ACGCAAC","indexing22":"GCATTGG","indexing23":"GATCTCG","indexing24":"CAATATG","indexing25":"TGACGTC","indexing26":"GATGCCA","indexing27":"CAATTAC","indexing28":"AGATAGG","indexing29":"CCGATTG","indexing30":"ATGCCGC","indexing31":"CAGTACT","indexing32":"AATAGTA","indexing33":"CATCCGG","indexing34":"TCATGGT","indexing35":"AGAACCG","indexing36":"TGGAATA","indexing37":"CAGGAGG","indexing38":"AATACCT","indexing39":"CGAATGC","indexing40":"TTCGCAA","indexing41":"AATTCAA","indexing42":"CGCGCAG","indexing43":"AAGGTCT","indexing44":"ACTGGAC","indexing45":"AGCAGGT","indexing46":"GTACCGG","indexing48":"AATGATG","indexing49":"AGTCAGA","indexing50":"AACTAGA","indexing51":"CTATGGC"}

#Dictionary of P7 indexes:
p7_dict = {"indexing1": "TCGCAGG" , "indexing2": "CTCTGCA" , "indexing3": "CCTAGGT" , "indexing4": "GGATCAA" , "indexing5": "GCAAGAT" , "indexing6": "ATGGAGA" , "indexing7": "CTCGATG" , "indexing8": "GCTCGAA" , "indexing9": "ACCAACT" , "indexing10": "CCGGTAC" , "indexing11": "AACTCCG" , "indexing12": "TTGAAGT" , "indexing13": "ACTATCA" , "indexing14": "TTGGATC" , "indexing15": "CGACCTG" , "indexing16": "TAATGCG" , "indexing17": "AGGTACC" , "indexing18": "TGCGTCC" , "indexing19": "GAATCTC" , "indexing20": "CATGCTC" , "indexing21": "ACGCAAC" , "indexing22": "GCATTGG" , "indexing23": "GATCTCG" , "indexing24": "CAATATG" , "indexing25": "TGACGTC" , "indexing26": "GATGCCA" , "indexing27": "CAATTAC" , "indexing28": "AGATAGG" , "indexing29": "CCGATTG" , "indexing30": "ATGCCGC" , "indexing31": "CAGTACT" , "indexing32": "AATAGTA" , "indexing33": "CATCCGG" , "indexing34": "TCATGGT" , "indexing35": "AGAACCG" , "indexing36": "TGGAATA" , "indexing37": "CAGGAGG" , "indexing38": "AATACCT" , "indexing39": "CGAATGC" , "indexing40": "TTCGCAA" , "indexing41": "AATTCAA" , "indexing42": "CGCGCAG" , "indexing43": "AAGGTCT" , "indexing44": "ACTGGAC" , "indexing45": "AGCAGGT" , "indexing46": "GTACCGG" , "indexing47": "GGTCAAG" , "indexing48": "AATGATG" , "indexing49": "AGTCAGA" , "indexing50": "AACTAGA" , "indexing51": "CTATGGC" , "indexing52": "CGACGGT" , "indexing53": "AACCAAG" , "indexing54": "CGGCGTA" , "indexing55": "GCAGTCC" , "indexing56": "CTCGCGC" , "indexing57": "CTGCGAC" , "indexing58": "ACGTATG" , "indexing59": "ATACTGA" , "indexing60": "TACTTAG" , "indexing61": "AAGCTAA" , "indexing62": "GACGGCG" , "indexing63": "AGAAGAC" , "indexing64": "GTCCGGC" , "indexing65": "TCAGCTT" , "indexing66": "AGAGCGC" , "indexing67": "GCCTACG" , "indexing68": "TAATCAT" , "indexing69": "AACCTGC" , "indexing70": "GACGATT" , "indexing71": "TAGGCCG" , "indexing72": "GGCATAG" , "indexing73": "TTCAACC" , "indexing74": "TTAACTC" , "indexing75": "TAGTCTA" , "indexing76": "TGCATGA" , "indexing77": "AATAAGC" , "indexing78": "AGCCTTG" , "indexing79": "CCAACCT" , "indexing80": "GCAGAAG" , "indexing81": "AGAATTA" , "indexing82": "CAGCATC" , "indexing83": "TTCTAGG" , "indexing84": "CCTCTAG" , "indexing85": "CCGGATA" , "indexing86": "GCCGCCT" , "indexing87": "AACGACC" , "indexing88": "CCAGCGG" , "indexing89": "TAGTTCC" , "indexing90": "TGGCAAT" , "indexing91": "CGTATAT" , "indexing92": "GCTAATC" , "indexing93": "GACTTCT" , "indexing94": "GTACTAT" , "indexing95": "CGAGATC" , "indexing96": "CGCAGCC" , "indexing97": "GAGAGGC" , "indexing98": "GCTTCAG" , "indexing99": "ATATCCA" , "indexing100": "GTTATAC" , "indexing101": "CCTTAAT" , "indexing102": "CGCCAAC" , "indexing103": "TACTCGC" , "indexing104": "AGCGCCA" , "indexing105": "TAAGTAA" , "indexing106": "TTGGTCA" , "indexing107": "GTTGCAT" , "indexing108": "ATCCTCT" , "indexing109": "GGCGGTC" , "indexing110": "ATATGAT" , "indexing111": "GGTACGC" , "indexing112": "AAGAACG" , "indexing113": "CCGTTGA" , "indexing114": "AGCAATC" , "indexing115": "GCTCCGT" , "indexing116": "CAACTCT" , "indexing117": "AGACTCC" , "indexing118": "CTATCTT" , "indexing119": "AAGCAGT" , "indexing120": "GTTACCG" , "indexing121": "CCTAACG" , "indexing122": "ATCATAA" , "indexing123": "TGATAAC" , "indexing124": "TCTCCTA" , "indexing125": "CAGAGCA" , "indexing126": "CGGCTGG" , "indexing127": "CGCTATT" , "indexing128": "GATCGTC" , "indexing129": "ACGGCAG" , "indexing130": "GACCGAT" , "indexing131": "ACTTGCG" , "indexing132": "GTAAGCC" , "indexing133": "GCCATGC" , "indexing134": "ATAACGT" , "indexing135": "CGGACGT" , "indexing136": "GCGAGTA" , "indexing137": "ACGCGGA" , "indexing138": "GTCTAAT" , "indexing139": "GAAGCGT" , "indexing140": "CGGTAAG" , "indexing141": "GGTAACT" , "indexing142": "AATATAG" , "indexing143": "CCTCGCC" , "indexing144": "TTAATAG" , "indexing145": "CCGAAGC" , "indexing146": "TCGTTAT" , "indexing147": "GGCTCTG" , "indexing148": "CCAAGTC" , "indexing149": "CTTGGAA" , "indexing150": "TGAAGCT" , "indexing151": "GGTTGAC" , "indexing152": "CCGCCAT" , "indexing153": "ACCAGAG" , "indexing154": "GCTGAGA" , "indexing155": "CCTTCGC" , "indexing156": "CTGGCCT" , "indexing157": "CGCAAGG" , "indexing158": "TGAGAGA" , "indexing159": "AAGATTC" , "indexing160": "ATCGGTT" , "indexing161": "ACGAGCC" , "indexing162": "TGGATAC" , "indexing163": "ATTCCAG" , "indexing164": "ACTCATT" , "indexing165": "ACCTCGT" , "indexing166": "AGCTTAT" , "indexing167": "GCGATCT" , "indexing168": "CTCCAGT" , "indexing169": "GAACTTA" , "indexing170": "CCAATAA" , "indexing171": "AGGCGAG" , "indexing172": "CTCAGAT" , "indexing173": "AAGACGA" , "indexing174": "ACCGCTC" , "indexing175": "AACTGAC" , "indexing176": "GCAACTG" , "indexing177": "GAGTAAC" , "indexing178": "GATTAGG" , "indexing179": "AGTCGCT" , "indexing180": "CATAGAC" , "indexing181": "ACTACGG" , "indexing182": "GGCTAGC" , "indexing183": "CCATGAG" , "indexing184": "GCTGGTT" , "indexing185": "AGAGTAG" , "indexing186": "TATCAAC" , "indexing187": "ATACGCG" , "indexing188": "GACGTAC" , "indexing189": "CTACTTC" , "indexing190": "TGGTCGG" , "indexing191": "AGACCAT" , "indexing192": "TCCGAAC" , "indexing193": "GAGGTTG" , "indexing194": "TATGAGT" , "indexing195": "CTTCGTT" , "indexing196": "CCGTCCG" , "indexing197": "AACGTTA" , "indexing198": "GCATATT" , "indexing199": "ACCTTCC" , "indexing200": "TTCCGAG"}

#Ensure directories exist and are in the right places
def ensure_directory_exists(dir) :
    try:
        os.mkdir(dir)
    except OSError as exception :
        if exception.errno != errno.EEXIST :
            raise

ensure_directory_exists('../trimmed')
ensure_directory_exists('../alignments')
ensure_directory_exists('../alignments/' + Genome_Name)
ensure_directory_exists('../alignments/' + Genome_Name + '/sai_files/')
ensure_directory_exists('../alignments/' + Genome_Name + '/bam_files/')
ensure_directory_exists('../alignments/' + Genome_Name + '/merged_bams/')
if un == "Y" :
	ensure_directory_exists('../alignments/' + Genome_Name + '/unaligned/')


#Get list of all fastq files in the current directory
fastq_files = ['.'.join(file.split('.')[0:-2]) for file in os.listdir(cwd) if file.endswith('q.gz')]

#Remove whitespace from read identifiers if index switch is on
def read_id_merge(file):

	try:
		test = subprocess.check_output('gunzip -c ' + file + '.fastq.gz | head -1 | grep "#"',shell=True, preexec_fn=default_sigpipe).rstrip('\n')
		print "Read identifiers for " + file + " already merged! Moving on..."
        except subprocess.CalledProcessError as grepexc :
		print "Merging read identifiers for " + file
		if end == "PE" :
	                call("gunzip -c " + file + ".fastq.gz | sed 's/ 1:N:0/#1+2:N:0/g' | sed 's/ 2:N:0/#1+2:N:0/g' | gzip >" + file + ".fastq.gz.temp",shell=True)
			call("mv " + file + ".fastq.gz.temp " + file + ".fastq.gz",shell=True)
		elif end == "SE":
	                call("gunzip -c " + file + ".fastq.gz | sed 's/ 1:N:0/#1:N:0/g'  | gzip >" + file + ".fastq.gz.temp",shell=True)
			call("mv " + file + ".fastq.gz.temp " + file + ".fastq.gz",shell=True)
		else :
        		sys.exit('Please enter valid sequencing strategy, values: PE, SE')

if index_check == "Y" :
	Parallel(n_jobs=parallel_jobs)(delayed(read_id_merge)(full_ID) for full_ID in fastq_files)
elif index_check == "N" :
	print "No index checks or fishing will be carried out. Turn the -fish switch on (Y) if you'd like this to be done, but make sure you delete trimmed FASTQs and bams without index information before rerunning"
else :
	sys.exit("Please enter valid option for index fish. values: Y, N")


#currently formatted for lane-seperated
sample_list = [file for file in fastq_files if len(file.split('-')) > 6] 
samples = sorted(set([line.split('-')[0] for line in fastq_files]))
# Print the samples (not fastq files) being processed

print samples

###########################################################
###############     ADAPTER TRIMMING          #############
###########################################################


#Set the trimming infix
if trimq == str('0'):
	infix= '.trimmed' + trimbp + 'bp' 
else :
	infix= '.trimmed' + trimbp + 'bp' + trimq + 'q'

#Single End Data


def cutadapt(full_ID):
       	try:
               	file = open('../trimmed/' + full_ID + infix + '.fastq.gz', 'r')
               	print "Trimmed file for " + full_ID + " already exists! Moving on..."
       	except IOError:
               	call('cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O1 -m' + trimbp + ' -q ' + trimq + ' ' + full_ID + '.fastq.gz -o ../trimmed/' + full_ID + infix + '.fastq.gz >../trimmed/' + full_ID + infix + '.cutadapt.log 2>../trimmed/' + full_ID + infix + '.cutadapt.error.log',shell=True)

#Uses minlength trimbp (default 25): get most data possible; then filter for length of 34 nts for any downstream analyses to try and minimise ref bias.
def adapter_removal(full_ID):
	try:
		file=open('../trimmed/'+full_ID+infix+'.collapsed.gz', 'r')
		print "Trimmed fastq files for "+full_ID+" already exist! Moving on..."
	except IOError :
		call('AdapterRemoval --threads 2 --collapse --minadapteroverlap 1 --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTACTATTA --minlength '+ trimbp + ' --minquality '+trimq+' --gzip --trimns --trimqualities --file1 ' + full_ID+'_R1.fastq.gz --file2 '+full_ID+'_R2.fastq.gz --basename ../trimmed/'+full_ID+infix+' 2>../trimmed/'+full_ID+infix+'.AdapterRemoval.log', shell=True)

#Trimming single ends

if end == "SE" :
	print  "Trimming single end adapters..."
	Parallel(n_jobs=parallel_jobs)(delayed(cutadapt)(full_ID) for full_ID in fastq_files)

elif end == "PE":
	print  "Trimming paired end adapters..."
	basenames= sorted(set([file.split('_R')[0] for file in fastq_files]))
	Parallel(n_jobs=parallel_jobs)(delayed(adapter_removal)(full_ID) for full_ID in basenames)

else :
	sys.exit("Please choose valid sequence read type (options: PE or SE)")

print "Finished adapter trimming!"

#####################################################
###############     BWA ALN             #############
#####################################################

#Change into Trimmed_Ends directory for files to align

os.chdir('../trimmed/')
cwd = os.getcwd()
#collect sai_files

#list of trimmed files to align
trimmed_files = sorted([file for file in os.listdir(cwd) if infix+'.' in file and (file.endswith('fastq.gz') or file.endswith('ed.gz'))])
print trimmed_files
def align_file(trimmed_file):
	if "fastq" in trimmed_file :
		output = ".".join(trimmed_file.split(".")[0:-2])
	else :
		output = ".".join(trimmed_file.split(".")[0:-1])
	try:
        	file=open('../alignments/'+Genome_Name+'/sai_files/'+output+'.sai', 'r')
        	print "Aligned file for " + trimmed_file + " to " + Genome_Name + " already exists! Moving on..."
	except IOError:
        	call('bwa aln -l 165000 -n 0.01 -o 2 -t '+ n_threads + ' ' + Path_to_genome + ' ' + trimmed_file +  '> ../alignments/'+ Genome_Name + '/sai_files/'+output+'.sai 2> ../alignments/'+Genome_Name+'/sai_files/'+output+'.bwa.log', shell=True)

#Run this:
print "Aligning Files..."
Parallel(n_jobs=1)(delayed(align_file)(trimmed_file) for trimmed_file in trimmed_files)
print "Alignment finished... Yipee!"


#############################################################
################    SAMSE/SAMPE: SAI->BAM       #############
#############################################################

def samse_sampe_sample(sam_file):

#Determine output names on the basis of file name (pair1 + pair2 = .paired.)

	if "pair1" in sam_file :
		output = ".".join(sam_file.replace("pair1.truncated","paired").split(".")[0:-1])
	elif "fastq" in sam_file:
		output = ".".join(sam_file.split(".")[0:-2])
	else:
		output = ".".join(sam_file.split(".")[0:-1])
	check = sam_file.split(".")[0]

#Get read groups. The user can manually override these inputs
#Flowcell and sequencer information is taken from the first read of the FASTQ file
	seq_id = subprocess.check_output("gunzip -c " + check + "*gz |  head -1 | awk '{print $1}' | cut -f1-4 -d: ",shell=True, preexec_fn=default_sigpipe).rstrip('\n')

#Detect whether old FASTQ ID format
	check_seqid = int(seq_id.split(':')[3])
	if  check_seqid > 16 :
		pu = ':'.join(seq_id.split(':')[0:2])
		try : 
			seq_id = user_seq_id + ":" + pu
		except IOError :
			print "No sequencer ID detected in read ids. Please input manually using -seqid"
	else :
		pu = ':'.join(seq_id.split(':')[2:4])

#Sample and Library information is taken from the file name. Format must be: ID-BP-EXT-LIB-P5_P7-PCR-RUN-LANE_RX.FASTQ.gz
	full_ID= ''.join(sam_file.split('.')[0:1])
	if library == "Y" :
		sample = ''.join(full_ID.split('-')[0:1])
		if len(full_ID.split('-')) > 6 :
			lib = '-'.join(full_ID.split('-')[1:6])

#Control library information based on control type - Library, PCR or other (Air, Water, Extraction)
		else: 
			print full_ID + " has less than 7 info columns ( '-' delimited). Treating as control. Please use -lib, -sam and -des switches to enter read group information manually." 
			if full_ID.startswith('L'):
				lib =   '-'.join(full_ID.split('-')[1:3])
			elif full_ID.startswith('PCR') :
				lib =   '-'.join(full_ID.split('-')[1:2])
			else : 
				lib =   '-'.join(full_ID.split('-')[1:5])
		if len(lib.split('-')) > 4 :
			if "-ss" in lib :
				strand = "SS"
			else :
				strand = "DS"
			if "-pUDG" in lib :
				dmg = "pUDG"
			elif "-UDG" in lib :
				dmg = "UDG"
			else :
				dmg = "nUDG"
		else :
			strand = "CTR"
			dmg = "CTR"
		description = data_type + "_" + dmg + "_" + strand + "_" + seq_type
	elif library == "N" :
		sample = ''.join((''.join(full_ID.split('-')[0:1])).split('_')[0:1])
		lib = "NK"
		description = "NK"


	if un == "Y" :
		try:
        		file=open('../alignments/' + Genome_Name + '/unaligned/'+output+'.unaligned.bam', 'r')
        		print "Raw BAM file (sorted and grouped) and unaligned reads BAM file for " + sam_file + " aligned to "+ Genome_Name + " already exists! Moving on..."
		except IOError:
			try:
				# check if aligned bams have already been made
				file=open('../alignments/' + Genome_Name + '/bam_files/'+output+'.sorted.grouped.bam', 'r')
				print ("Raw aligned BAM for " + sam_file + "aligned to " + Genome_Name + " already exists! Creating unmapped BAM...")
				if "paired" in output:
					call('bwa sampe -r "@RG\\tID:' + seq_id + ':' + full_ID + '\\tLB:' + lib  + '\\tPL:ILLUMINA'  + '\\tSM:' + sample + '\\tPU:' + pu + '\\tCN:' +  center  + '\\tDS:' + description + '\\tPM:' + platform + '" ' + Path_to_genome + '  ../alignments/' + Genome_Name + '/sai_files/' + full_ID + infix + '.pair1.truncated.sai ../alignments/' + Genome_Name + '/sai_files/' + full_ID + infix + '.pair2.truncated.sai ' + full_ID + infix + '.pair1.truncated.gz   ' + full_ID + infix + '.pair2.truncated.gz 2> ../alignments/' + Genome_Name + '/sai_files/' + output  + '.unaligned.sampe.log | samtools view -Sb -f4 - > ../alignments/' +Genome_Name + '/unaligned/' + output + '.unaligned.bam 2>../alignments/' + Genome_Name + '/sai_files/'+ output + '.viewf4.log', shell=True)
				else:
					call('bwa samse -r "@RG\\tID:' + seq_id + ':' + full_ID + '\\tLB:' + lib  + '\\tPL:ILLUMINA'  + '\\tSM:' + sample + '\\tPU:' + pu + '\\tCN:' +  center  + '\\tDS:' + description + '\\tPM:' + platform + '" ' + Path_to_genome + '  ../alignments/' + Genome_Name + '/sai_files/' + output + '.sai ' + sam_file + ' 2> ../alignments/' + Genome_Name + '/sai_files/' + output  + '.unaligned.samse.log | samtools view -Sb -f4 - > ../alignments/' + Genome_Name + '/unaligned/' + output + '.unaligned.bam 2> ../alignments/' + Genome_Name + '/sai_files/' + output + '.viewf4.log', shell=True)
			except IOError:
				# make all raw bams 
			
				if "paired" in output :
			        	call('bwa sampe -r "@RG\\tID:' + seq_id + ':' + full_ID + '\\tLB:' + lib  + '\\tPL:ILLUMINA'  + '\\tSM:' + sample + '\\tPU:' + pu + '\\tCN:' +  center  + '\\tDS:' + description + '\\tPM:' + platform + '" ' + Path_to_genome + '  ../alignments/' + Genome_Name + '/sai_files/' + full_ID + infix + '.pair1.truncated.sai ../alignments/' + Genome_Name + '/sai_files/' + full_ID + infix + '.pair2.truncated.sai ' + full_ID + infix + '.pair1.truncated.gz   ' + full_ID + infix + '.pair2.truncated.gz 2> ../alignments/' + Genome_Name + '/sai_files/' + output  + ".sampe.log | bash -c 'tee >(samtools view -Sb -f4 - >../alignments/" + Genome_Name + '/unaligned/' + output + '.unaligned.bam 2>../alignments/' + Genome_Name + '/sai_files/'+ output + ".viewf4.log)' | samtools view -Sb -F4 -f2 -  2>../alignments/" + Genome_Name + '/sai_files/'+ output  + '.viewF4f2.log | samtools sort - -o ../alignments/' + Genome_Name + '/bam_files/' + output  + '.sorted.grouped.bam 2>../alignments/' + Genome_Name + '/sai_files/' + output +'.sort.log', shell=True)
				else:
       		 			call('bwa samse -r "@RG\\tID:' + seq_id + ':' + full_ID + '\\tLB:' + lib  + '\\tPL:ILLUMINA'  + '\\tSM:' + sample + '\\tPU:' + pu + '\\tCN:' +  center  + '\\tDS:' + description + '\\tPM:' + platform + '" ' + Path_to_genome + '  ../alignments/' + Genome_Name + '/sai_files/' + output + '.sai ' + sam_file + ' 2> ../alignments/' + Genome_Name + '/sai_files/' + output  + ".samse.log | bash -c 'tee >(samtools view -Sb -f4 - >../alignments/" + Genome_Name + '/unaligned/' + output + '.unaligned.bam 2>../alignments/' + Genome_Name + '/sai_files/'+ output + ".viewf4.log)' | samtools view -Sb -F4 -  2>../alignments/" + Genome_Name + '/sai_files/'+ output + '.viewF4.log | samtools sort - -o ../alignments/' + Genome_Name + '/bam_files/' + output + '.sorted.grouped.bam 2>../alignments/' + Genome_Name + '/sai_files/' + output + '.sort.log', shell=True)				
				call("samtools flagstat ../alignments/" + Genome_Name + "/bam_files/" + output + ".sorted.grouped.bam >../alignments/" + Genome_Name + "/bam_files/" + output + ".sorted.grouped.flagstat",shell=True)

	else:
		try:
        		file=open('../alignments/' + Genome_Name + '/bam_files/'+output+'.sorted.grouped.bam', 'r')
        		print "Raw BAM file (sorted and grouped) for " + sam_file + " aligned to "+ Genome_Name + " already exists! Moving on..."
		except IOError:
			if "paired" in output :
			        call('bwa sampe -r "@RG\\tID:' + seq_id + ':' + full_ID + '\\tLB:' + lib  + '\\tPL:ILLUMINA'  + '\\tSM:' + sample + '\\tPU:' + pu + '\\tCN:' +  center  + '\\tDS:' + description + '\\tPM:' + platform + '" ' + Path_to_genome + '  ../alignments/' + Genome_Name + '/sai_files/' + full_ID + infix + '.pair1.truncated.sai ../alignments/' + Genome_Name + '/sai_files/' + full_ID + infix + '.pair2.truncated.sai ' + full_ID + infix + '.pair1.truncated.gz   ' + full_ID + infix + '.pair2.truncated.gz 2> ../alignments/' + Genome_Name + '/sai_files/' + output  + '.sampe.log | samtools view -Sb -F4 -f2 -  2>../alignments/' + Genome_Name + '/sai_files/'+ output  + '.viewF4f2.log | samtools sort - -o ../alignments/' + Genome_Name + '/bam_files/' + output  + '.sorted.grouped.bam 2>../alignments/' + Genome_Name + '/sai_files/' + output +'.sort.log', shell=True)
			else:
       		 		call('bwa samse -r "@RG\\tID:' + seq_id + ':' + full_ID + '\\tLB:' + lib  + '\\tPL:ILLUMINA'  + '\\tSM:' + sample + '\\tPU:' + pu + '\\tCN:' +  center  + '\\tDS:' + description + '\\tPM:' + platform + '" ' + Path_to_genome + '  ../alignments/' + Genome_Name + '/sai_files/' + output + '.sai ' + sam_file + ' 2> ../alignments/' + Genome_Name + '/sai_files/' + output  + '.samse.log | samtools view -Sb -F4 -  2>../alignments/' + Genome_Name + '/sai_files/'+ output + '.viewF4.log | samtools sort - -o ../alignments/' + Genome_Name + '/bam_files/' + output + '.sorted.grouped.bam 2>../alignments/' + Genome_Name + '/sai_files/' + output + '.sort.log', shell=True)
			call("samtools flagstat ../alignments/" + Genome_Name + "/bam_files/" + output + ".sorted.grouped.bam >../alignments/" + Genome_Name + "/bam_files/" + output + ".sorted.grouped.flagstat",shell=True)





print "Creating Bams..."
sam_files=sorted([file for file in trimmed_files if 'pair2' not in file])
Parallel(n_jobs=parallel_jobs)(delayed(samse_sampe_sample)(full_ID) for full_ID in sam_files)
print "Samse/Sampe complete!"


#Change directory to Bam files and make list of bam files from sam_file list:
os.chdir('../alignments/' + Genome_Name + '/bam_files/')
cwd = os.getcwd()
original_bams = sorted([file.replace("pair1.truncated","paired").replace( "fastq.gz" , "sorted.grouped.bam").replace("gz","sorted.grouped.bam") for file in sam_files if "ingleton" not in file and "discard" not in file and "collapsed.truncat" not in file])


################################################
#############  INDEX FISHING    ################
################################################

###Fishes index matches from original raw BAM files and outputs them into new directory

def index_fish(file):
	file = '.'.join(file.split('.')[0:-1])
	full_ID= ''.join(file.split('.')[0:1])
	info = []

#Get index and name info for samples and controls
	if len(full_ID.split('-')) > 6 :
		index = '-'.join(full_ID.split('-')[4:5])
		match_id = str('-'.join(full_ID.split('-')[0:6]))
	else :
		if full_ID.startswith("L") :
			index = ''.join(full_ID.split('-')[1:2])
			match_id = str('-'.join(full_ID.split('-')[0:3]))
		elif full_ID.startswith("P") :
			index = ''.join(full_ID.split('-')[1:2])
			match_id = str('-'.join(full_ID.split('-')[0:2]))
		else :
			index = ''.join(full_ID.split('-')[2:3])
			match_id = str('-'.join(full_ID.split('-')[0:4]))

#Get sequence info for SE or PE reads
	if index.startswith("0") :
		index.replace("0","")
	if end == "SE" :
		check = "indexing" + str(index)
	        p7_sequence = p7_dict.get(check,"P7_index_not_found")
		sequence = p7_sequence
		seq_check = sequence

        elif end == "PE" :
                if "_" in index :
                        P7= index.split('_')[0]
                        P5= index.split('_')[1]
                        p7_check = "indexing" + str(P7)
                        p5_check = "indexing" + str(P5)
                        p7_sequence = p7_dict.get(p7_check,"P7_index_not_found")
                        p5_sequence = p5_dict.get(p5_check,"P5_index_not_found")
                        sequence = p7_sequence + "+" + p5_sequence
                        seq_check = sequence.replace("+","\+")
                else :
                        check = "indexing" + str(index)
                        sequence = p7_dict.get(check,"P7_index_not_found")
                        seq_check = sequence

#Create fished bams
	if fish == "Y" :
		try:
        		file_name=open(file + ".fish.bam", 'r')
        		print "index fishing already carried out for " + file + "! Moving on..."
		except IOError:
			call("samtools view -h " +  file + ".bam | awk '{if ($1 ~ /@/ || $1 ~ /:" + p7_sequence + "/) print $0}' | samtools view -b - >" + file + ".fish.bam",shell=True)
			call("samtools flagstat " + file + ".fish.bam >" + file + ".fish.flagstat",shell=True)
#get match counts
		try :
			file_name=open('../../../raw_data/' + match_id + '.p7match.txt', 'r')
		except IOError:
			if end == "SE" :
				call('gunzip -c ../../../raw_data/' + match_id + '*gz | grep "@" | grep "#"  | cut -f2 -d"#" | cut -f4 -d: | grep ' + p7_sequence + ' | wc -l  >../../../raw_data/' + match_id + '.p7match.txt' , shell=True)
			elif end == "PE" :
				call('gunzip -c ../../../raw_data/' + match_id + '*_R1*gz | grep "@" | grep "#"  | cut -f2 -d"#" | cut -f4 -d: | grep ' + p7_sequence + ' | wc -l >../../../raw_data/' + match_id + '.p7match.txt', shell=True)


	info.append(sequence)
	if match_id not in index_match_info :
		dict = {}
		dict[match_id] = info
		index_queue.put(dict)

index_queue = Queue()
index_match_info = {}

if index_check == "Y" :
	Parallel(n_jobs=parallel_jobs,prefer="threads")(delayed(index_fish)(file) for file in original_bams)
	while not index_queue.empty():
        	index_match_info.update(index_queue.get())

fish_bams = sorted([file.replace("pair1.truncated","paired").replace( "fastq.gz" , "sorted.grouped.fish.bam").replace("gz","sorted.grouped.fish.bam") for file in sam_files if "ingleton" not in file and "discard" not in file and "collapsed.truncat" not in file])



#####################################################################################
######################## BAM  MERGES AND DUPLICATE REMOVAL ##########################
#####################################################################################


#merge command for a list of files. The input is the merge label, which is dependent on the level of the merge (lanes, adapter removal output files, extractions etc.) and the list of files to be merged in picard input format
#following the merge duplicates are removed

def merge_duprm(merge_label,file_string):
	merge_label =  merge_label
	picard_merge= "java -Xmx" + Java_Mem + 'g -jar ' + path_to_picard +  ' MergeSamFiles '

	picard_duprm="java -Xmx" + Java_Mem + 'g -jar  ' + path_to_picard + ' MarkDuplicates OPTICAL_DUPLICATE_PIXEL_DISTANCE=' + op_dup + ' REMOVE_DUPLICATES=TRUE  I='

	if fish == "Y" :
		suffix=".sorted.grouped.fish"
	else:
		suffix=".sorted.grouped"

	try:
		file = open('../merged_bams/' + merge_label + infix + suffix +'.duprm.bam.bai', 'r')
		print merge_label + " already merged and deduped. Moving on..."
        except IOError:
		try:
			file = open('../merged_bams/' + merge_label + infix + suffix+ '.bam', 'r')
			print merge_label + " already merged. Deduping..."
			call(picard_duprm +  '../merged_bams/' + merge_label + infix + suffix + '.bam O=../merged_bams/' + merge_label  + infix + suffix + '.duprm.bam M=../merged_bams/' + merge_label + infix + suffix + '.duprm.stats 2>../merged_bams/' + merge_label + infix + suffix +'.duprm.log',shell = True)
			call('samtools index ../merged_bams/' + merge_label + infix + suffix + '.duprm.bam',shell=True)

		except IOError:
                	merge_command = picard_merge + ' I=' + file_string + ' O=../merged_bams/' + merge_label + infix + suffix + '.bam 2>../merged_bams/' + merge_label + infix + suffix+ '.merge.log'
                	call(merge_command, shell=True)
			call(picard_duprm + '../merged_bams/' + merge_label + infix + suffix + '.bam O=../merged_bams/' + merge_label  + infix + suffix + '.duprm.bam M=../merged_bams/' + merge_label + infix + suffix + '.duprm.stats 2>../merged_bams/' + merge_label + infix + suffix + '.duprm.log',shell = True)
			call('samtools index ../merged_bams/' + merge_label + infix + suffix + '.duprm.bam',shell=True)

###################
#SETTING UP MERGES#
###################

#Two types of merges.
#Sample level for downstream analysis
#Extraction level non-UDG merges for mapDamage

for_mapDamage = []
full_merges = []


#Get list of files in directory for each sample (with matching trimming parameters). We exclude singletons, discarded and truncated collapsed reads for PE data.
for i in samples :
	if fish == "Y" :
	        files = [file for file in fish_bams if file.startswith(i + "-") or file.startswith(i + ".")]
		suffix=".sorted.grouped.fish.duprm.bam"
	else :
		files = [file for file in original_bams if file.startswith(i + "-") or file.startswith(i + ".")]
		suffix=".sorted.grouped.duprm.bam"

#Get a list of BP tubes, Extactions, etc. for each sample
	if len(files[0].split('-')) > 6 :
	        BP_Tubes = sorted(set(['-' + line.split('-')[1] for line in files]))
		Extractions = sorted(set(['-' + line.split('-')[2] for line in files]))
		Libraries = sorted(set(['-' + line.split('-')[3] for line in files]))
		PCRs = sorted(set(['-' + '-'.join(line.split('-')[4:6]) for line in files]))
		Lanes = sorted(set(['-' + ('-'.join(line.split('-')[6:])).split('.')[0] for line in files]))
        	seq_run = str(sorted(set([('-'.join(line.split('-')[6:7])).split('.')[0] for line in files]))[0])
	elif len(files[0].split('-')) > 5 :
	        BP_Tubes = []
		Extractions = sorted(set(['-' + line.split('-')[1] for line in files]))
		Libraries = sorted(set(['-' + line.split('-')[2] for line in files]))
		PCRs = sorted(set(['-' + '-'.join(line.split('-')[3:5]) for line in files]))
		Lanes = sorted(set(['-' + ('-'.join(line.split('-')[5:])).split('.')[0] for line in files]))
        	seq_run = str(sorted(set([('-'.join(line.split('-')[5:6])).split('.')[0] for line in files]))[0])
	elif len(files[0].split('-')) > 3 :
		BP_Tubes = []
		Extractions = []
                Libraries = []
                PCRs = sorted(set(['-' + '-'.join(line.split('-')[1:3]) for line in files]))
                Lanes = sorted(set(['-' + ('-'.join(line.split('-')[3:])).split('.')[0] for line in files]))
                seq_run = str(sorted(set(['-'.join(line.split('-')[3:4]).split('.')[0] for line in files]))[0])
	else :
		BP_Tubes = []
		Extractions = []
		Libraries = []
		PCRs = sorted(set(['-' + '-'.join(line.split('-')[1:2]) for line in files]))
                Lanes = sorted(set(['-' + ('-'.join(line.split('-')[2:])).split('.')[0] for line in files]))
                seq_run = str(sorted(set(['-'.join(line.split('-')[2:3]).split('.')[0] for line in files]))[0])



#Set up list of files for merging with picard
	Merge_Connector = ' I='
	file_string = Merge_Connector.join(files)

#merge in a hierarchical manner, starting with BP merge.
	if len(BP_Tubes) > 1 :
		#Merge label is sample name
                merge_label= i + '-' + seq_run
		merge_duprm(merge_label,file_string)
		full_merges.append(merge_label + infix + suffix)

#If there is only one BP tube now move to extracts for full merges. Merge label will be sample and BP tube.
        elif len(Extractions) > 1 :

#for controls (only samples have BP tube)
		if len(BP_Tubes) > 0:
                	merge_label = i + str(BP_Tubes[0]) + '-' + seq_run
		else :
			merge_label = i + '-' + seq_run
         	merge_duprm(merge_label,file_string)
		full_merges.append(merge_label + infix + suffix)


#If there is only one BP tube and extract now move to libraries for full merges. Merge label will be sample, BP tube and extract.
        elif len(Libraries) > 1 :
		if len(BP_Tubes) > 0 :
	                merge_label = i + str(BP_Tubes[0]) + str(Extractions[0]) + '-' + seq_run
#for controls
		else :
			 merge_label = i +  str(Extractions[0]) + '-' + seq_run
                merge_duprm(merge_label,file_string)
		full_merges.append(merge_label + infix + suffix)

#Multiple PCRs for one library, extract and BP tube
        elif len(PCRs) > 1 :
		if len(BP_Tubes) > 0:
	                merge_label =  i + str(BP_Tubes[0]) + str(Extractions[0]) + str(Libraries[0]) + '-' + seq_run
#for controls
		elif len(Extractions) > 0 :
                        merge_label = i + str(Extractions[0]) + str(Libraries[0]) + str(PCRs[0]) + '-' + seq_run
                else :
                        merge_label = i  + str(PCRs[0]) + '-' + seq_run
	        merge_duprm(merge_label,file_string)
		full_merges.append(merge_label + infix + suffix)

        elif len(Lanes) > 1 :
                if len(BP_Tubes) > 0:
                        merge_label = i + str(BP_Tubes[0]) + str(Extractions[0]) + str(Libraries[0]) + str(PCRs[0]) + '-' + seq_run
                elif len(Extractions) > 0 :
                        merge_label = i + str(Extractions[0]) + str(Libraries[0]) + str(PCRs[0]) + '-' + seq_run
                else :
                        merge_label = i  + str(PCRs[0]) + '-' + seq_run
	        merge_duprm(merge_label,file_string)
		full_merges.append(merge_label + infix + suffix)

#file level of merging should occur for all paired end data, given multi-file output of adapter removal

	elif len(files) > 1 :
		merge_label = str(sorted(set([line.split('.')[0] for line in files]))[0])
		merge_duprm(merge_label,file_string)
		full_merges.append(merge_label + infix + suffix)
#if only one file exists for the sample (occurs only for single end data)

        else :
                merge_label = str(sorted(set(['.'.join(line.split('.')[0:-1]) for line in files]))[0])
		try:
                	file = open('../merged_bams/'+ merge_label +  '.duprm.bam.bai', 'r')
                	print merge_label + " does not require merging and already deduped. Moving on..."
        	except IOError:
                	print merge_label + " does not require merging. Deduping..."
			picard_duprm="java -Xmx" + Java_Mem + 'g -jar  ' + path_to_picard + ' MarkDuplicates OPTICAL_DUPLICATE_PIXEL_DISTANCE=' + op_dup +' REMOVE_DUPLICATES=TRUE  I='
			call(picard_duprm + merge_label + '.bam O=../merged_bams/' + merge_label + '.duprm.bam M=../merged_bams/' + merge_label + '.duprm.stats 2>../merged_bams/' + merge_label + '.duprm.log',shell = True)
			call('samtools index ../merged_bams/' + merge_label + '.duprm.bam',shell=True)
		full_merges.append(merge_label + '.duprm.bam')



################################################
###########    MAPPING QUALITY    ##############
################################################

def map_quality(bam_file) :
	if int(MQ) > 0 :
		try:
			bam_file = '.'.join(bam_file.split('.')[0:-1])
			file = open('../merged_bams/' + bam_file + '.mq' + MQ + '.bam.bai', 'r')
			print "Filtered (mq" + MQ + ")  BAM file aligned to " + Genome_Name + " for " + bam_file + " exists! Moving On..."
		except IOError:
			call('samtools view -b -q' + MQ + ' ../merged_bams/' + bam_file + '.bam >../merged_bams/' + bam_file + '.mq' + MQ + '.bam',shell = True) 
			call('samtools index ../merged_bams/' + bam_file + '.mq' + MQ + '.bam',shell = True) 
			call('samtools flagstat ../merged_bams/' + bam_file + '.mq' + MQ + '.bam >../merged_bams/' + bam_file + '.mq' + MQ + '.flagstat',shell = True)
	else :
		print "No MQ filter set (add using the -mq switch). Moving On..."


#run this function
print "Filtering for mapping quality..."
Parallel(n_jobs=parallel_jobs)(delayed(map_quality)(bam_file) for bam_file in full_merges)

####################################################
############        MAPDAMAGE         ##############
####################################################

#mapDamage. Creates another directory to store logs and output stats files.

def mapDamage(file) :
    name = '.'.join(file.split('.')[0:-1])
    ensure_directory_exists('../mapDamage')
    try:
        file=open('../mapDamage/results_' + name +'/' + name + '.mapDamage.log', 'r')
        print "mapDamage already carried out for " + name + "! Moving on..."
    except IOError:
        call('mapDamage -i ../merged_bams/' + file + ' -d ' + '../mapDamage/results_' + name +' -r ' + Path_to_genome + ' >../mapDamage/' + name + '.mapDamage.log 2>../mapDamage/'  + name + '.mapDamage_error.log' , shell=True)
        call('mv ../mapDamage/' + name  + '*mapDamage*log ../mapDamage/results_' + name +'/',shell=True)


for_mapDamage = [file for file in full_merges if 'LIB' in file]
if len(for_mapDamage) > 0 and md == 'Y' :
        print "Running mapDamage on non-UDG libraries for all files merged to library level or lower. Manually run mapDamage for extract level merges or higher. Note mapDamage has been run on files with no MQ filter"
        Parallel(n_jobs=parallel_jobs)(delayed(mapDamage)(file) for file in for_mapDamage)
        print "MapDamage complete!"

else:
        print "No files for mapDamage... moving on!"


################################################
#############       STATS       ################
################################################

#Definitions for getting read numbers and endogenous contents and adding each result to a queue (q) for parallel addition to our list of results
def get_stats(lib):

#Get column IDs for each file
	sample = ''.join(lib.split('-')[0])
	
#Samples
	if len(lib.split('-')) > 5 :
		BP_tube = '-'.join(lib.split('-')[1:2])
		extract = '-'.join(lib.split('-')[2:3])
		library = '-'.join(lib.split('-')[3:4])
		index = '-'.join(lib.split('-')[4:5])
		PCR = '-'.join(lib.split('-')[5:6])
		full_lib = BP_tube + "-" + extract + "-" + library + "-" + index + "-" + PCR #for grepping from MarkDup stats file

	else:
		BP_tube='-'
		if lib.startswith("L") :
			index = ''.join(lib.split('-')[1:2])
			extract = '-'
			library = '-'
			PCR = ''.join(lib.split('-')[2:3])
			full_lib = index + "-" + PCR #for grepping from MarkDup stats file

		elif lib.startswith("P") :
			index = ''.join(lib.split('-')[1:2])
			extract = '-'
			library = '-'
			PCR = '-'
			full_lib = index #for grepping from MarkDup stats file

		else :
			extract = ''.join(lib.split('-')[1:2])
			library = ''.join(lib.split('-')[2:3])
			index = ''.join(lib.split('-')[3:4])
			PCR = ''.join(lib.split('-')[4:5])
			full_lib = extract + "-" + library + "-" + index + "-" + PCR #for grepping from MarkDup stats file
	match_id = sample + "-" + full_lib

#Calculate read numbers and raw endogenous contents. We sum across lanes:

	if end == 'SE' :
		raw_reads = int(subprocess.check_output('cat ../../../trimmed/' + match_id + "*" + infix + '.cutadapt.log' + "| grep 'Total reads processed' |  awk '{print $4}' | sed 's/,//g' | awk '{SUM+=$1}END{print SUM}'",shell=True))
		trimmed_reads = int(subprocess.check_output('cat ../../../trimmed/' + match_id + "*" + infix + '.cutadapt.log' + "| grep 'Reads written' |  awk '{print $5}' | sed 's/,//g' | awk '{SUM+=$1}END{print SUM}' ",shell=True))
		average_length = round(float(subprocess.check_output('cat ../../../trimmed/' + match_id + "*" + infix + '.cutadapt.log' + "| grep 'Total written (filtered)' | awk '{print $4}' | sed 's/,//g' | awk '{SUM+=$1}END{print SUM}' ",shell=True))/trimmed_reads,2)

		#get average aligned read length from first file of that file
		if fish == 'Y' :
			aligned_files = sorted([file for file in os.listdir(cwd) if file.startswith(match_id) and infix+'.' in file and file.endswith('grouped.fish.bam')])
		else:
			aligned_files = sorted([file for file in os.listdir(cwd) if file.startswith(match_id) and infix+'.' in file and file.endswith('grouped.bam')])
		check_file = aligned_files[0]
		try :
			average_aligned_length = float(subprocess.check_output("samtools view " + check_file + " | awk '{print length($10)}' | awk '{ total += $1 } END { print total/NR }'",shell=True))
		except subprocess.CalledProcessError :
			average_aligned_length = 0
		aligned_reads = int(subprocess.check_output("head -1 " + match_id + "*" + infix + "*sorted.grouped.flagstat | grep total | awk '{print $1}' | awk '{SUM+=$1}END{print SUM}'",shell=True))
		if fish == 'Y' :
			fished_reads = int(subprocess.check_output("head -1 " + match_id + "*" + infix + "*sorted.grouped.fish.flagstat | grep total | awk '{print $1}' | awk '{SUM+=$1}END{print SUM}'",shell=True))

	else :
		raw_reads = int(subprocess.check_output('cat ../../../trimmed/' +  match_id + "*" + infix + ".*settings | grep 'Total number of read pairs' |  awk '{print $6}' | sed 's/,//g'  | awk '{SUM+=$1}END{print SUM}' ",shell=True))
		#ONLY COUNT PAIR1 PAIR2 and COLLAPSED AS PASSING TRIMMING FILTERS. WE ONLY COUNT ALIGNED READS FROM THESE FILES
	   	total_trim = int(subprocess.check_output('cat ../../../trimmed/' +  match_id + "*" + infix + ".*settings | grep 'Number of retained reads' |  awk '{print $5}' | sed 's/,//g'  | awk '{SUM+=$1}END{print SUM}' ",shell=True))
	   	collapsed_trim = int(subprocess.check_output('cat ../../../trimmed/' +  match_id + "*" + infix + ".*settings | grep 'Number of full-length collapsed pairs' |  awk '{print $6}' | sed 's/,//g'  | awk '{SUM+=$1}END{print SUM}' ",shell=True))
	   	collapsed_trun_trim = int(subprocess.check_output('cat ../../../trimmed/' +  match_id + "*" + infix + ".*settings | grep 'Number of truncated collapsed pairs' |  awk '{print $6}' | sed 's/,//g'  | awk '{SUM+=$1}END{print SUM}' ",shell=True))
	   	single1_trim = int(subprocess.check_output('cat ../../../trimmed/' +  match_id + "*" + infix + ".*settings | grep 'Number of singleton mate 1 reads' |  awk '{print $7}' | sed 's/,//g'  | awk '{SUM+=$1}END{print SUM}' ",shell=True))
	   	single2_trim = int(subprocess.check_output('cat ../../../trimmed/' +  match_id + "*" + infix + ".*settings | grep 'Number of singleton mate 2 reads' |  awk '{print $7}' | sed 's/,//g'  | awk '{SUM+=$1}END{print SUM}' ",shell=True))
		pair_trim = (total_trim - collapsed_trim - collapsed_trun_trim - single1_trim - single2_trim)/2
		trimmed_reads = collapsed_trim + pair_trim
 		average_length = round(float(subprocess.check_output('cat ../../../trimmed/' +  match_id +  "*" + infix + ".*settings | grep 'Average length of retained reads' |  awk '{print $6}' | sed 's/,//g' | awk '{ sum += $1 } END { if (NR > 0) print sum / NR }' ",shell=True)),2)

		if fish == 'Y' :
			p_aligned_files = sorted([file for file in os.listdir(cwd) if file.startswith(match_id) and infix+'.' in file and file.endswith('paired.sorted.grouped.fish.bam')])
			c_aligned_files = sorted([file for file in os.listdir(cwd) if file.startswith(match_id) and infix+'.' in file and file.endswith('collapsed.sorted.grouped.fish.bam')])
		else :
			p_aligned_files = sorted([file for file in os.listdir(cwd) if file.startswith(match_id) and infix+'.' in file and file.endswith('paired.sorted.grouped.bam')])
			c_aligned_files = sorted([file for file in os.listdir(cwd) if file.startswith(match_id) and infix+'.' in file and file.endswith('collapsed.sorted.grouped.bam')])

		#CHECK AVERAGE ALIGNED LENGTH USING ONE FILE ONLY (ignore multiple lanes)
		p_check_file = p_aligned_files[0]
		c_check_file = c_aligned_files[0]
		p_weight = float(subprocess.check_output("samtools view " + p_check_file + " |  wc -l ",shell=True))
		c_weight = float(subprocess.check_output("samtools view " + c_check_file + " | wc -l ",shell=True))
		try  :
			if p_weight > 0 : 
				p_average_aligned_length = float(subprocess.check_output("samtools view " + p_check_file + " | awk '{if ($9 > 0) print $9}' | awk '{ total += $1 } END { print total/NR }'",shell=True))
				c_average_aligned_length = float(subprocess.check_output("samtools view " + c_check_file + " | awk '{print length($10)}' | awk '{ total += $1 } END { print total/NR }'",shell=True))
				average_aligned_length = round(((p_average_aligned_length * p_weight) + (c_average_aligned_length * c_weight))/(c_weight + p_weight) ,2)
			else :
				average_aligned_length = float(subprocess.check_output("samtools view " + c_check_file + " | awk '{print length($10)}' | awk '{ total += $1 } END { print total/NR }'",shell=True))
		except subprocess.CalledProcessError :                        
                	average_aligned_length = 0

		aligned_collapsed = int(subprocess.check_output("head -1 " + match_id + "*" + infix + "*collapsed.sorted.grouped.flagstat | grep total | awk '{print $1}' | awk '{SUM+=$1}END{print SUM}'",shell=True))
		aligned_paired = int(subprocess.check_output("head -1 " + match_id + "*" + infix + "*paired.sorted.grouped.flagstat | grep total | awk '{print $1}' | awk '{SUM+=$1}END{print SUM}'",shell=True))/2
		aligned_reads = aligned_collapsed + aligned_paired
		if fish == 'Y' :
			fished_collapsed = int(subprocess.check_output("head -1 " + match_id + "*" + infix + "*collapsed.sorted.grouped.fish.flagstat | grep total | awk '{print $1}' | awk '{SUM+=$1}END{print SUM}'",shell=True))
			fished_paired = int(subprocess.check_output("head -1 " + match_id + "*" + infix + "*paired.sorted.grouped.fish.flagstat | grep total | awk '{print $1}' | awk '{SUM+=$1}END{print SUM}'",shell=True))/2
			fished_reads = fished_collapsed + fished_paired

# Get alignment and complexity information from MarkDup stats
	if fish == 'Y' :
		pcr_duplicates=int(subprocess.check_output('cat ../merged_bams/' + sample +  "-*" + infix + "*sorted.grouped.fish.duprm*stats  | awk '{if ($1 ~ /" + full_lib + "/) print $0}' | sort -u | awk '{print $6+$7}'",shell=True))
		mq_files = [file for file in os.listdir("../merged_bams/") if file.endswith(".fish.duprm.bam") and (file.startswith(sample+"-") or file.startswith(sample+".")) and infix in file]
		try :
			dup_rate = round(float(pcr_duplicates/fished_reads),6)
		except ZeroDivisionError  :
			dup_rate = 0
		dup_ec =  round(float( (fished_reads - pcr_duplicates)/trimmed_reads),6)
	else: 
		pcr_duplicates=int(subprocess.check_output('cat ../merged_bams/' + sample +  "-*" + infix + "*sorted.grouped.duprm*stats  | awk '{if ($1 ~ /" + full_lib + "/) print $0}' | sort -u | awk '{print $6+$7}'",shell=True))
		mq_files = [file for file in os.listdir("../merged_bams/") if file.endswith(".grouped.duprm.bam") and (file.startswith(sample+"-") or file.startswith(sample+".")) and infix in file]
		try :
			dup_rate = round(float(pcr_duplicates/aligned_reads),6)
		except ZeroDivisionError  :
                	dup_rate = 0
		dup_ec =  round(float((aligned_reads - pcr_duplicates)/trimmed_reads),6)


	if len(mq_files) == 1 :
		mq_file = str(mq_files[0])
	else:
		mq_file = str([file for file in mq_files if "LIB" not in file][0])	
	mq20_reads = float(subprocess.check_output("samtools view -F128 -q20 ../merged_bams/" + mq_file + " | grep " + match_id + " | wc -l" ,shell=True))

	#library_size=int(subprocess.check_output('cat ' + sample +  "*" + infix + "*stats  |  awk '{if ($1 ~ /" + full_lib + "/) print $0}' | sort -u | awk '{print $10}'",shell=True))

    #Create row of sample values and add to list of lists.
	if index_check is "Y" :
		sequence = index_match_info.get(match_id,"not_found")[0]
	else :
		sequence = "NA"
	if fish == "Y" :
		fished_raw = int(subprocess.check_output('cat ../../../raw_data/' + match_id + '.p7match.txt',shell=True))
		index_match_raw = round(float(fished_raw/raw_reads),4)
		index_match_aligned = round(float(fished_reads/aligned_reads),4)
		final_aligned = fished_reads
	else :
		index_match_raw = "NA"
		index_match_aligned = "NA"
		final_aligned = aligned_reads
	try:
		Sample_Row = [sample,BP_tube,extract,library,index,PCR,raw_reads,trimmed_reads,final_aligned,round(((final_aligned)/trimmed_reads),4),pcr_duplicates,dup_ec,dup_rate,mq20_reads,round(((mq20_reads)/trimmed_reads), 4), average_length,average_aligned_length,round((average_aligned_length*mq20_reads)/Genome_Length,6),sequence,index_match_raw,index_match_aligned,end,cycles]
		q.put(Sample_Row)
	except ZeroDivisionError:
		Sample_Row = [sample,BP_tube,extract,library,index,PCR,raw_reads,trimmed_reads,final_aligned,"Division_Error",pcr_duplicates,"Division_Error",dup_rate,mq20_reads,"Division_Error", average_length,average_aligned_length,round((average_aligned_length*mq20_reads)/Genome_Length,6),sequence,index_match_raw,index_match_aligned,end,cycles]
     		q.put(Sample_Row)


list_of_lists = [['Sample_ID','BP_Tube','Extraction','Library','Index','PCR','Raw_Reads','Trimmed_Reads_(' + trimbp + 'bp_' + trimq + 'q)','Aligned_Reads','Raw_EC','PCR_Duplicates','Duprm_EC','Duplication_Rate','Final_Reads_(mq20)','Final_EC','Avg_Read_Length_Raw','Avg_Read_Length_Aligned','Estimated_Coverage','Index_Sequence(s)','Index_Match_Raw','Index_Match_Aligned','Seq_Type','Cycles']]

#Create Queue  to be shared across different processes

q = Queue()

#Runs stats definitions. Make sure the file does not already exist.

print "Calculating Stats..."
# get list of unique libraries to run stats on
libs = sorted(set(['-'.join(line.split('-')[0:-1]) for line in original_bams]))

try:
    file=open('../../../' + output + '.csv', 'r')
    print "Stats already calculated! Moving on..."

except IOError:

    Parallel(n_jobs=parallel_jobs,prefer="threads")(delayed(get_stats)(lib) for lib in libs)
    while not q.empty():
        list_of_lists.append(q.get())
    #Place stats list of lists in a csv file
    with open('../../../' + output + '.csv', 'wb') as csvfile:
        stats_writer = csv.writer(csvfile, delimiter = ',')
        stats_writer.writerows(list_of_lists)
        csvfile.close()

#Check All Bams for corruption

Bad_Bams = []
for file in os.listdir(cwd) :
    if file.endswith("bam") :
        try:
            subprocess.check_output('samtools quickcheck -v ' + file, shell=True)
        except  subprocess.CalledProcessError as e:
            Bad_Bams.append(e.output)

if len(Bad_Bams) < 1 :
    print "No Corrupt Bams identified with samtools quickcheck"
else :
    print "Problem Bams detected! See printed error messages and check bad_bams.txt file"
    output_file = open("bad_bams.txt",'w')
    for i in Bad_Bams :
        print>>output_file, i
    output_file.close()


print "Novascreen complete! Please check error logs before moving on."
