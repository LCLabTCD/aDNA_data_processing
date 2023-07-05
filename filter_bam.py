"""

filter_bam.py 

07.06.20

Lara Cassidy

usage: filter_bam.py [-h] [-mq MAPPING_QUAL] [-bp READ_BP_LENGTH]
                     [-clip SOFTCLIP] [-R REF] [-mem JAVA_MEM] [-p PARALLEL]
                     [-op OP_DUP] [-i BAMLIST] [-o OUTDIR]

===========
DESCRIPTION
===========
 
The program applies mapping quality, read length and softclipping filters, which can be switched on or off by the user.

==================
SETUP INSTRUCTIONS
==================

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

parser.add_argument("-mq","--mapping_qual",help="Mapping quality [25]",type=int,default=25)
parser.add_argument("-minbp","--min_read_bp_length",help="Read length filter [34]",type=int,default=34)
parser.add_argument("-maxbp","--max_read_bp_length",help="Read length filter [off]",type=int,default=1000)
parser.add_argument("-clip","--softclip",help="Reduce base qualities at the ends of reads [default 2bp] to a phred score of 2 ",type=int,default=2)
parser.add_argument("-p","--parallel",help="Number of samples to run in parallel - includes bwa samse/sampe but not aln [10]",type=int,default=10)

#REQUIRED NO DEFALUTS
parser.add_argument("-i","--bamlist",help="list of bams",type=str)

#Parse Arguments
args = parser.parse_args()
parallel_jobs = args.parallel
MQ = str(args.mapping_qual)
maxbp = str(args.max_read_bp_length)
clip = str(args.softclip)
minbp = str(args.min_read_bp_length)
bamlist  =  open(args.bamlist, 'r')

#Get list of bam files
all_files = [line.rstrip() for line in bamlist]
bams = [ '.'.join(line.split(".")[0:-1])  for line in all_files]

def filter(full_ID):
	if int(MQ) > 0 :
		file1 = full_ID + '.mq' + MQ 
	else :
		file1 = full_ID 
	if int(clip) > 0 :
		file2 = file1 + '.clip' + clip + 'bp'
	else :
		file2 = file1

        if int(minbp) > 25 and int(maxbp)  < 1000:
                file3 = file2 + '.' + minbp + '-' + maxbp + 'bp'   

        elif int(minbp) > 25 and int(maxbp) > 999 :
                file3 = file2 + '.' + minbp + 'bp' 

	elif int(minbp) < 26 and int(maxbp) < 1000 :
		file3 = file2 + '.25-' + maxbp + 'bp'
	else :
		file3  = file2 + '.25bp'

	try:
        	file = open(file3 + '.bam', 'r')
                print "Filtering for " + full_ID + " already complete!"
        except IOError:	
		samfile = pysam.AlignmentFile(full_ID + ".bam", "rb")			
		print "Filtering " + full_ID + "..."
		output = pysam.AlignmentFile(file3 +  '.bam', "wb", template=samfile) 

		for read in samfile :
			if read.mapping_quality >= int(MQ) and read.query_length >= int(minbp) and read.query_length <= int(maxbp):
				clipped = "#" * int(clip)
				internal_pos =  read.query_length - int(clip)
				internal_string = str(read.qual[int(clip):internal_pos])
				soft_clipped_quality_string =  clipped + internal_string + clipped
				read.qual = soft_clipped_quality_string
				output.write(read)
		output.close()
		call('samtools index ' + file3  + '.bam',shell=True)

#Processing

Parallel(n_jobs=parallel_jobs)(delayed(filter)(full_ID) for full_ID in bams)	

