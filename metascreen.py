#!/usr/bin/python

"""
metascreen.py: for all your metagenomics needs

Script to process and filter unaligned reads for metagenomic analysis. 

Input: Unaligned BAM files from output of seqpro_v7.py run with -un Y flag

Run this in the unaligned directory within the seqpro_v7 directory structure.

Only uses *collapsed* reads if data is paired-end

Steps:

1. Filters unaligned reads for exact p7 index match. If -p5 Y is set, also filters for exact p5 match
2. Convert unaligned BAM to FASTQ format
3. Deduplicate FASTQ using Prinseq++ (currently only installed in Iseult's conda env prinseq-pp: need to install system-wide)
4. Merge deduplicated FASTQ to sample level
5. Taxonomic profiling with Kraken2. Uses /st1/hdd/pg/metagenomic_resources/kracken2_standard_db_mar2020/ as default. 
        ****NOTE: if you change this, you need to run bracken-build for the new database.****
6. Refining taxonomic assignments using Bracken. Default settings: readlength of 65, threshold 50 reads. 
        Threshold should be adjusted for data with higher sequencing depths (e.g. HiSeq, NovaSeq- use a 100 read threshold)
        Readlength of 100 is also possible, but not recommended for ancient data, so I haven't given the option to change it here: however, Bracken takes about half a minute to run, so this can be done after the process is complete.
        If other average readlengths are wanted, you will need to run bracken-build for these readlengths

Usage: 


metascreen.py [-h] [-db DATABASE] [-pe PAIRED] [-t THREADS]
                     [-th THRESHOLD] [-p7 P7_FISH] [-p5 P5_FISH]
                     [-j PARALLEL_JOBS]

optional arguments:
  -h, --help            show this help message and exit
  -db DATABASE, --database DATABASE
                        Kraken/Bracken Database
  -pe PAIRED, --paired PAIRED
                        Are files paired-end? Y or N
  -t THREADS, --threads THREADS
                        Number of threads to use
  -th THRESHOLD, --threshold THRESHOLD
                        Threshold number of reads for Bracken reassingments
  -p7 P7_FISH, --p7_fish P7_FISH
                        Fish for exact p7 index matches
  -p5 P5_FISH, --p5_fish P5_FISH
                        Fish for exact p5 index matches
  -j PARALLEL_JOBS, --parallel_jobs PARALLEL_JOBS
                        Number of jobs to be run in parallel (NOT Kraken)

Iseult Jackson 07.12.2021

"""

#########################################################
#########################################################
####                    SETUP                        ####
#########################################################
#########################################################


# import statements
from __future__ import division
import sys
import os 
import argparse
import pysam
import errno
from joblib import Parallel, delayed
from subprocess import call, check_output

## get cwd: where unaligned BAM files live

cwd = os.getcwd()

## parser module
parser = argparse.ArgumentParser()

## Arguments
parser.add_argument("-db","--database",help="Kraken/Bracken Database",default="/st1/hdd/pg/metagenomic_resources/kracken2_standard_db_mar2020/")
parser.add_argument("-pe","--paired",help="Are files paired-end? Y or N",default="Y")
parser.add_argument("-t","--threads",help="Number of threads to use",default=5)
parser.add_argument("-th","--threshold",help="Threshold number of reads for Bracken reassingments", default=50)
parser.add_argument("-p7","--p7_fish", help="Fish for exact p7 index matches", default="Y")
parser.add_argument("-p5","--p5_fish", help="Fish for exact p5 index matches", default="N")
parser.add_argument("-j","--parallel_jobs",help="Number of jobs to be run in parallel (NOT Kraken)", default=5)
parser.add_argument("-bracken","--bracken",help="path to bracken", default="/Software/Bracken-2.5/bracken")
parser.add_argument("-prinseq","--prinseq",help="path to prinseq++", default="/home/iseult/miniconda3/envs/prinseq-pp/bin/prinseq++")

args = parser.parse_args()

threshold = str(args.threshold)
db = str(args.database)
paired = str(args.paired)
threads = str(args.threads)
p7_fish = str(args.p7_fish)
p5_fish = str(args.p5_fish)
parallel_jobs = int(args.parallel_jobs)
path_to_bracken=str(args.bracken)
path_to_prinseq=str(args.prinseq)
## Set name of database: use this to create directory within kraken_output for these assignments
db_name = db.split('/')[-2]

## Index dictionaries for index fishing

p5_dict = {"indexing1": "TCGCAGG", "indexing2": "CTCTGCA", "indexing3": "CCTAGGT", "indexing4": "GGATCAA", "indexing5": "GCAAGAT", "indexing6": "ATGGAGA", "indexing7": "CTCGATG", "indexing8": "GCTCGAA", "indexing9": "ACCAACT", "indexing10": "CCGGTAC", "indexing11": "AACTCCG", "indexing12": "TTGAAGT", "indexing13": "ACTATCA","indexing14":"TTGGATC","indexing15":"CGACCTG","indexing16":"TAATGCG","indexing19":"GAATCTC","indexing20":"CATGCTC","indexing21":"ACGCAAC","indexing22":"GCATTGG","indexing23":"GATCTCG","indexing24":"CAATATG","indexing25":"TGACGTC","indexing26":"GATGCCA","indexing27":"CAATTAC","indexing28":"AGATAGG","indexing29":"CCGATTG","indexing30":"ATGCCGC","indexing31":"CAGTACT","indexing32":"AATAGTA","indexing33":"CATCCGG","indexing34":"TCATGGT","indexing35":"AGAACCG","indexing36":"TGGAATA","indexing37":"CAGGAGG","indexing38":"AATACCT","indexing39":"CGAATGC","indexing40":"TTCGCAA","indexing41":"AATTCAA","indexing42":"CGCGCAG","indexing43":"AAGGTCT","indexing44":"ACTGGAC","indexing45":"AGCAGGT","indexing46":"GTACCGG","indexing48":"AATGATG","indexing49":"AGTCAGA","indexing50":"AACTAGA","indexing51":"CTATGGC"}

p7_dict = {"indexing1": "TCGCAGG" , "indexing2": "CTCTGCA" , "indexing3": "CCTAGGT" , "indexing4": "GGATCAA" , "indexing5": "GCAAGAT" , "indexing6": "ATGGAGA" , "indexing7": "CTCGATG" , "indexing8": "GCTCGAA" , "indexing9": "ACCAACT" , "indexing10": "CCGGTAC" , "indexing11": "AACTCCG" , "indexing12": "TTGAAGT" , "indexing13": "ACTATCA" , "indexing14": "TTGGATC" , "indexing15": "CGACCTG" , "indexing16": "TAATGCG" , "indexing17": "AGGTACC" , "indexing18": "TGCGTCC" , "indexing19": "GAATCTC" , "indexing20": "CATGCTC" , "indexing21": "ACGCAAC" , "indexing22": "GCATTGG" , "indexing23": "GATCTCG" , "indexing24": "CAATATG" , "indexing25": "TGACGTC" , "indexing26": "GATGCCA" , "indexing27": "CAATTAC" , "indexing28": "AGATAGG" , "indexing29": "CCGATTG" , "indexing30": "ATGCCGC" , "indexing31": "CAGTACT" , "indexing32": "AATAGTA" , "indexing33": "CATCCGG" , "indexing34": "TCATGGT" , "indexing35": "AGAACCG" , "indexing36": "TGGAATA" , "indexing37": "CAGGAGG" , "indexing38": "AATACCT" , "indexing39": "CGAATGC" , "indexing40": "TTCGCAA" , "indexing41": "AATTCAA" , "indexing42": "CGCGCAG" , "indexing43": "AAGGTCT" , "indexing44": "ACTGGAC" , "indexing45": "AGCAGGT" , "indexing46": "GTACCGG" , "indexing47": "GGTCAAG" , "indexing48": "AATGATG" , "indexing49": "AGTCAGA" , "indexing50": "AACTAGA" , "indexing51": "CTATGGC" , "indexing52": "CGACGGT" , "indexing53": "AACCAAG" , "indexing54": "CGGCGTA" , "indexing55": "GCAGTCC" , "indexing56": "CTCGCGC" , "indexing57": "CTGCGAC" , "indexing58": "ACGTATG" , "indexing59": "ATACTGA" , "indexing60": "TACTTAG" , "indexing61": "AAGCTAA" , "indexing62": "GACGGCG" , "indexing63": "AGAAGAC" , "indexing64": "GTCCGGC" , "indexing65": "TCAGCTT" , "indexing66": "AGAGCGC" , "indexing67": "GCCTACG" , "indexing68": "TAATCAT" , "indexing69": "AACCTGC" , "indexing70": "GACGATT" , "indexing71": "TAGGCCG" , "indexing72": "GGCATAG" , "indexing73": "TTCAACC" , "indexing74": "TTAACTC" , "indexing75": "TAGTCTA" , "indexing76": "TGCATGA" , "indexing77": "AATAAGC" , "indexing78": "AGCCTTG" , "indexing79": "CCAACCT" , "indexing80": "GCAGAAG" , "indexing81": "AGAATTA" , "indexing82": "CAGCATC" , "indexing83": "TTCTAGG" , "indexing84": "CCTCTAG" , "indexing85": "CCGGATA" , "indexing86": "GCCGCCT" , "indexing87": "AACGACC" , "indexing88": "CCAGCGG" , "indexing89": "TAGTTCC" , "indexing90": "TGGCAAT" , "indexing91": "CGTATAT" , "indexing92": "GCTAATC" , "indexing93": "GACTTCT" , "indexing94": "GTACTAT" , "indexing95": "CGAGATC" , "indexing96": "CGCAGCC" , "indexing97": "GAGAGGC" , "indexing98": "GCTTCAG" , "indexing99": "ATATCCA" , "indexing100": "GTTATAC" , "indexing101": "CCTTAAT" , "indexing102": "CGCCAAC" , "indexing103": "TACTCGC" , "indexing104": "AGCGCCA" , "indexing105": "TAAGTAA" , "indexing106": "TTGGTCA" , "indexing107": "GTTGCAT" , "indexing108": "ATCCTCT" , "indexing109": "GGCGGTC" , "indexing110": "ATATGAT" , "indexing111": "GGTACGC" , "indexing112": "AAGAACG" , "indexing113": "CCGTTGA" , "indexing114": "AGCAATC" , "indexing115": "GCTCCGT" , "indexing116": "CAACTCT" , "indexing117": "AGACTCC" , "indexing118": "CTATCTT" , "indexing119": "AAGCAGT" , "indexing120": "GTTACCG" , "indexing121": "CCTAACG" , "indexing122": "ATCATAA" , "indexing123": "TGATAAC" , "indexing124": "TCTCCTA" , "indexing125": "CAGAGCA" , "indexing126": "CGGCTGG" , "indexing127": "CGCTATT" , "indexing128": "GATCGTC" , "indexing129": "ACGGCAG" , "indexing130": "GACCGAT" , "indexing131": "ACTTGCG" , "indexing132": "GTAAGCC" , "indexing133": "GCCATGC" , "indexing134": "ATAACGT" , "indexing135": "CGGACGT" , "indexing136": "GCGAGTA" , "indexing137": "ACGCGGA" , "indexing138": "GTCTAAT" , "indexing139": "GAAGCGT" , "indexing140": "CGGTAAG" , "indexing141": "GGTAACT" , "indexing142": "AATATAG" , "indexing143": "CCTCGCC" , "indexing144": "TTAATAG" , "indexing145": "CCGAAGC" , "indexing146": "TCGTTAT" , "indexing147": "GGCTCTG" , "indexing148": "CCAAGTC" , "indexing149": "CTTGGAA" , "indexing150": "TGAAGCT" , "indexing151": "GGTTGAC" , "indexing152": "CCGCCAT" , "indexing153": "ACCAGAG" , "indexing154": "GCTGAGA" , "indexing155": "CCTTCGC" , "indexing156": "CTGGCCT" , "indexing157": "CGCAAGG" , "indexing158": "TGAGAGA" , "indexing159": "AAGATTC" , "indexing160": "ATCGGTT" , "indexing161": "ACGAGCC" , "indexing162": "TGGATAC" , "indexing163": "ATTCCAG" , "indexing164": "ACTCATT" , "indexing165": "ACCTCGT" , "indexing166": "AGCTTAT" , "indexing167": "GCGATCT" , "indexing168": "CTCCAGT" , "indexing169": "GAACTTA" , "indexing170": "CCAATAA" , "indexing171": "AGGCGAG" , "indexing172": "CTCAGAT" , "indexing173": "AAGACGA" , "indexing174": "ACCGCTC" , "indexing175": "AACTGAC" , "indexing176": "GCAACTG" , "indexing177": "GAGTAAC" , "indexing178": "GATTAGG" , "indexing179": "AGTCGCT" , "indexing180": "CATAGAC" , "indexing181": "ACTACGG" , "indexing182": "GGCTAGC" , "indexing183": "CCATGAG" , "indexing184": "GCTGGTT" , "indexing185": "AGAGTAG" , "indexing186": "TATCAAC" , "indexing187": "ATACGCG" , "indexing188": "GACGTAC" , "indexing189": "CTACTTC" , "indexing190": "TGGTCGG" , "indexing191": "AGACCAT" , "indexing192": "TCCGAAC" , "indexing193": "GAGGTTG" , "indexing194": "TATGAGT" , "indexing195": "CTTCGTT" , "indexing196": "CCGTCCG" , "indexing197": "AACGTTA" , "indexing198": "GCATATT" , "indexing199": "ACCTTCC" , "indexing200": "TTCCGAG"}

#########################################################
#########################################################
####                    FUNCTIONS                    ####
#########################################################
#########################################################

# Create correct dir structure

def ensure_directory_exists(dir) :
    try:
        os.mkdir(dir)
    except OSError as exception :
        if exception.errno != errno.EEXIST :
            raise

# Fish for exact index matches. Output is a BAM file

def index_fish(bamfile, p7_fish,p5_fish):
        # Define outfile name
        outbam = ".".join(bamfile.split('.')[0:-2]) + '.fished.bam'
        outdir = 'filtered_reads/fished/'
        # Define p5 and p7 indices for bamfile
        index = bamfile.split('-')[4]
        #ensure formatted correctly
        if index.startswith('0'):
                index.replace("0","")
        P7 = index.split("_")[0]
        p7_check = "indexing" + str(P7)
        try:
                P5 = index.split("_")[1]
                p5_check = "indexing" + str(P5)
        except IndexError:
                P5 = ''
                p5_fish = "N"
                print('No p5 info available for this file; p5 fishing turned off')

        # If p5 fishing is on, check that data is actually paired end
        if p5_fish == "Y":
                if paired == "N":
                        sys.exit("p5 fishing turned on for SE data: please enter valid combination of parameters.")
                # Fish for both p5 and p7
                elif p7_fish == "Y":
                        # index lookup
                        p7_seq = p7_dict.get(p7_check, "P7_index_not_found")
                        if p7_seq == "P7_index_not_found":
                                sys.exit("P7 not in dictionary: please check all file names and formats")
                        p5_seq = p5_dict.get(p5_check, "P5_index_not_found")
                        if p5_seq=="P5_index_not_found":
                             sys.exit("P5 not in dictionary: please re-try with a custom config file")   
                        try:
                                f = open(outdir + outbam,'r')
                                print("Fished bam file for " + bamfile + " already exists! Moving on...")
                        except IOError:
                                # open bam file for reading
                                bfile = pysam.AlignmentFile(bamfile,'rb')
                                # create outfile (bam format)
                                bout = pysam.AlignmentFile(outdir + outbam,"wb", template = bfile)
                                for line in bfile.fetch(until_eof=True):
                                        read = str(line).split("\t")[0]
                                        if p7_seq + '+' + p5_seq in read:
                                                bout.write(line)
                                bout.close()
                                bfile.close()
                # Don't allow fishing for just p5s?
                else:
                        sys.exit("If p5_fish is on, p7_fish must also be on. Please re-run with correct parameters")
        elif p5_fish == "N" and p7_fish =="Y":
                # index lookup
                p7_seq = p7_dict.get(p7_check, "P7_index_not_found")
                if p7_seq == "P7_index_not_found":
                        sys.exit("P7 not in dictionary: please check all file names and formats")
                try:
                        f = open(outdir + outbam,'r')
                        print("Fished bam file for " + bamfile + " already exists! Moving on...")
                except IOError:
                        # get exact index matches
                        # open bam file for reading
                        bfile = pysam.AlignmentFile(bamfile,'rb')
                        # create outfile (bam format)
                        bout = pysam.AlignmentFile(outdir + outbam,"wb", template = bfile)
                        for line in bfile.fetch(until_eof=True):
                                read = str(line).split("\t")[0]
				if paired == "Y":
                                	if p7_seq + '+' in read:
                                        	bout.write(line)
				elif paired == "N":
					if read.split(":")[-1] == p7_seq:
						bout.write(line)
                        bout.close()
                        bfile.close()
        elif p5_fish =="N" and p7_fish =="N":
                print("No index fishing to be carried out: continuing...")
        else:
                sys.exit("Urecognised options for index fishing: please re-run with values Y or N for p7_fish and p5_fish")

# Convert bam to fastq

def bam_to_fq(bamfile, index_fish):
        if index_fish == "Y":
                # If index fishing is on, convert fished bam to fastq file
                inbam = ".".join(bamfile.split('.')[0:-2]) + '.fished.bam'
                outfq = ".".join(bamfile.split('.')[0:-2]) + '.fished.fastq.gz'
                fishdir = 'filtered_reads/fished/'
                inbam = fishdir + inbam
                outfq = fishdir + outfq
                try: 
                        f = open(outfq,'r')
                        print("Fished fastq for " + bamfile + " exists! Moving on...")
                except IOError:
                        pysam.bam2fq('-n','-0',outfq,inbam)
        else:
                outfq = ".".join(bamfile.split('.')[0:-2]) + '.fastq.gz'
                try: 
                        f = open(outfq,'r')
                        print("Unaligned fastq for " + bamfile + " exists! Moving on...")
                except IOError:
                        pysam.bam2fq('-n','-0', outfq,bamfile)

## Deduplication with prinseq++

def dedup(fq):
        outdir = 'filtered_reads/deduped/'
        out_prefix = ".".join(fq.split('/')[-1].split('.')[0:-2]) + ".prinseq"
        try:
                f=open(outdir + out_prefix + "_good_out.fastq.gz",'r')
                print("Fastq " + str(fq) + " already de-duped! Moving on...")
        except IOError:
                call(path_to_prinseq + " -fastq " + fq + " -out_format 0 -out_name " + outdir + out_prefix + " -derep -threads " + threads + " > " + outdir + out_prefix +".stdout 2> "+ outdir + out_prefix +".stderr", shell=True)
                # bgzip filtered fastqs: space!
                call('bgzip ' + outdir + out_prefix + "_good_out.fastq", shell=True)
                call('bgzip ' + outdir + out_prefix + "_bad_out.fastq", shell=True)

## Merge to sample level

def sample_merge(sample):
        outdir = 'filtered_reads/merged/'
        outfile = sample + '.trimmed25bp25q.fished.dedup.concat.fastq.gz'
        try:
                f = open(outdir + outfile, 'r')
                print("Sample merges for " + sample + " complete! Moving on...")
        except IOError:
                call('cat filtered_reads/deduped/' + sample + '-*good_out.fastq.gz > ' + outdir + outfile, shell=True)

## Kraken2

def kraken(sample):
        input = 'filtered_reads/merged/'  + sample + '.trimmed25bp25q.fished.dedup.concat.fastq.gz'
        output = 'tax_assignment/kraken_output/' + db_name +'/' + sample + '.trimmed25bp25q.fished.dedup.concat.kraken'
        try:
                f = open(output + '.report','r')
                print("Kraken2 already run on " + input + "! Moving on...")
        except IOError:
                call('kraken2 --use-names --threads ' + threads + ' --db ' + db + ' --report ' + output + '.report --output ' + output + '.out ' + input, shell = True)

## Bracken

def bracken(sample):
        input = 'tax_assignment/kraken_output/' + db_name +'/'+ sample + '.trimmed25bp25q.fished.dedup.concat.kraken.report'
        output = 'tax_assignment/bracken_output/'  + db_name +'/'+ sample + '.trimmed25bp25q.fished.dedup.concat.r65t' + threshold + '.bracken'
        try:
                f = open(output + '.bracken','r')
                print("Bracken already run for " + sample + "! Moving on...")
        except IOError:
                call(path_to_bracken + " -d " + db + " -i " + input + " -o " + output + ' -r 65 -t ' + threshold, shell=True)
		# move bracken reports to bracken output folder
		call('mv tax_assignment/kraken_output/' + db_name +'/'+ sample + '*bracken* tax_assignment/bracken_output/' + db_name +'/', shell=True)


#########################################################
#########################################################
####                    PROCESSING                   ####
#########################################################
#########################################################


# Ensure directories exist

dirs = ['filtered_reads','filtered_reads/fished','filtered_reads/deduped','filtered_reads/merged','tax_assignment','tax_assignment/kraken_output','tax_assignment/kraken_output/'+db_name,'tax_assignment/bracken_output','tax_assignment/bracken_output/' + db_name]

for d in dirs:
        ensure_directory_exists(d)

# Get list of files to be processed in the current directory
# If data is paired end, this should be collapsed data only *collapsed.unaligned.bam
# If data is single end, this is all data(*unaligned.bam)
if paired == "Y":
        bams = ['.'.join(file.split('.'))for file in os.listdir(cwd) if file.endswith("collapsed.unaligned.bam")]
elif paired == "N":
        bams = ['.'.join(file.split('.'))for file in os.listdir(cwd) if file.endswith("unaligned.bam")]
else:
        sys.exit('Please enter valid value for --paired, values: Y,N')

# Sample list 
samples = sorted(set([line.split('-')[0] for line in bams]))

print(samples)

# Set single index fish variable to decide how to process reads
if p7_fish == "Y":
        index_fish =="Y"
elif p5_fish =="Y":
        index_fish == "Y"
else:
        index_fish == "N"

if index_fish == "Y":
        # Index fishing
        Parallel(n_jobs=parallel_jobs)(delayed(index_fish)(bam)for bam in bams)

# for all bams, need to convert to fastq format
Parallel(n_jobs=parallel_jobs)(delayed(bam_to_fq)(bam, index_fish) for bam in bams)

# Dedup first fastqs: only do one at a time because multi-threaded
if index_fish == "Y":
        fished_fastqs = []
        for bamfile in bams:
                 outfq = ".".join(bamfile.split('.')[0:-2]) + '.fished.fastq.gz'
                 fq = 'filtered_reads/fished/' + outfq
                 fished_fastqs.append(fq)
        
        # Dedup
        Parallel(n_jobs=1)(delayed(dedup)(fq)for fq in fished_fastqs)
else:
        first_fqs = []
        for bamfile in bams:
                 fq = ".".join(bamfile.split('.')[0:-2]) + '.fastq.gz'
                 first_fqs.append(fq)
        # Dedup
        Parallel(n_jobs=1)(delayed(dedup)(fq)for fq in first_fqs)

# Merge
Parallel(n_jobs=parallel_jobs)(delayed(sample_merge)(sample) for sample in samples)



#########################################################
#########################################################
####                    PROFILING                    ####
#########################################################
#########################################################

# Run kraken: one by one, not in parallel
Parallel(n_jobs = 1)(delayed(kraken)(sample)for sample in samples)

# Run Bracken in parallel

Parallel(n_jobs = parallel_jobs)(delayed(bracken)(sample) for sample in samples)

## Done! 

print("Metagenomics processing finished! Please delete any intermediate files and check logs before moving onto the next steps.")
