#!/usr/bin/Rscript
# edit_dist.R: inspect edit distance distribution of aligned bam files

# Iseult Jackson 11.05.2021

## Input: CSV file of counts of edit distance from aligned bam file

## Get this with cmd: samtools view ${bam} | perl -ne 'if ($_ =~ m/NM:i:(\d+)/) { print $1, "\n"}'  | sort | uniq -c | sort -k2n > ${out}.csv
## sorting final csv numerically by edit distance is important

# usage: Rscript edit_dist.R input.csv

# prints output to stdout 
args <- commandArgs(trailingOnly=TRUE)
input <- args[1]
df <- read.csv(input,header=FALSE,sep="")
colnames(df) <- c("Count","Distance")
#Calculate -
deltas <- integer(length(df$Distance)-1)
for (i in seq(2,length(df$Distance))) {
  deltas[i-1] <- (df$Count[i] - df$Count[i-1])
}
#Metric for edit distribution
neg_delt <- sum(abs(deltas[deltas < 0]))/sum(abs(deltas))

print(input)
print(neg_delt)
