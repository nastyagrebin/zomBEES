#####################################################
# This script is meant to generate the first figure #
# submitted for pset4. However, we do not recommend #
# running it as it takes a while to unpack the raw  #
# fastq files. Instead, use the .bam file generated #
# by this script and uploaded to the GitHub repo to #
# reproduce it.                                     #
# Run this whole script from the code directory.    #
#####################################################
import subprocess

# define the bee for which the script iis run
BEE = "SRR13397491"

# download the .sra file to ~/ncbi/public/sra/ (will create directory if not present)
print("Prefetching SRR file")
prefetch = "prefetch " + BEE
subprocess.call(prefetch, shell=True)

# extract the .sra file from above into a folder named 'fastq'
# this uses the fasterq-dump command
print("Extracting fastq from sra")
extract = "fasterq-dump --outdir fastq ./" + BEE + "/" + BEE + ".sra"
subprocess.call(extract, shell=True)

# create a bowtie index for aligning using a reference genome in the data folder
# Note that this only needs to be done once. If the directory already exists,
# comment out thhis entire block of code.
ref_genome = 'GCA_000002195.1_Amel_4.5_genomic.fna'
# make a temporary directory
cmnd = "mkdir ../data/alignments"
subprocess.call(cmnd, shell=True)
# invoke bowtie2
cmnd = "bowtie2-build " + ref_genome + " honeybee"
subprocess.call(cmnd, shell=True)

# align the bee data to the reference index
print(f"Aligning bee {BEE}")
# pull out fwd and rev read directions
fwd_file = f'../fastq/{BEE}_1.fastq'
rev_file = f'../fastq/{BEE}_2.fastq'

# -1 the forward filtered read file
# -2 the reverse filtered read file
# --no-unal ignore unaligned reads
# -p 8 use no more than 8 cores
# -S name of the .sam file to save the results to
cmnd = f"bowtie2 -x honeybee -1 {fwd_file} -2 {rev_file} --no-unal -p 8 -S {BEE}.sam"
subprocess.call(cmnd, shell=True)

# use samtools command to convert the sample to a .bam file
cmnd = f"samtools view -S -b ../data/alignments/{BEE}.sam > ../data/alignments/{BEE}.bam"
subprocess.call(cmnd, shell=True)

# use samtools to sort the reads
cmnd = f"samtools sort ../data/alignments/{BEE}.bam -o ../data/alignments/{BEE}_sorted.bam"
subprocess.call(cmnd, shell=True)

# re-index the bam file for downstream viisualization
cmnd = f"samtools faidx ../data/GCA_000002195.1_Amel_4.5_genomic.fna"
subprocess.call(cmnd, shell=True)
cmnd = f"samtools index ../data/alignments/{BEE}_sorted.bam"
subprocess.call(cmnd, shell=True)
