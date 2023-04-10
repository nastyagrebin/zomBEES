import subprocess

# The point of this script is to generate a bowtie2 reference
# index for our bee genome
# This only needs to be done once
# The index can then be used by TopHat2 to align RNA-seq reads
# to the reference genome.
#
# The reference genome used in this project is found at https://www.ncbi.nlm.nih.gov/assembly/GCF_000002195.4/
# It is also available on the GitHub repo for this project

# IMPORTANT: run this script from the data directory!!!!

# path to the reference
ref_file = 'GCA_000002195.1_Amel_4.5_genomic.fna'

# create a new directory
cmnd = "mkdir alignments"
subprocess.call(cmnd, shell=True)


# Invoke bowtie2
cmnd = "bowtie2-build " + ref_file + " honeybee"
subprocess.call(cmnd, shell=True)