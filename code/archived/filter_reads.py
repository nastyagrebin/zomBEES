import glob
import subprocess
# The goal of this script is to align fastq reads to 
# a reference genome using TopHat after sorting reads
# based on the quality score

# all the numbers for the 20 bee brain samples
sra_numbers = [
    "SRR13397490", 
    "SRR13397491",
    "SRR13397492", 
    "SRR13397493", 
    "SRR13397494", 
    "SRR13397495", 
    "SRR13397496", 
    "SRR13397497", 
    "SRR13397498", 
    "SRR13397499", 
    "SRR13397500",
    "SRR13397501",
    "SRR13397502",
    "SRR13397503",
    "SRR13397504",
    "SRR13397505",
    "SRR13397506",
    "SRR13397507",
    "SRR13397508",
    "SRR13397509",
    ]

# all the fastq files
path = '../data/fastq/'
files = glob.glob(path + '*.fastq')
outpath = '../data/fastq_filtered/'

# use fastx_quality_filter to remove low-quality reads
# for now, before we use clusters or parallelization, only for one bee brain
for file in ['../data/fastq/SRR13397490_1.fastq', '../data/fastq/SRR13397490_2.fastq']:
    # filter and save the filtered reads
    outfile = outpath + file[13+1:]
    print ("Currently quality-filtering: " + file)
    fastq_quality_filter = f"fastq_quality_filter -v -q 20 -p 95 -i {file} -o {outfile}"
    print ("The command used was: " + fastq_quality_filter)
    subprocess.call(fastq_quality_filter, shell=True)
    
    # gzip the unfiltered file to save memory
    zip_cmnd = "gzip " + file
    print("Compressing " + file)
    subprocess.call(zip_cmnd, shell=True)
    