import subprocess

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

# this will download the .sra files to ~/ncbi/public/sra/ (will create directory if not present)
for sra_id in sra_numbers:
    print ("Currently downloading: " + sra_id)
    prefetch = "prefetch " + sra_id
    print ("The command used was: " + prefetch)
    subprocess.call(prefetch, shell=True)

# this will extract the .sra files from above into a folder named 'fastq'
for sra_id in sra_numbers:
    print ("Generating fastq for: " + sra_id)
    fastq_dump = "fasterq-dump --outdir fastq ./" + sra_id + "/" + sra_id + ".sra"
    print ("The command used was: " + fastq_dump)
    subprocess.call(fastq_dump, shell=True)