import subprocess
# The goal of this script is to convert the the .sam alignment
# output file into a more manageable .bam file.

# IMPORTANT: RUN THIS SCRIPT FROM THE CODE DIRECTORY

# specify the bee sample which is being converted
bee = "SRR13397491"

# use samtools command to convert the sample to a .bam file
cmnd = f"samtools view -S -b ../data/alignments/{bee}.sam > ../data/alignments/{bee}.bam"
subprocess.call(cmnd, shell=True)

# use samtools to sort the reads
cmnd = f"samtools sort ../data/alignments/{bee}.bam -o ../data/alignments/{bee}_sorted.bam"
subprocess.call(cmnd, shell=True)

# re-index the bam file for downstream viisualization
cmnd = f"samtools faidx ../data/GCA_000002195.1_Amel_4.5_genomic.fna"
subprocess.call(cmnd, shell=True)

cmnd = f"samtools index ../data/alignments/{bee}_sorted.bam"
subprocess.call(cmnd, shell=True)

