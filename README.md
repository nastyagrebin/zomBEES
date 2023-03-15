OVERVIEW
- This repository contains the code and figures generated for the 20.440 project zomBEEs. The goal of the project is to evaluate whether infection of the common honeybee with the deformed wing virus (DWV) causes a change in the transcriptional pattern of genes related to hygienic behavior. Notably, since we are using very large datasets and GitHub has limits on the amount of data that can be stored, the datasets have to be downloaded locally into a data/ folder for the code to properly run. Currently, the repo contains code needed to align and process the paired-end raw RNA-seq reads from a single bee brain and a Jupyter notebook required for generating a quality score vs read length graph for that alignment.

SCRIPT CITATIONS 
- fastq_download.py
  - Blog post by Eric Liu: https://erilu.github.io/python-fastq-downloader/

DATA
- Boutin, S., Alburaki, M., Mercier, PL. et al. Differential gene expression between hygienic and non-hygienic honeybee (Apis mellifera L.) hives. BMC Genomics 16, 500 (2015). https://doi.org/10.1186/s12864-015-1714-y
  - RNA-seq was performed on brains of honeybees (Apis mellifera) from hygenic and non-hygenic colonies (8 colonies total).

- Pizzorno MC, Field K, Kobokovich AL, Martin PL, Gupta RA, Mammone R, Rovnyak D, Capaldi EA. Transcriptomic Responses of the Honey Bee Brain to 
Infection with Deformed Wing Virus. Viruses. 2021 Feb 12;13(2):287. doi: 10.3390/v13020287. PMID: 33673139; PMCID: PMC7918736.
  - RNA-seq was performed on the brains of 20 honeybees (Apis mellifera) experimentally infected with deformed wing virus (DWV)
  - Data accession: https://www.ncbi.nlm.nih.gov/sra/?term=Transcriptomic+Responses+of+the+Honey+Bee+Brain+to+Infection+with+Deformed+Wing+Virus 

- Wallberg, A. et al. A hybrid de novo genome assembly of the honeybee, Apis mellifera, with chromosome-length scaffolds. BMC Genomics 20, 275 (2019).
  - The genome of Apis mellifera was assembled using hyper long-read sequences and chromosome-length scaffolds.
  - Data accession: https://www.ncbi.nlm.nih.gov/assembly/GCF_003254395.2/



FOLDER STRUCTURE
- data: All 20 FASTQ files from the RNA-seq of DWV-infected bee brains, alignment outputs such as indices and .bam files, .csv tables of known hygiene-related genes in Apis mellifera.
  - Note: this folder is not currently in our repository! In order to run the code to generate figure 1, the user must download this repo, create a data/ folder in the root directory, download the reference genome GCF_003254395.2 file, and run the generate_figure.py script from the code/ directory. This wiill automatically pull the fastq files for all thhe RNA-seq samples from the SRA database. In the future, our scriipt will also pull the reference genome automatically from the web. To generate the final figure, the user then simply needs to run the Figure.ipynb file -- it will be automatically saved to the figures/ directory.
- code: Code for figure generation, downloading data, alignment, .sam file conversion, and read filtering.
- figures: Figure 1, submitted for pset4

INSTALLATION
- Clone the repo to your local machine
- Download the reference genome at https://www.ncbi.nlm.nih.gov/assembly/GCF_003254395.2/
- Install the following packages using conda or another high-level package manager
  -SRA-toolkit=3.0.1
  -fastx-toolkit=0.0.13
  -git-lfs (as up to date as possible)
  -bowtie2=2.4.5
  -samtools=1.9
    -This is not the most recent version! It is important to install this particular version, later versions crash the .sam conversion!
  -bokeh=2.4.3
    -This is a visualization package. You can use matplotlib as well, but the figure generated was made using Bokeh, which provides more interactivity.

