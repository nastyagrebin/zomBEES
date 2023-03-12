OVERVIEW
- Repo contains all data analysis and figure generating code for our 20.440 project

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

FOLDER STRUCTURE
- data: All 20 FASTQ files from the RNA-seq of DWV-infected bee brains
- code: Code for figure generation 

INSTALLATION
- Download SRA toolkit (https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit)
  - Our machines used 3.0.1 release for Mac OS X.
