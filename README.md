# cs123-final-project

### The code folder contains the python code and a folder for the setup scripts.

# alignment.py
# -- file for the alignment, traceback, and genotyping functions

# file_io.py
# -- file that handles many input and output functions used
#    during map-reduce.

# mr_job.py
# -- contains the MRAlignment class used for submitting to EMR

# annovar.sh
# -- contains the script for running ANNOVAR in order to annotate
#    the SNPs with disease risk factors.
# -- ANNOVAR program can be downloaded from 
#    http://annovar.openbioinformatics.org/en/latest/user-guide/download/

# align.c
# -- originally we wrote the alignment algorithm in C before switching
#    to python due to complications in the input data format. 

### The testing_data folder contains small sample data sets that were
### used to test the functions.

### The results folder contains the results from the data obtained through
### BWA alignment and SAMtools. 

# SRR710115.flt.hq.vcf
# - contains the genotyping alignment and genotyping results obtained through
#   BWA alignment and SAMtools. **Note that these are not the results generated
#   by our own functions.**

# SRR710115.annotation.txt
# - contains the annotated SNPs obtained through our ANNOVAR script
# - using SRR710115.flt.hq.vcf as an input.

# plots folder
# -contains the plots obtained by analyzing the SNPs with plots.py code

### setup_scripts

# download_ref_seq.sh
# - script to download and unzip reference sequence chromosome files
# from NCBI.

# setup.sh
# - script used to set up environment for EC2 instance.

# transfer_files.sh
# - script used to transfer files from Virtual Box to EC2 instance.












