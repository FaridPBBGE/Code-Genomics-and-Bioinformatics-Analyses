#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=fastqc_maize           # job name
#SBATCH --time=04:00:00             # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=8           # CPUs (threads) per command
#SBATCH --mem=5G                    # total memory per node
#SBATCH --output=/scratch/user/farid-bge/RNAseq_pipeline_test/stdout/stdout.%j          # save stdout to file
#SBATCH --error=/scratch/user/farid-bge/RNAseq_pipeline_test/stderr/stderr.%j           # save stderr to file

module purge
module load FastQC/0.11.9-Java-11

<<README
- FASTQC homepage: http://www.bioinformatics.babraham.ac.uk/projects/fastqc
- FASTQC manual: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help
README

################################### VARIABLES ##################################
# TODO Edit these variables as needed:
########## INPUTS ##########
####pe1_1='/scratch/data/bio/GCATemplates/data/miseq/c_dubliniensis/DR34_R1.fastq.gz'
####pe1_2='/scratch/data/bio/GCATemplates/data/miseq/c_dubliniensis/DR34_R2.fastq.gz'
input_file='/scratch/user/farid-bge/RNAseq_pipeline_test/maize_Raw_RNAdata/trimmed_sample'
######## PARAMETERS ########
threads=$SLURM_CPUS_PER_TASK

########## OUTPUTS #########
output_dir='/scratch/user/farid-bge/RNAseq_pipeline_test/maize_Raw_RNAdata/fastqc_output/fastqc_after_trimming'

################################### COMMANDS ###################################
# use -o <directory> to save results to <directory> instead of directory where reads are located
#   <directory> must already exist before using -o <directory> option
# --nogroup will calculate average at each base instead of bins after the first 50 bp
# fastqc runs one thread per file; using 20 threads for 2 files does not speed up the processing
for f in $input_file/*.fq.gz
do
zcat ${f} | fastqc -t $threads -o $output_dir ${f}
done

####comment: took 26 mins for 8 files in HPRC
################################################################################
<<CITATION
- Acknowledge TAMU HPRC: https://hprc.tamu.edu/research/citations.html

- FastQC: http://www.bioinformatics.babraham.ac.uk/projects/fastqc
CITATION