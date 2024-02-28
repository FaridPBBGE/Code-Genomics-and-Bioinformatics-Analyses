#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=trimmomatic      # job name
#SBATCH --time=01:00:00           # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=8           # CPUs (threads) per command
#SBATCH --mem=5G                    # total memory per node
#SBATCH --output=/scratch/user/farid-bge/RNAseq_pipeline_test/stdout/stdout.%x.%j       # save stdout to file
#SBATCH --error=/scratch/user/farid-bge/RNAseq_pipeline_test/stderr/stderr.%x.%j        # save stderr to file

module purge
module load Trimmomatic/0.39-Java-11

<<README
- Trimmomatic manual:
http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
README

################################### VARIABLES ##################################
# TODO Edit these variables as needed:

########## INPUTS ##########
pe1_1='/scratch/user/farid-bge/RNAseq_pipeline_test/maize_Raw_RNAdata/maize_raw_rnaseq_datafiles/SRR8197402_Transcriptome_dynamics_of_nucellus_in_early_maize_seed_1.fastq.gz'
pe1_2='/scratch/user/farid-bge/RNAseq_pipeline_test/maize_Raw_RNAdata/maize_raw_rnaseq_datafiles/SRR8197402_Transcriptome_dynamics_of_nucellus_in_early_maize_seed_2.fastq.gz'

######## PARAMETERS ########
cpus=$SLURM_CPUS_PER_TASK
min_length=36
quality_format='-phred33'       # -phred33, -phred64    # see https://en.wikipedia.org/wiki/FASTQ_format#Encoding
adapter_file="$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa"
# available adapter files:
#   Nextera:      NexteraPE-PE.fa
#   GAII:         TruSeq2-PE.fa,   TruSeq2-SE.fa
#   HiSeq,MiSeq:  TruSeq3-PE-2.fa, TruSeq3-PE.fa, TruSeq3-SE.fa

########## OUTPUTS #########
prefix='SRR8197402'

################################### COMMANDS ###################################

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar \
PE -threads $cpus $quality_format $pe1_1 $pe1_2 \
${prefix}_1_paired.fq.gz ${prefix}_1_unpaired.fq.gz \
${prefix}_2_paired.fq.gz ${prefix}_2_unpaired.fq.gz \
ILLUMINACLIP:$adapter_file:2:30:10:2:True \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$min_length

################################################################################
<<CITATIONS
- Acknowledge TAMU HPRC: https://hprc.tamu.edu/research/citations.html

- Trimmomatic:
Anthony M. Bolger1,2, Marc Lohse1 and Bjoern Usadel. Trimmomatic: A flexible trimmer for Illumina Sequence Data.
Bioinformatics. 2014 Aug 1;30(15):2114-20. doi: 10.1093/bioinformatics/btu170.
CITATIONS