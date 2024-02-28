#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=trimmomatic      # job name
#SBATCH --time=05:00:00           # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=8           # CPUs (threads) per command
#SBATCH --mem=5G                    # total memory per node
#SBATCH --output=/scratch/user/farid-bge/RNAseq_pipeline_test/stdout/stdout.%x.%j       # save stdout to file
#SBATCH --error=/scratch/user/farid-bge/RNAseq_pipeline_test/stderr/stderr.%x.%j        # save stderr to file

module purge
module load Trimmomatic/0.39-Java-11

# Define an array of sample prefixes
samples=("SRR8197399" "SRR8197400" "SRR8197401" "SRR8197402")

# Set common parameters
cpus=$SLURM_CPUS_PER_TASK
min_length=36
quality_format='-phred33'
adapter_file="$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa"

# Loop through the samples and process each one
for prefix in "${samples[@]}"
do
file_1="/scratch/user/farid-bge/RNAseq_pipeline_test/maize_Raw_RNAdata/trimmotic_loop_test/${prefix}_Transcriptome_dynamics_of_nucellus_in_early_maize_seed_1.fastq.gz"
file_2="/scratch/user/farid-bge/RNAseq_pipeline_test/maize_Raw_RNAdata/trimmotic_loop_test/${prefix}_Transcriptome_dynamics_of_nucellus_in_early_maize_seed_2.fastq.gz"

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar \
PE -threads $cpus $quality_format $file_1 $file_2 \
${prefix}_1_paired.fq.gz ${prefix}_1_unpaired.fq.gz \
${prefix}_2_paired.fq.gz ${prefix}_2_unpaired.fq.gz \
ILLUMINACLIP:$adapter_file:2:30:10:2:True \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$min_length
done
