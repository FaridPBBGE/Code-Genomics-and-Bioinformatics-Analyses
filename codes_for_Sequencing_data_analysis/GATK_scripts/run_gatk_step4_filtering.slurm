#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=gatk_step4           # job name
#SBATCH --time=48:00:00             # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=20         # CPUs (threads) per command
#SBATCH --mem=300G                   # total memory per node
#SBATCH --output=/scratch/user/farid-bge/stdout/stdout_gatk_base_recalib.%j          # save stdout to file
#SBATCH --error=/scratch/user/farid-bge/stderr/stderr_gatk_base_recalib.%j           # save stderr to file


### Loading modules
module purge
### Loading GATK module
module load GCCcore/12.3.0
module load GATK/4.5.0.0-Java-17

### Step-4: Filter and prepare variants for downstream analysis. For this, there are several methods/protocols for filtering vcf files
###in gatk. The Variant Quality Score Recalibration or the VQSR strategy is not suitable for all organinsms because of required for so many known information. 
### So, here, I have used Convolutional Neural Net based protocol to annotate and filter the VCF file. Based on genome and available informatino, suitable method needs to be selected.

## Method name: Convolutional Neural Net. It has two steps.
## step_1: CNNScoreVariants
gatk --java-options "-Xmx200g" CNNScoreVariants \
-R /scratch/user/farid-bge/maize_ref_genome_index/Maize_ref_genome/maize_v5_Ref_genome/assembly/Zmays_833_Zm-B73-REFERENCE-NAM-5.0.fa \
-V ERR200_203_125_466_chr1.vcf.gz \
-O ERR200_203_125_466_chr1.cnns.vcf





