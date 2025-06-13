#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=gatk           # job name
#SBATCH --time=40:00:00             # max job run time dd-hh:mm:ss
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

### Step-3: Base Quality recalibration
## Before running BaseRecalibrator function for base quality recalibration, dictionary and index files needs to be created for Ref. genome fasta file.
## Creating the fasta seq. dictionary file
#gatk CreateSequenceDictionary \
#-R /scratch/user/farid-bge/maize_ref_genome_index/Maize_ref_genome/maize_v5_Ref_genome/assembly/Zmays_833_Zm-B73-REFERENCE-NAM-5.0.fa \
#-O /scratch/user/farid-bge/maize_ref_genome_index/Maize_ref_genome/maize_v5_Ref_genome/assembly/Zmays_833_Zm-B73-REFERENCE-NAM-5.0.dict
## org.broadinstitute.hellbender.tools.IndexFeatureFile done. Elapsed time: 1.20 minutes.


### Creating index file for known SNP file in vcf format
#gatk IndexFeatureFile \
#-I /scratch/user/farid-bge/maize_ref_genome_index/Maize_ref_genome/maize_v5_Ref_genome/annotation/Zmays_known_SNP_Ref_v5.vcf


#####Base quality recalibration has two steps to complete.
#### Building the model-Step-1

#for file in *_sort.dup.bam; do
#base_name="${file%%_sort.dup.bam}"

#gatk --java-options "-Xmx200g" BaseRecalibrator \
#-I "$file" \
#-R /scratch/user/farid-bge/maize_ref_genome_index/Maize_ref_genome/maize_v5_Ref_genome/assembly/Zmays_833_Zm-B73-REFERENCE-NAM-5.0.fa \
#--known-sites /scratch/user/farid-bge/maize_ref_genome_index/Maize_ref_genome/maize_v5_Ref_genome/annotation/Zmays_known_SNP_Ref_v5.vcf \
#-O "${base_name}_recal_data.table"
#done
 ##Elapsed time: 143.64 minutes.
 
#### Apply the model to adjust the base quality scores-Step-2
 
#for file1 in *_sort.dup.bam;do 
 #base_name="${file1%%_sort.dup.bam}"
 #file2="${base_name}_recal_data.table"
  
 #gatk --java-options "-Xmx200g" ApplyBQSR \
 #-I "$file1" \
 #-R /scratch/user/farid-bge/maize_ref_genome_index/Maize_ref_genome/maize_v5_Ref_genome/assembly/Zmays_833_Zm-B73-REFERENCE-NAM-5.0.fa \
 #--bqsr-recal-file "$file2" \
 #-O "${base_name}_sort.dup.bqsr.bam"
 #done
 
## Elapsed time: 164.73 minutes.

##### For Variant Calling, there are multiple steps to be performed
### step-1:Applying HaplotypeCaller
#gatk --java-options "-Xmx200g" HaplotypeCaller \
#-R /scratch/user/farid-bge/maize_ref_genome_index/Maize_ref_genome/maize_v5_Ref_genome/assembly/Zmays_833_Zm-B73-REFERENCE-NAM-5.0.fa \
#-I ERR10286466_sort.dup.bqsr.bam \
#-O ERR10286466.g.vcf.gz \
#-ERC GVCF
##org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCaller done. Elapsed time: 3,539.07 minutes/58 hrs for just one sample whole genome.

### Step-2: Applying CombineGVCFs 
#gatk --java-options "-Xmx200g" CombineGVCFs \
#-R /scratch/user/farid-bge/maize_ref_genome_index/Maize_ref_genome/maize_v5_Ref_genome/assembly/Zmays_833_Zm-B73-REFERENCE-NAM-5.0.fa \
#-V ERR10235200.g.vcf.gz \
#-V ERR10235203.g.vcf.gz \
#-V ERR10238125.g.vcf.gz \
#-V ERR10286466.g.vcf.gz \
#-L 2 \
#-O ERR200_203_125_466_chr2.g.vcf.gz
## Elapsed time: 22.15 minutes for single chromosome.
### Note: There is another alternative to CombineGVCFs is GenomicsDBImport command, which is suitable for large sample numbers. For few samples, CombineGVCFs is suitable.

### Step-3: Applying GenotypeGVCFs
gatk --java-options "-Xmx200g" GenotypeGVCFs \
-R /scratch/user/farid-bge/maize_ref_genome_index/Maize_ref_genome/maize_v5_Ref_genome/assembly/Zmays_833_Zm-B73-REFERENCE-NAM-5.0.fa \
-V ERR200_203_125_466_chr1.g.vcf.gz \
-L 1 \
-O ERR200_203_125_466_chr1.vcf.gz


