#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=vcftools       # job name
#SBATCH --time=10:00:00             # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=16           # CPUs (threads) per command
#SBATCH --mem=200G                    # total memory per node
#SBATCH --output=/scratch/user/farid-bge/stdout/stdout_vcftools.%j      # save stdout to file
#SBATCH --error=/scratch/user/farid-bge/stderr/stderr_vcftools.%j           # save stderr to file

module purge
module load GCC/13.2.0
module load VCFtools/0.1.16

## reading samples
sample_1="Zn_bottom.txt"
sample_2="Zn_top.txt"

## running vcftools for Fst
#vcftools --vcf final_AD_gwas_snp.vcf --weir-fst-pop $sample_1 --weir-fst-pop $sample_2 --fst-window-size 100000 --fst-window-step 50000 --out "100K_step50K_Anox_intol_vs_tol"

### Getting fst result chromosome-wise using for loop
for Chr in {1..5}
do
vcftools --vcf final_AD_gwas_snp.vcf --chr $Chr --weir-fst-pop $sample_1 --weir-fst-pop $sample_2 --fst-window-size 100000 --fst-window-step 50000 --out "chr_${Chr}_100K_step50K_Zn_bottom_vs_top"
done

## running for pi
#vcftools --vcf final_AD_gwas_snp.vcf --window-pi 100000 --window-pi-step 100000 --keep "$sample_1" --out "100K_Anox_intol__pi"
#vcftools --vcf final_AD_gwas_snp.vcf --window-pi 100000 --window-pi-step 50000 --keep "$sample_1" --out "100K_step50K_Anox_intol__pi"

### Extracting vcf file for each chromosome using file name in loop
#for sample_file in *.txt
#do
#base_name="${sample_file%.txt}"

#for Chr in {1..5}
#do
#vcftools --vcf final_AD_gwas_snp.vcf --chr $Chr --keep "$sample_file" --recode --out "${base_name}_chr${Chr}"

#mv "${base_name}_chr${Chr}.recode.vcf" "${base_name}_chr${Chr}.vcf"
#done
#done

