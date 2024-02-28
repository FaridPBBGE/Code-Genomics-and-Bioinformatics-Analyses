#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=star_indexing        # job name
#SBATCH --time=04:00:00             # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=8           # CPUs (threads) per command
#SBATCH --mem=15G                    # total memory per node
#SBATCH --output=/scratch/user/farid-bge/RNA_mapping_pipeline/stdout/stdout.%j      # save stdout to file
#SBATCH --error=/scratch/user/farid-bge/RNA_mapping_pipeline/stderr/stderr.%j           # save stderr to file

module purge
module load STAR/2.7.9a-GCC-11.2.0

STAR --runMode genomeGenerate \
--runThreadN 8 \
--genomeDir /scratch/user/farid-bge/RNA_mapping_pipeline/medicago_rhizo_RNA_mapping/medi_rhizo_v2_annot_CDS_genomeDir \
--genomeFastaFiles /scratch/user/farid-bge/RNA_mapping_pipeline/medicago_rhizo_RNA_mapping/ref_genome/merged_medtr.R108__sm2011_ref_genome.fasta \
--sjdbGTFfile /scratch/user/farid-bge/RNA_mapping_pipeline/medicago_rhizo_RNA_mapping/ref_genome/merged_medtr.R108_rhizo_382_1197.gff3 \
--sjdbOverhang 100 \
--sjdbGTFtagExonParentGene ID \
--genomeSAindexNbases 10 \
--sjdbGTFfeatureExon CDS





