#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=star_mapping        # job name
#SBATCH --time=10:00:00             # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=8           # CPUs (threads) per command
#SBATCH --mem=40G                    # total memory per node
#SBATCH --output=/scratch/user/farid-bge/RNAseq_pipeline_test/stdout/stdout.%j      # save stdout to file
#SBATCH --error=/scratch/user/farid-bge/RNAseq_pipeline_test/stderr/stderr.%j           # save stderr to file

module purge
module load STAR/2.7.9a-GCC-11.2.0

for f in $(ls *1_paired.fq.gz | sed 's/1_paired.fq.gz//')
do
STAR --runMode alignReads \
--runThreadN 8 \
--genomeDir /scratch/user/farid-bge/RNAseq_pipeline_test/maize_genomeDir \
--readFilesIn ${f}1_paired.fq.gz ${f}2_paired.fq.gz \
--readFilesCommand zcat \
--outSAMtype BAM Unsorted \
--outReadsUnmapped Fastx \
--quantMode GeneCounts \
--sjdbGTFfile /scratch/user/farid-bge/RNAseq_pipeline_test/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.55.renamed_chr_1.gff3 \
--sjdbGTFfeatureExon CDS \
--sjdbGTFtagExonParentGene ID \
--outFileNamePrefix ${i} \
--alignIntronMax 3000
done