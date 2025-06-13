#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=picard_mark_duplicate           # job name
#SBATCH --time=30:00:00             # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=20         # CPUs (threads) per command
#SBATCH --mem=300G                   # total memory per node
#SBATCH --output=/scratch/user/farid-bge/stdout/stdout_picard_m_dup.%j          # save stdout to file
#SBATCH --error=/scratch/user/farid-bge/stderr/stderr_picard_m_dup.%j           # save stderr to file


### Loading modules
module purge
### Loading module
module load picard/2.25.1-Java-11


### Step-2: Mark duplicate reads
for file in *_sort.bam; do
base_name="${file%%_sort.bam}"

java -Xmx150g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
I="$file" \
O="${base_name}_sort.dup.bam" \
METRICS_FILE=marked_dup_metrics.txt
done

## picard.sam.markduplicates.MarkDuplicates done. Elapsed time: 308.98 minutes. _file size is 60 GB

### Evaluating processed bqsr bam file by observing insert size before forwarding variant calling (not confirmed yet)
#java -Xmx100g -Dpicard.useLegacyParser=false -jar $EBROOTPICARD/picard.jar CollectMultipleMetrics \
--REFERENCE_SEQUENCE /scratch/user/farid-bge/maize_ref_genome_index/Maize_ref_genome/maize_v5_Ref_genome/assembly/Zmays_833_Zm-B73-REFERENCE-NAM-5.0.fa \
--INPUT ERR10235200_sort.dup.bqsr.bam \
--OUTPUT ERR10235200_sort.dup.bqsr.bam.collectmultiplemetrics \
--EXTRA_ARGUMENT "CollectAlignmentSummaryMetrics::--CREATE_HISTOGRAM false" \
--EXTRA_ARGUMENT "CollectInsertSizeMetrics::--CREATE_HISTOGRAM false" \
--EXTRA_ARGUMENT "QualityScoreDistribution::--CREATE_HISTOGRAM false" \
--EXTRA_ARGUMENT "MeanQualityByCycle::--CREATE_HISTOGRAM false" \
--EXTRA_ARGUMENT "CollectBaseDistributionByCycle::--CREATE_HISTOGRAM false" \
--EXTRA_ARGUMENT "CollectGcBiasMetrics::--CREATE_HISTOGRAM false" \
--EXTRA_ARGUMENT "RnaSeqMetrics::--CREATE_HISTOGRAM false"