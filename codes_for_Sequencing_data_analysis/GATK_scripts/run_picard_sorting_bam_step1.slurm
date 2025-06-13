#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=picard           # job name
#SBATCH --time=24:00:00             # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=20         # CPUs (threads) per command
#SBATCH --mem=300G                   # total memory per node
#SBATCH --output=/scratch/user/farid-bge/stdout/stdout_picard.%j          # save stdout to file
#SBATCH --error=/scratch/user/farid-bge/stderr/stderr_picard.%j           # save stderr to file


### Loading modules
module purge
### Loading module
module load picard/2.25.1-Java-11

### Step-1: Sorting the BAM files
# Sorting bam files in tamu HPRC (for single sample)
java -Xmx80g -jar $EBROOTPICARD/picard.jar SortSam \
I=ERR10235200_bwa_mem2_with_tem_dir.bam \
O=ERR10235200_sort.bam \
VALIDATION_STRINGENCY=LENIENT \
SORT_ORDER=coordinate \
MAX_RECORDS_IN_RAM=3000000 \
CREATE_INDEX=True 

## Total Elapsed time for 59GB file size sample genome: 233.59 minutes.

# Sorting bam files in tamu HPRC (for looping for multiple samples)

for file in *_tem_dir.bam; do
base_name="${file%%_bwa_mem2_with_tem_dir.bam}"

java -Xmx150g -jar $EBROOTPICARD/picard.jar SortSam \
I="$file" \
O="${base_name}_sort.bam" \
VALIDATION_STRINGENCY=LENIENT \
SORT_ORDER=coordinate \
MAX_RECORDS_IN_RAM=3000000 \
CREATE_INDEX=True 
done