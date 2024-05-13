## DNA_Mapping Pipeline
For DNA mapping, I used BWA tool. At first, I used BWA 0.7.17 version with samtools 1.17 for conversting sam files to sorted bam files in the same script file. Because of slow processing of BWA 0.7.17, I used then BWA last version, BWA-mem2 tool with samtools 1.17.
~~~
#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=bwa_mapping           # job name
#SBATCH --time=12:00:00             # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=20          # CPUs (threads) per command
#SBATCH --mem=300G                   # total memory per node
#SBATCH --output=/scratch/user/farid-bge/stdout/stdout.%j          # save stdout to file
#SBATCH --error=/scratch/user/farid-bge/stderr/stderr.%j           # save stderr to file
## For BWA-mem2
module purge
module load GCC/11.3.0
module load SAMtools/1.17
### Loading module
module load bwa-mem2/2.2.1-Linux64
################################### VARIABLES ##################################
# TODO Edit these variables as needed:
########## INPUTS ##########
###pe1_1='/scratch/data/bio/GCATemplates/e_coli/seqs/SRR10561103_1.fastq.gz'
###pe1_2='/scratch/data/bio/GCATemplates/e_coli/seqs/SRR10561103_2.fastq.gz'
###look for already indexed genome here /scratch/data/bio/genome_indexes/ucsc/
ref_genome='/scratch/user/farid-bge/DNA_mapping/ref_genome_maize/ref_index_bwa-mem2/Zmays_833_Zm-B73-REFERENCE-NAM-5.0.fa.gz'
######## PARAMETERS ########
cpus=$SLURM_CPUS_PER_TASK
readgroup='maize_sra'
library='pe'
sample='maize_mapping_test'
platform='ILLUMINA'
########## OUTPUTS #########
###output_bam="${sample}_sorted.bam"
################################### COMMANDS ###################################
# NOTE index genome only if not using already indexed genome
##if [ ! -f ${ref_genome}.bwt ]; then
##bwa index $ref_genome
##fi
### looping the bwa command
##for f in $(ls *1_paired.fq.gz | sed 's/1_paired.fq.gz//')
##do
#### For single sample running###
samp='ERR10235200'
##bwa-mem2 mem -M -t $cpus -R "@RG\tID:$readgroup\tLB:$library\tSM:$sample\tPL:$platform" $ref_genome ${samp}_1_paired.fq.gz ${samp}_2_paired.fq.gz > /scratch/user/farid-bge/DNA_mapping/mapping/bwa-mem2_mapping_result${samp}.sam 
bwa-mem2 mem -M -t $cpus -R "@RG\tID:$readgroup\tLB:$library\tSM:$sample\tPL:$platform" $ref_genome ${samp}_1_paired.fq.gz ${samp}_2_paired.fq.gz -o $TMPDIR/${samp}_aln.sam
samtools sort $TMPDIR/${samp}_aln.sam -o /scratch/user/farid-bge/DNA_mapping/mapping/bwa-mem2_mapping_result/${samp}_sorted_with_tem_dir_sam.bam  -m 5G -@ $cpus -T $TMPDIR/tmp4sort
~~~
