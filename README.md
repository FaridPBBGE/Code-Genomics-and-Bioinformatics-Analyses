# Postdoc_IHA-CS
## TAMU-HPRC Job Script Structure
~~~
#!/bin/bash
## Necessary Job Specifications
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=fastqc_maize           # job name
#SBATCH --time=04:00:00             # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=8           # CPUs (threads) per command
#SBATCH --mem=5G                    # total memory per node
#SBATCH --output=/scratch/user/farid-bge/RNAseq_pipeline_test/stdout/stdout.%j          # save stdout to file
#SBATCH --error=/scratch/user/farid-bge/RNAseq_pipeline_test/stderr/stderr.%j           # save stderr to file

## Loading modules
module purge
module load FastQC/0.11.9-Java-11
~~~
### Job submission and monitoring job
* Submit a job: **sbatch** [script_name]
* Cancel/Kill a job: **scancel** [job_id]
* Check status of a single job: **squeue** --job [job_id]
### Searching Software
To search a software in the grace cluster, this [HPRC Available Software webpage](https://hprc.tamu.edu/kb/Software/) needs to be visited and type the software name. The following codes need to be used often to search, load, and unload the software in the cluster terminal:
* Finding most available software on Grace: **module avail**
* Searching particular software (if I know the software name): **module spider**
* Loading the software: **module load** *software name*
* To know how many softwares or modules are loaded in the current terminal: **module list**
* To remove loaded modules: **module purge**
> It is always recommended to use module purge before using another modules in the same terminal session 
## RNA_Mapping Pipeline
### Reference genome indexing
~~~
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
~~~
### Quality check of fasta files
~~~
#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=fastqc_maize           # job name
#SBATCH --time=04:00:00             # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=8           # CPUs (threads) per command
#SBATCH --mem=5G                    # total memory per node
#SBATCH --output=/scratch/user/farid-bge/RNAseq_pipeline_test/stdout/stdout.%j          # save stdout to file
#SBATCH --error=/scratch/user/farid-bge/RNAseq_pipeline_test/stderr/stderr.%j           # save stderr to file

module purge
module load FastQC/0.11.9-Java-11
################################### VARIABLES ##################################
# TODO Edit these variables as needed:
########## INPUTS ##########
####pe1_1='/scratch/data/bio/GCATemplates/data/miseq/c_dubliniensis/DR34_R1.fastq.gz'
####pe1_2='/scratch/data/bio/GCATemplates/data/miseq/c_dubliniensis/DR34_R2.fastq.gz'
input_file='/scratch/user/farid-bge/RNAseq_pipeline_test/maize_Raw_RNAdata/trimmed_sample'
######## PARAMETERS ########
threads=$SLURM_CPUS_PER_TASK
########## OUTPUTS #########
output_dir='/scratch/user/farid-bge/RNAseq_pipeline_test/maize_Raw_RNAdata/fastqc_output/fastqc_after_trimming'

################################### COMMANDS ###################################
# use -o <directory> to save results to <directory> instead of directory where reads are located
#   <directory> must already exist before using -o <directory> option
# --nogroup will calculate average at each base instead of bins after the first 50 bp
# fastqc runs one thread per file; using 20 threads for 2 files does not speed up the processing
for f in $input_file/*.fq.gz
do
zcat ${f} | fastqc -t $threads -o $output_dir ${f}
done
####comment: took 26 mins for 8 files in HPRC
~~~
### Viewing MultiQC fasta files together
~~~
## load the module where I like to run the command
module load MultiQC/1.14-foss-2022b
##un the command in the directory all fastqc result files are located and add full path to save the result 
multiqc /scratch/user/farid-bge/RNAseq_pipeline_test/maize_Raw_RNAdata/fastqc_output/fastqc_after_trimming
~~~
### 
