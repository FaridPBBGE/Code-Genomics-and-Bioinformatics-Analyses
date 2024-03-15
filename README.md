# Bioinformatics Pipelines for Postdoctoral Research @IHA-TAMU
This repository showcases bioinformatics analysis pipelines developed during my postdoctoral research at the Institute for Advancing Health Through Agriculture(IHA), Texas A&M AgriLife Research, College Station, TX.These pipelines cover the following tasks of bioinformatics analyses:
* **RNA and DNA sequencing analysis:** Mapping reads, identifying differentially expressed genes, and exploring alternative splicing events.
* **SNP calling:** Identifying and characterizing genetic variations.
* **Quality control of raw sequencing data:** Assessing data integrity and filtering low-quality reads.

These pipelines utilize scripting languages, including Bash, R, and Python. Scripts designed specifically for execution on the TAMU HPRC cluster are included, along with tools commands.

**For further details and application of these pipelines within my research, please refer to my publications:**

* [Link to Publication 1]
* [Link to Publication 2]

**Feel free to contact me for any inquiries.**

## TAMU-HPRC Job Script Structure
~~~
#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=fastqc_maize           # job name
#SBATCH --time=04:00:00             # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=8           # CPUs (threads) per command
#SBATCH --mem=5G                    # total memory per node
#SBATCH --output=stdout.%j          # save stdout to file
#SBATCH --error=stderr.%j           # save stderr to file
## These are the parameters for TAMU HPRC slurm script.
~~~
### Job submission and monitoring job
* Submit a job: **sbatch** [script_name]
* Cancel/Kill a job: **scancel** [job_id]
* Check status of a single job: **squeue** --job [job_id]
### Searching Software
To search a software in the grace cluster, this [HPRC Available Software webpage](https://hprc.tamu.edu/kb/Software/) needs to be visited and type the software name. The following codes need to be used often to search, load, and unload the software in the HPRC terminal:
* Finding most available software on Grace: **module avail**
* Searching particular software (if I know the software name): **module spider**
* Loading the software: **module load** *software name*
* To know how many softwares or modules are loaded in the current terminal: **module list**
* To remove loaded modules: **module purge**
> It is always recommended to use module purge before using another modules in the same terminal session 
## RNA_Mapping Pipeline
Before mapping sample reads to reference genome, at first, I will index reference genome of crop of interest by using STAR/2.7.9a tool (This is the update version available in TAMU HPRC when I am using the tool). 
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
#### *Code explanation*
To save the indexed genome in a desired directory, full path of that directory was given to --genomeDir command

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
### Trimming the raw fasta files using Trimmomatic
I ran Trimmomatic with loop and without loop. Scripts for both format are given below.
~~~
### With loop
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
~~~
The below one is without loop script
~~~
#!/bin/bash
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=trimmomatic      # job name
#SBATCH --time=01:00:00           # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=8           # CPUs (threads) per command
#SBATCH --mem=5G                    # total memory per node
#SBATCH --output=/scratch/user/farid-bge/RNAseq_pipeline_test/stdout/stdout.%x.%j       # save stdout to file
#SBATCH --error=/scratch/user/farid-bge/RNAseq_pipeline_test/stderr/stderr.%x.%j        # save stderr to file

module purge
module load Trimmomatic/0.39-Java-11

################################### VARIABLES ##################################
# TODO Edit these variables as needed:

########## INPUTS ##########
pe1_1='/scratch/user/farid-bge/RNAseq_pipeline_test/maize_Raw_RNAdata/maize_raw_rnaseq_datafiles/SRR8197402_Transcriptome_dynamics_of_nucellus_in_early_maize_seed_1.fastq.gz'
pe1_2='/scratch/user/farid-bge/RNAseq_pipeline_test/maize_Raw_RNAdata/maize_raw_rnaseq_datafiles/SRR8197402_Transcriptome_dynamics_of_nucellus_in_early_maize_seed_2.fastq.gz'

######## PARAMETERS ########
cpus=$SLURM_CPUS_PER_TASK
min_length=36
quality_format='-phred33'       # -phred33, -phred64    # see https://en.wikipedia.org/wiki/FASTQ_format#Encoding
adapter_file="$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa"
# available adapter files:
#   Nextera:      NexteraPE-PE.fa
#   GAII:         TruSeq2-PE.fa,   TruSeq2-SE.fa
#   HiSeq,MiSeq:  TruSeq3-PE-2.fa, TruSeq3-PE.fa, TruSeq3-SE.fa

########## OUTPUTS #########
prefix='SRR8197402'

################################### COMMANDS ###################################

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar \
PE -threads $cpus $quality_format $pe1_1 $pe1_2 \
${prefix}_1_paired.fq.gz ${prefix}_1_unpaired.fq.gz \
${prefix}_2_paired.fq.gz ${prefix}_2_unpaired.fq.gz \
ILLUMINACLIP:$adapter_file:2:30:10:2:True \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$min_length
~~~
### Mapping reads to the indexed reference genome
~~~
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
~~~
### Counting transcript
So far I used TPMCalculator and Salmon.
#### Salmon
At first, Reference genome needs to be indexed by Salmon. The below code is for indexing:
~~~
salmon index -t Zmays_493_APGv4.fa.gz -i maizev4_index
~~~
Next, 
~~~
salmon quant -i transcripts_index -l <LIBTYPE> -1 reads1.fq -2 reads2.fq --validateMappings -o transcripts_quant
~~~

