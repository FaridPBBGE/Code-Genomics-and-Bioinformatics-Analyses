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

#### *Notes*
* To save the indexed genome in a desired directory, full path of that directory was given to --genomeDir command.
  

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
#### *Notes*
* The full path to the directory where sample's fastq data are located need to be given in the "input_file" variable
* Also, the output of this script can be saved to a specific directory by giving full path to the desired directory in the "output_dir" variable.
### Viewing MultiQC fasta files together
~~~
## load the module in the directory where the fastqc results files are located
module load MultiQC/1.14-foss-2022b
## run the command in the directory all fastqc result files are located and add full path to save the result 
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
### Getting count summary result of all mapping files
* At first, a bash script file needs to be made. For that, a txt file using **nano** text editor or **">"** or **touch** <file name>.sh will be created.
* The following bash code will be copied and pasted into the newly created file.
~~~
for f in $(grep -l 'Uniquely mapped reads number' *Log.final.out); do 
export percent=$(sed -n '/Uniquely mapped reads %/p' $f |cut -f2)
export num=$(sed -n '/Uniquely mapped reads number/p' $f| cut -f2)
export file=$(paste <(echo "$f") <(echo "$num") <(echo "$percent"))
echo "$file" >> compiled_STAR_stats.txt
done
~~~
* The file will be saved and closed
* If the script does not run for execute permission issue, the following command needs to be run to change the file permission.
~~~
chmod u+x <file name>.sh
~~~
* The script written here needs to be placed in the same directory where all the STAR output files are located; otherwise this script will not find the files.
### Counting transcript
So far I used TPMCalculator and Salmon.
#### Salmon
At first, Reference genome needs to be indexed by Salmon. The below code is for indexing:
~~~
salmon index -t Zmays_493_APGv4.fa.gz -i maizev4_index
~~~
Next, 
~~~
#### Loading modules and software
module purge
### Loading GCC as requirements to run Salmon in HPRC
module load GCC/11.3.0
### Loading Salmon :
module load Salmon/1.10.1

### setting parameters
cpus=$SLURM_CPUS_PER_TASK
##file1="/scratch/user/farid-bge/Transcripts_counting/maize/samples/trimmed_samples/SRR8197457_1_paired.fq.gz"
##file2="/scratch/user/farid-bge/Transcripts_counting/maize/samples/trimmed_samples/SRR8197457_2_paired.fq.gz"

##samp='SRR8197457' ## For single sample running in case need
file='/scratch/user/farid-bge/RNA_mapping_pipeline/maize_Ref_gen_Raw_data/trimmed_sample'

### Transcripts quantification by Salmon with loop to run all samples together
for samp in $(ls $file/*1_paired.fq.gz | sed 's/1_paired.fq.gz//') 
do
salmon quant -i /scratch/user/farid-bge/Transcripts_counting/maize/Ref_genome/maizev4_index -l A -1 ${samp}1_paired.fq.gz -2 ${samp}2_paired.fq.gz -p $cpus --validateMappings -o /scratch/user/farid-bge/Transcripts_counting/maize/counting/${samp}quant
done
~~~
#### TPMCalculator
