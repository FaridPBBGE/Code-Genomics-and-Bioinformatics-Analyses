//Enable DSL2 syntax
nextflow.enable.dsl= 2

params.reads = "$baseDir/data/ggal/*_{1,2}.fq"
params.transcriptome = "$baseDir/data/ggal/transcriptome.fa"
params.genomeDir ="$baseDir/data/genomeDir"
params.gtf ="$baseDir/data/file.gff3"
params.outdir = "results"

// Workflow: Connecting the Steps to be used
workflow {

// Creating a channel that pair the files and extracts the smaple ID/name
read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

// Quality Control and Trimming
TRIMMOMATIC(read_pairs_ch)
FASTQC(TRIMMOMATIC.out.trimmed_reads)

// Indexing the transcriptome for Salmon
SALMON_INDEX(params.transcriptome)

// Dual-Analysis path
// A. Transcriptome Quantification by Salmon
SALMON_QUANT(SALMON_INDEX.out, TRIMMOMATIC.out.trimmed_reads)
// B. Genomic Alignment by STAR
ALIGN_STAR(TRIMMOMATIC.out.trimmed_reads, params.genomeDir, params.gtf)

// Collecting all outputs from fastqc and salmon and send them to MultiQC
MULTIQC(
FASTQC.out.collect(),
SALMON_QUANT.out.collect(),
ALIGN_STAR.out.log.collect()
)
}

// -----Process-----: Trimming
process TRIMMOMATIC {
    tag "Trimming $sample_id"
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    input:
    tuple val(sample_id), path(reads)
    output:
    tuple val(sample_id), path("${sample_id}_*{1,2}_paired.fq.gz"), emit: trimmed_reads
    script:
    """
    trimmomatic PE -phred33 ${reads[0]} ${reads[1]} \\
    ${sample_id}_1_paired.fq.gz ${sample_id}_1_unpaired.fq.gz \\
    ${sample_id}_2_paired.fq.gz ${sample_id}_2_unpaired.fq.gz \\
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

// Process: Quality check
process FASTQC {
tag "QC on $sample_id"
publishDir "${params.outdir}/fastqc", mode: 'copy'

input:
tuple val(sample_id), path(reads)

output:
path "*.{html,zip}"

script:
"""
fastqc ${reads[0]} ${reads[1]}
"""
}

// Process: Salmon Indexing

process SALMON_INDEX {
tag "Indexing $transcriptome.simpleName"

input:
path transcriptome

output:
path 'salmon_index'

script:
"""
salmon index --threads $task.cpus -t $transcriptome -i salmon_index
"""
}

//Process: Quantification

process SALMON_QUANT {
tag "Quantifying $sample_id"
publishDir "${params.outdir}/quant", mode: 'copy'

input:
path index_dir
tuple val(sample_id), path(reads)

output:
path "${sample_id}_quant"

script:
"""
salmon quant -i $index_dir -l A -1 ${reads[0]} -2 ${reads[1]} -o ${sample_id}_quant
"""
}

// Process: Mapping by STAR

process ALIGN_STAR {
    tag "STAR Align: $sample_id"
    publishDir "${params.outdir}/star_alignment", mode: 'copy'
    input:
    tuple val(sample_id), path(reads)
    path genomeDir
    path gtf
    output:
    path "${sample_id}*"
    path "${sample_id}Log.final.out", emit: log
    script:
    """
    STAR --runMode alignReads --runThreadN 8 --genomeDir $genomeDir \\
         --readFilesIn ${reads[0]} ${reads[1]} --readFilesCommand zcat \\
         --outSAMtype BAM Unsorted --quantMode GeneCounts --sjdbGTFfile $gtf \\
         --outFileNamePrefix ${sample_id}
    """
}

//Process: MultiQC visualization 

process MULTIQC {
tag "MultiQC:Summarizing all results"
publishDir "${params.outdir}/multiqc", mode: 'copy'

input:
path fastqc_files
path salmon_files

output:
path "multiqc_report.html"

script:
"""
multiqc .
"""
}
