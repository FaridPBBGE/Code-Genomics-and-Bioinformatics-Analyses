RNA-Seq Pipeline with Nextflow & Docker
This repository contains a complete, automated pipeline for RNA-Seq analysis. It is built using Nextflow (DSL2) for workflow management and Docker for complete environment reproducibility.

🧬 Pipeline Workflow
The workflow integrates standard tools to process raw sequencing data into quantified expression levels:

Trimming: Adapter and quality trimming using Trimmomatic.

Quality Control: Post-trimming assessment with FastQC.

Alignment: Genomic mapping using STAR.

Quantification: Transcript-level abundance estimation using Salmon.

Reporting: Interactive summary of all steps using MultiQC.

🐳 Docker Hub Registry
All tools used in this pipeline are containerized in custom Docker images to ensure results are identical across any machine.

🚀 Getting Started
Prerequisites
Nextflow (>=22.10.0)
Docker

How to Run
Navigate to the RNAseq_test_sample directory.

Update the params in main.nf to point to your local data.

Execute the pipeline with the Docker profile:
