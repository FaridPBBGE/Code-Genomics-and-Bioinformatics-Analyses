# RNA-Seq Pipeline with Nextflow & Docker

This repository contains a professional-grade Nextflow pipeline for processing RNA-Seq data. It integrates quality control, read trimming, genomic alignment, and transcript quantification into a single, reproducible workflow.



## 🧬 Pipeline Workflow
The pipeline follows bioinformatics best practices by executing the following steps:
1. **Raw QC:** Initial quality assessment using **FastQC**.
2. **Trimming:** Removal of adapters and low-quality bases via **Trimmomatic**.
3. **Alignment:** Genomic mapping using **STAR**.
4. **Quantification:** Transcript-level abundance estimation using **Salmon**.
5. **Reporting:** Aggregating all logs and metrics into a single interactive **MultiQC** report.

## 🐳 Docker Environment
To ensure 100% reproducibility across different computing environments, every tool is containerized. All images are built on Ubuntu 22.04 and hosted on my [Docker Hub Registry](https://hub.docker.com/u/faridpbbge).

| Tool | Version | Docker Hub Link | Pull Command |
| :--- | :--- | :--- | :--- |
| **STAR** | 2.7.10b | [Link](https://hub.docker.com/r/faridpbbge/star) | `docker pull faridpbbge/star:2.7.10b` |
| **Salmon** | 1.10.0 | [Link](https://hub.docker.com/r/faridpbbge/salmon) | `docker pull faridpbbge/salmon:1.10.0` |
| **FastQC** | 0.12.1 | [Link](https://hub.docker.com/r/faridpbbge/fastqc) | `docker pull faridpbbge/fastqc:0.12.1` |
| **Trimmomatic** | 0.40 | [Link](https://hub.docker.com/r/faridpbbge/trimmomatic) | `docker pull faridpbbge/trimmomatic:0.40` |
| **MultiQC** | 1.21 | [Link](https://hub.docker.com/r/faridpbbge/multiqc) | `docker pull faridpbbge/multiqc:1.21` |

## 🚀 Usage

### 1. Prerequisites
* **Nextflow** (>=22.10.0)
* **Docker** installed and running.

### 2. Configuration
Input data paths are managed in the `params` section of `main.nf`. Ensure your paths for `reads`, `genomeDir`, and `gtf` are correct.

### 3. Execution
Run the pipeline using the pre-configured Docker profile:

```bash
nextflow run main.nf -profile docker
